package main

import (
	"bufio"
	"flag"
	"fmt"
	"io/ioutil"
	"os"
	"sort"
	"strconv"
	"strings"
)

// global variables
var Debug bool = true
//var ExpandRange int = 2000000 // 1MB FNA file
var ExpandRange int = 150
var SizeRangeRow int = 3

var Filter1_ConsensusLength int = 100
var Filter2_only_TR_shared = false

var Flag1_print_consensus = false
var Flag2_print_TR_with_flanking = true

func main() {
	// Step 1: parse parameters
	pFNA_ziheng := flag.String("FNA_ziheng", "", "ziheng's Filename of the FNA file")
	pDirTrf_ziheng := flag.String("DirTRF_ziheng", "", "ziheng's Dir of the TRF output (full of html files)")
	pOutDir_ziheng := flag.String("DirOut_ziheng", "", "ziheng's Filename of the output Dir")

	pFNA_harish := flag.String("FNA_harish", "", "harish's Filename of the FNA file")
	pDirTrf_harish := flag.String("DirTRF_harish", "", "harish's Dir of the TRF output (full of html files)")
	pOutDir_harish := flag.String("DirOut_harish", "", "harish's Filename of the output Dir")

	pFNA_candace := flag.String("FNA_candace", "", "candace's Filename of the FNA file")
	pDirTrf_candace := flag.String("DirTRF_candace", "", "candace's Dir of the TRF output (full of html files)")
	pOutDir_candace := flag.String("DirOut_candace", "", "candace's Filename of the output Dir")

	flag.Parse()

	vecName := []string{"candace", "harish", "ziheng"}
	vecFNA := []string{*pFNA_candace, *pFNA_harish, *pFNA_ziheng}
	vecTRF := []string{*pDirTrf_candace, *pDirTrf_harish, *pDirTrf_ziheng}
	vecOut := []string{*pOutDir_candace, *pOutDir_harish, *pOutDir_ziheng}

	vecMapScaffold := []map[string]*Scaffold{nil,nil,nil}
	// Step 2: obtain line range of each Scaffold, list html files
	// Step 3: read trf files and record info of TR
	for i, _ := range vecName {
		mapScaffold := getScaffolds(vecFNA[i], vecTRF[i])
		if mapScaffold == nil {
			fmt.Println("Fatal: Nil mapScaffold, exit")
			return
		}
		vecMapScaffold[i] = mapScaffold
	}


	// Step 4: process all mapScaffolds, compare entries between three dolphins, and only remain interesting entries
	processMapScaffold(vecMapScaffold, vecName)

	// Step 5: one thread to receive string and write into output file
	for i, _ := range vecName {
		generateOutput(vecMapScaffold[i], vecFNA[i], vecOut[i])
	}

}

type Scaffold struct {
	appearOrder int

	index   int    // example: 1002
	strName string // example: "scaffold_1002"
	head    string

	// the actual data is in [intHeadLine+1, intEndLine]
	intHeadLine int // example: 3411802
	intEndLine  int // example: 3412020
	intTotalBp  int // number of base pairs

	// FRS
	// Note that "XXX.s1.XXX.1.html" may be for any Scaffold
	vecFrsInfo      []*FrsInfo
	vecAllFrsEntry []*frsEntry  // vecAllFrsEntry is just the combination of all entries in vecFrsInfo
	vecGoodFrsEntry []*frsEntry
	vecSubSequence  []*SubSequence

}

type FrsInfo struct {
	ptrScaffold *Scaffold
	fileName    string // example: "GCA_011762515.1_mTurTru1.pat_genomic.fna.s62.2.5.7.80.10.50.2000.1.html"
	vecEntry    []*frsEntry
}

// SubSequence is for each subsequence in the compressed FNA. One scaffold may generate multiple subsequences
type SubSequence struct {
	ptrEntry []*frsEntry

	intID       int
	indiceRange indexRange
	strHead     string
	strBP       string
}

type frsEntry struct {
	ptrScaffold *Scaffold

	ptrSubSequence *SubSequence // map the frsEntry to the subsequence

	indiceRange indexRange
	expandRange indexRange
	//intPeriodSize          int
	//floatCopyNumber        float64
	//intConsensusSize       int
	//intPercentMatches      int
	//intPercentIndels       int
	//intScore               int
	//intA, intC, intG, intT int
	//floatEntropy           float64
	intPeriodSize          string
	floatCopyNumber        string
	intConsensusSize       string
	intPercentMatches      string
	intPercentIndels       string
	intScore               string
	intA, intC, intG, intT string
	floatEntropy           string

	strConsensus string
}

func generateOutput(m map[string]*Scaffold, FNA string, outDir string) {
	// Step 5.1: generate compressed FNA file

	// open FNA file
	fna, err := os.Open(FNA)
	if err != nil {
		fmt.Println(err)
		return
	}
	defer fna.Close()

	scannerFNA := bufio.NewScanner(fna)

	// create output file
	outFNA, err := os.Create(outDir + "/compressedFNA.fna")
	if err != nil {
		fmt.Println(err)
		panic(1)
	}
	writerFNA := bufio.NewWriter(outFNA)

	vecScaffold := []*Scaffold{}
	for _, scaffold := range m {
		vecScaffold = append(vecScaffold, scaffold)
	}
	sort.SliceStable(vecScaffold, func(i, j int) bool {
		return vecScaffold[i].appearOrder < vecScaffold[j].appearOrder
	})

	rowFNA := 0
	for i, scaffold := range vecScaffold {

		for _, subS := range scaffold.vecSubSequence {
			var strSubS string

			if Flag1_print_consensus == false {
				headReplace := "_" + strconv.Itoa(subS.intID) + ", "
				strSubSHead := strings.Replace(scaffold.head, ", ", headReplace, 1)
				strSubSHead += ". Note: the first BP in this subScaffold is BP NO." + strconv.Itoa(subS.indiceRange.begin) + " in the original scaffold"
				strSubS += strSubSHead + "\n"

				begin := subS.indiceRange.begin
				rowBegin := scaffold.intHeadLine + (begin-1)/80 + 1

				if rowBegin <= rowFNA {
					fmt.Println("Warning: rowBegin:", rowBegin, "\trowFNA:", rowFNA, "\tof scaffold:", scaffold.index)
				}
				found := false
				for scannerFNA.Scan() {
					rowFNA++
					if rowFNA == rowBegin {
						found = true
						strSubS += scannerFNA.Text() + "\n"
						break
					}
				}
				if found == false {
					fmt.Println("Warning: Didn't find the row of a subsequence, whose BP begin at:", begin, "\trowBegin at:", rowBegin, "\tof scaffold:", scaffold.index)
					panic(0)
				}

				var endRow int
				if Flag2_print_TR_with_flanking == false {
					endRow = rowBegin + SizeRangeRow
				} else {
					end := subS.indiceRange.end
					endRow = scaffold.intHeadLine + end / 80
				}

				for scannerFNA.Scan() {
					rowFNA++
					if rowFNA < endRow {
						strSubS += scannerFNA.Text() + "\n"
					} else if rowFNA == endRow {
						strSubS += scannerFNA.Text() + "\n"
						break
					}
				}
			} else {
				// just write the subscaffold using the consensus
				headReplace := "_" + strconv.Itoa(subS.intID) + ", "
				strSubSHead := strings.Replace(scaffold.head, ", ", headReplace, 1)
				strSubSHead += ". Note: the subscaffold is just one repeat of the consensus"
				strSubS += strSubSHead + "\n"

				if len(subS.ptrEntry) != 1 {
					fmt.Printf("A subscaffold has %d entries when Flag1 is on!\n", len(subS.ptrEntry))
					panic(1)
				}
				consensus := subS.ptrEntry[0].strConsensus
				x := (len(consensus) - 1) / 80
				j := 0
				for j < x {
					strSubS += consensus[j*80: j*80+80] + "\n"
					j++
				}
				strSubS += consensus[j*80:] + "\n"
			}



			_, err := writerFNA.WriteString(strSubS)
			if err != nil {
				fmt.Println("Got error while writing to a file. Err: %s", err.Error())
				panic(1)
			}
		}

		if i%20 == 0 {
			writerFNA.Flush()
		}
	}

	// Step 5.2: generate TR table file
	outTR, err := os.Create(outDir + "/TR_Table.txt")
	if err != nil {
		fmt.Println(err)
		panic(1)
	}
	writerTR := bufio.NewWriter(outTR)
	// format: each TR has one line. Each line is like below
	//scaffoldNum,subScaffoldNum,,offsetToHeadOfSubScaffold,indicesBegin,indicesEnd,periodSize,copyNumber,consensusSize,percentMatches,percentIndels,score,A,C,G,T,entropy
	writerTR.WriteString("scaffoldNum,subScaffoldNum,offsetToHeadOfSubScaffold,indicesBegin,indicesEnd,periodSize,copyNumber,consensusSize,percentMatches,percentIndels,score,A,C,G,T,entropy,consensusPattern" + "\n")
	for i, scaffold := range vecScaffold {
		for _, entry := range scaffold.vecGoodFrsEntry {
			str := strconv.Itoa(scaffold.index) + ","
			str += strconv.Itoa(entry.ptrSubSequence.intID) + ","
			str += strconv.Itoa(entry.indiceRange.begin-entry.ptrSubSequence.indiceRange.begin) + ","
			str += strconv.Itoa(entry.indiceRange.begin) + ","
			str += strconv.Itoa(entry.indiceRange.end) + ","
			str += entry.intPeriodSize + "," + entry.floatCopyNumber + "," + entry.intConsensusSize + "," + entry.intPercentMatches + "," + entry.intPercentIndels + ","
			str += entry.intScore + "," + entry.intA + "," + entry.intC + "," + entry.intG + "," + entry.intT + "," + entry.floatEntropy + ","
			str += entry.strConsensus
			str += "\n"
			if len(str) < 10 {
				fmt.Println("str less than 10 in scaffold", scaffold.index, " entry:", entry.indiceRange.begin, " subsequence:", &entry.ptrSubSequence.indiceRange.begin)
			}
			writerTR.WriteString(str)
		}
		if i%20 == 0 {
			writerTR.Flush()
		}
	}
	writerTR.Flush()
	return
}

func processMapScaffold(vecMap []map[string]*Scaffold, vecName []string) {
	mapBadConsensus := make(map[string]struct{})
	vecMapLongEntry := []map[string]*frsEntry{make(map[string]*frsEntry),make(map[string]*frsEntry),make(map[string]*frsEntry)}
	// decide what entries should be remained
	countRepeatConsensusInSameDolphin := 0
	for i, mapScaffold := range vecMap {
		for _, scaffold := range mapScaffold {
			for _, frsInfo := range scaffold.vecFrsInfo {
				for _, entry := range frsInfo.vecEntry {
					scaffold.vecAllFrsEntry = append(scaffold.vecAllFrsEntry, entry)
					if len(entry.strConsensus) >= Filter1_ConsensusLength {

						if _, exist := vecMapLongEntry[i][entry.strConsensus]; exist {
							// this consensus repeats multiple times even in one dolphin
							countRepeatConsensusInSameDolphin ++
							mapBadConsensus[entry.strConsensus] = struct{}{}
						}

						vecMapLongEntry[i][entry.strConsensus] = entry
					}
				}
			}
		}
	}

	// display some statistics of entries
	fmt.Println("========Statistics of length of consensus")
	countLongerThanFilter1 := 0
	longest := 0
	for i, mapScaffold := range vecMap {
		name := vecName[i]
		var count1, count10, count100, count1000, countRest int
		for _, scaffold := range mapScaffold {
			for _, entry := range scaffold.vecAllFrsEntry {
				lenConsensus := len(entry.strConsensus)
				if lenConsensus == 1 {
					count1++
				} else if lenConsensus < 10 {
					count10++
				} else if lenConsensus < 100 {
					count100++
				} else if lenConsensus < 1000 {
					count1000++
				} else  {
					countRest++
				}
				if lenConsensus > Filter1_ConsensusLength {
					countLongerThanFilter1++
				}
				if lenConsensus > longest {
					longest = lenConsensus
				}
			}
		}
		fmt.Println("===Dolphin of", name)
		fmt.Printf("\tLen 1:%d\tLen 1 to 10:%d\tLen 10 to 100:%d\tLen 100 to 1000:%d\tLen 1000+:%d\n", count1, count10, count100, count1000, countRest)
	}
	fmt.Printf("==========\n\tSum up all dolphins, there are %d TR whose consensus is longer than %d " +
		"(%d TRs longer than this), but is not unique even in the same dolphin\n", countRepeatConsensusInSameDolphin, Filter1_ConsensusLength, countLongerThanFilter1)
	fmt.Printf("==========\n\tThe greatest consensus length is %d\n", longest)

	// find entries shared by three dolphins
	mapGoodEntry := make(map[*frsEntry]struct{})
	for i, mapLongEntry := range vecMapLongEntry {
		for key, entry := range mapLongEntry {
			if Filter2_only_TR_shared {
				countFound := 0
				for j, mapOther := range vecMapLongEntry {
					if i == j {
						continue
					}
					_, exist := mapOther[key]
					if exist {
						countFound ++
					}
				}
				if countFound == 2 {
					// good entry
					if _, exist2 := mapBadConsensus[entry.strConsensus]; exist2 {
						// this consensus repeats multiple times even in one dolphin
						continue
					}
					mapGoodEntry[entry] = struct{}{}
				}
			} else {
				if _, exist2 := mapBadConsensus[entry.strConsensus]; exist2 {
					// this consensus repeats multiple times even in one dolphin
					continue
				}
				mapGoodEntry[entry] = struct{}{}
			}
		}
	}


	countEventualTR := 0
	for _, mapScaffold := range vecMap {
		for _, scaffold := range mapScaffold {
			for _, entry := range scaffold.vecAllFrsEntry {
				if _, exist := mapGoodEntry[entry]; exist {
					if Flag1_print_consensus == true {
						// see if this entry's one consensus overlaps with an existing entry
						found_overlap := false
						for _, other_entry := range scaffold.vecGoodFrsEntry {
							if other_entry.indiceRange.begin <= entry.indiceRange.begin && entry.indiceRange.begin <= other_entry.indiceRange.begin + len(other_entry.strConsensus) + 80 {
								found_overlap = true
								break
							}
						}
						if found_overlap {
							continue
						}
					}
					if Flag2_print_TR_with_flanking == true {
						// see if this entry's extended range overlaps with an existing entry
						found_overlap := false
						for _, other_entry := range scaffold.vecGoodFrsEntry {
							if isOverlap(other_entry.expandRange, entry.expandRange) {
								found_overlap = true
								break
							}
						}
						if found_overlap {
							continue
						}
					}
					countEventualTR ++
					scaffold.vecGoodFrsEntry = append(scaffold.vecGoodFrsEntry, entry)
				}
			}
			sort.SliceStable(scaffold.vecGoodFrsEntry, func(i, j int) bool {
				return scaffold.vecGoodFrsEntry[i].expandRange.end-scaffold.vecGoodFrsEntry[i].expandRange.begin > scaffold.vecGoodFrsEntry[j].expandRange.end-scaffold.vecGoodFrsEntry[j].expandRange.begin
			})


			if Flag1_print_consensus == false {
				if Flag2_print_TR_with_flanking == false {
					for _, entry := range scaffold.vecGoodFrsEntry {
						eBegin := entry.expandRange.begin
						eEnd := entry.expandRange.end
						// for each subsequence, see if there is one that overlaps with this entry
						// if there is, update the range of the subsequence
						// if not, create one
						for _, subS := range scaffold.vecSubSequence {
							sBegin := subS.indiceRange.begin
							sEnd := subS.indiceRange.end

							if isIn(eBegin, sBegin, sEnd) || isIn(eEnd, sBegin, sEnd) || isIn(sBegin, eBegin, eEnd) || isIn(sEnd, eBegin, eEnd) {
								var newBegin, newEnd int
								if eBegin > sBegin {
									newBegin = sBegin
								} else {
									newBegin = eBegin
								}
								if eEnd < sEnd {
									newEnd = sEnd
								} else {
									newEnd = eEnd
								}
								subS.indiceRange = indexRange{
									begin: newBegin,
									end:   newEnd,
								}

								entry.ptrSubSequence = subS
								subS.ptrEntry = append(subS.ptrEntry, entry)
								break
							}
						}

						if entry.ptrSubSequence == nil {
							newSubSequence := &SubSequence{
								intID:       -1,
								indiceRange: indexRange{begin: eBegin, end: eEnd},
							}
							scaffold.vecSubSequence = append(scaffold.vecSubSequence, newSubSequence)

							entry.ptrSubSequence = newSubSequence
							newSubSequence.ptrEntry = append(newSubSequence.ptrEntry, entry)
						}
					}
				} else {
					// Flag2: make the subscaffold include the TR and its left and right flanking sequence
					for _, entry := range scaffold.vecGoodFrsEntry {
						newSubSequence := &SubSequence{
							intID:       -1,
							indiceRange: indexRange{begin: entry.expandRange.begin, end: entry.expandRange.end},
						}
						scaffold.vecSubSequence = append(scaffold.vecSubSequence, newSubSequence)

						entry.ptrSubSequence = newSubSequence
						newSubSequence.ptrEntry = append(newSubSequence.ptrEntry, entry)
					}
				}

			} else {
				// Flag1: just make the subscaffold the same as one consensus of the entry
				for _, entry := range scaffold.vecGoodFrsEntry {
					newSubSequence := &SubSequence{
						intID:       -1,
						indiceRange: indexRange{begin: entry.indiceRange.begin, end: entry.indiceRange.begin + len(entry.strConsensus)},
					}
					scaffold.vecSubSequence = append(scaffold.vecSubSequence, newSubSequence)

					entry.ptrSubSequence = newSubSequence
					newSubSequence.ptrEntry = append(newSubSequence.ptrEntry, entry)
				}
			}

			intSubSID := 1
			for _, subS := range scaffold.vecSubSequence {
				subS.intID = intSubSID
				intSubSID++
			}

			// code to make the sub-sequences not overlap
			if Flag1_print_consensus == false {
				mapSubSequence := make(map[int]*SubSequence)
				for _, subS := range scaffold.vecSubSequence {
					mapSubSequence[subS.intID] = subS
				}

				detected := true
			detect_loop:
				for detected {
					detected = false
					for _, subS := range mapSubSequence {
						for _, other := range mapSubSequence {
							if other == subS {
								continue
							}

							if other.indiceRange.begin >= subS.indiceRange.begin && other.indiceRange.begin <= subS.indiceRange.begin + SizeRangeRow * 80 + 80 {
								detected = true
								for _, entry := range other.ptrEntry {
									subS.ptrEntry = append(subS.ptrEntry, entry)
									entry.ptrSubSequence = subS
								}
								delete(mapSubSequence, other.intID)
								continue detect_loop
							}
						}
					}
				}
				scaffold.vecSubSequence = nil
				for _, subS := range mapSubSequence {
					scaffold.vecSubSequence = append(scaffold.vecSubSequence, subS)
				}

				for _, subS := range scaffold.vecSubSequence {
					for _, other := range scaffold.vecSubSequence {
						if other == subS {
							continue
						}

						if other.indiceRange.begin >= subS.indiceRange.begin && other.indiceRange.begin <= subS.indiceRange.begin + SizeRangeRow * 80 + 80 {
							fmt.Println("Fatal: Still has duplicate subsequence!")
							panic(2)
						}
					}
				}
			} else {

			}


			indexNew := 1
			sort.SliceStable(scaffold.vecSubSequence, func(i, j int) bool {
				return scaffold.vecSubSequence[i].indiceRange.begin < scaffold.vecSubSequence[j].indiceRange.begin
			})
			for _, subS := range scaffold.vecSubSequence {
				subS.intID = indexNew
				indexNew++
			}

		}
	}

	fmt.Printf("============\n\tFound %d TR that is: (1) longer than %d; (2) its consensus doesn't appear " +
		"in other TRs of the same dolphin; (3) all three dolphins have the exact same consensus\n", countEventualTR/3, Filter1_ConsensusLength)


	return
}

func getScaffolds(fnaFileName string, dirTrf string) (mapScaffold map[string]*Scaffold) {
	defer func() {
		print()
	}()
	appearOrderCount := 0
	mapScaffold = make(map[string]*Scaffold)

	// Step 1.1: parse FNA file
	fna, err := os.Open(fnaFileName)
	if err != nil {
		fmt.Println(err)
		return nil
	}
	defer fna.Close()

	scannerFNA := bufio.NewScanner(fna)
	rowFNA := 0

	numRowHead := 0

	// Find lines like "scaffold_X, ", and record Scaffold
	var lastScaffold *Scaffold
	var strRow, strLastRow string
	for scannerFNA.Scan() {
		rowFNA++
		strRow = scannerFNA.Text()
		if strings.HasPrefix(strRow, ">") {
			numRowHead ++
			iSca := strings.Index(strRow, "scaffold_") // ziheng and harish's format
			iScaOffSet := 9
			if iSca == -1 {
				iSca = strings.Index(strRow, "scaffold") // candace's format
				iScaOffSet = 8
			}

			if iSca > -1 {
				iComma := strings.Index(strRow, ", ")
				if iComma == -1 || iComma < iSca+iScaOffSet {
					if Debug {
						fmt.Println("Warning: A line in FNA has strange iComma\n\t" + strRow)
						continue
					}
				}
				strIndex := strRow[iSca+iScaOffSet : iComma]
				intIndex, err := strconv.Atoi(strIndex)
				if err != nil && Debug {
					fmt.Println("Warning: A line in FNA has strange format when converting index\n\t" + strRow)
					continue
				}

				// found a new scaffold
				if lastScaffold != nil {
					lastScaffold.intEndLine = rowFNA - 1
					lastScaffold.intTotalBp = (lastScaffold.intEndLine-lastScaffold.intHeadLine-1)*80 + len(strLastRow)
				}

				strName := "scaffold_" + strIndex
				newScaffold := &Scaffold{
					appearOrder:     appearOrderCount,
					index:           intIndex,
					strName:         strName,
					head:            strRow,
					intHeadLine:     rowFNA,
					intEndLine:      0,
					vecFrsInfo:      nil,
				}
				appearOrderCount++
				if lastScaffold != nil {
					lastScaffold.intEndLine = rowFNA - 1
				}
				mapScaffold[strName] = newScaffold

				lastScaffold = newScaffold
			} else {
				if Debug {
					fmt.Println("Warning: A line in FNA has > but no scaffold_\n\t" + strRow)
					continue
				}
			}
		}
		strLastRow = strRow
	}

	// when the scan ends, write intEndLine to the lastScaffold
	lastScaffold.intEndLine = rowFNA
	lastScaffold.intTotalBp = (lastScaffold.intEndLine-lastScaffold.intHeadLine-1)*80 + len(strLastRow)

	if err := scannerFNA.Err(); err != nil {
		fmt.Println(err)
		return nil
	}

	//fmt.Println("Number of lines starting with >:", numRowHead, "Number of scaffolds:", len(mapScaffold))

	// Step 1.2: list dirTrf
	vecTrfFile, err := ioutil.ReadDir(dirTrf)
	if err != nil {
		fmt.Println(err)
		return nil
	}

	for _, trfFileInfo := range vecTrfFile {

		// ignore files ending with ".txt.html" or ".summary.html", or files not ending with ".html"
		if strings.HasSuffix(trfFileInfo.Name(), ".txt.html") || strings.HasSuffix(trfFileInfo.Name(), ".summary.html") || !strings.HasSuffix(trfFileInfo.Name(), ".html") {
			continue
		}

		var newFrsInfo *FrsInfo

		trfFile, err := os.Open(dirTrf + "/" + trfFileInfo.Name())
		if err != nil {
			fmt.Println(err)
			return nil
		}

		mapEntry := make(map[indexRange]*frsEntry)

		scannerTrf := bufio.NewScanner(trfFile)
		rowTrf := 0

		for scannerTrf.Scan() {
			rowTrf++
			strRowTrf := scannerTrf.Text()

			if strings.HasPrefix(strRowTrf, "Sequence:") {
				iSca := strings.Index(strRowTrf, "scaffold_") // ziheng and harish's format
				iScaOffSet := 9
				if iSca == -1 {
					iSca = strings.Index(strRowTrf, "scaffold") // candace's format
					iScaOffSet = 8
				}
				if iSca > -1 {
					iComma := strings.Index(strRowTrf, ", ")
					if iComma == -1 || iComma < iSca+iScaOffSet {
						if Debug {
							fmt.Println("Warning: A line in file in TRF Dir has strange iComma\n\t" + strRowTrf)
							fmt.Println("\tFile location:", trfFileInfo.Name())
							panic(3)
						}
					}
					strIndex := strRowTrf[iSca+iScaOffSet : iComma]
					_, err := strconv.Atoi(strIndex)
					if err != nil && Debug {
						fmt.Println("Warning: A line in file in TRF Dir has strange format when converting index\n\t" + strRowTrf)
						fmt.Println("\tFile location:", trfFileInfo.Name())
						panic(3)
					}
					strName := "scaffold_" + strIndex
					scaffold := mapScaffold[strName]
					newFrsInfo = &FrsInfo{
						fileName:    dirTrf + "/" + trfFileInfo.Name(),
						vecEntry:    nil,
						ptrScaffold: scaffold,
					}
					scaffold.vecFrsInfo = append(scaffold.vecFrsInfo, newFrsInfo)
				} else {
					if Debug {
						fmt.Println("Warning: A line in file in TRF Dir has Sequence: but no scaffold_\n\t" + strRowTrf)
						fmt.Println("\tFile location:", trfFileInfo.Name())
						panic(3)
					}
				}
			}

			if strings.HasPrefix(strRowTrf, "<TR><TD><CENTER><A HREF=") {
				if newFrsInfo == nil {
					fmt.Println("Warning: A line in file in TRF Dir has table head but newFrsInfo is nil\n\t" + strRowTrf)
					fmt.Println("\tFile location:", trfFileInfo.Name())
					panic(3)
				}

				iTxt := strings.Index(strRowTrf, "\">")
				iDash := strings.LastIndex(strRowTrf, "--")
				iBeforePeriod := strings.Index(strRowTrf, "</A></CENTER>")

				if iTxt == -1 || iDash == -1 || iBeforePeriod == -1 || iTxt > iDash || iDash > iBeforePeriod {
					fmt.Println("Warning: A line in file in TRF Dir has table head but the format is strange for range\n\t" + strRowTrf)
					fmt.Println("\tFile location:", trfFileInfo.Name())
					panic(3)
				}
				strBegin := strRowTrf[iTxt+2 : iDash]
				strEnd := strRowTrf[iDash+2 : iBeforePeriod]
				intBegin, err1 := strconv.Atoi(strBegin)
				intEnd, err2 := strconv.Atoi(strEnd)
				if err1 != nil {
					fmt.Println(err1)
					panic(3)
				} else if err2 != nil {
					fmt.Println(err2)
					panic(3)
				}
				r := indexRange{
					begin: intBegin,
					end:   intEnd,
				}
				enpandR := expandRange(r, newFrsInfo.ptrScaffold)
				entry := &frsEntry{}
				entry.indiceRange = r
				entry.expandRange = enpandR

				str := strRowTrf[iBeforePeriod+13:]
				for i := 0; i < 11; i++ {
					iBegin := strings.Index(str, "</TD><TD><CENTER>")
					iEnd := strings.Index(str, "</CENTER>")
					if iBegin == -1 || iEnd == -1 || iEnd <= iBegin {
						fmt.Println("Warning: A line in file in TRF Dir has table head but the format is strange later than range\n\t" + strRowTrf)
						fmt.Println("\tFile location:", trfFileInfo.Name())
						panic(3)
					}
					content := str[iBegin+17 : iEnd]

					/*
					//vecErr := [11]error{}
					//if i == 0 {
					//	entry.intPeriodSize, vecErr[i] = strconv.Atoi(content)
					//} else if i == 1 {
					//	entry.floatCopyNumber, vecErr[i] = strconv.ParseFloat(content, 32)
					//} else if i == 2 {
					//	entry.intConsensusSize, vecErr[i] = strconv.Atoi(content)
					//} else if i == 3 {
					//	entry.intPercentMatches, vecErr[i] = strconv.Atoi(content)
					//} else if i == 4 {
					//	entry.intPercentIndels, vecErr[i] = strconv.Atoi(content)
					//} else if i == 5 {
					//	entry.intScore, vecErr[i] = strconv.Atoi(content)
					//} else if i == 6 {
					//	entry.intA, vecErr[i] = strconv.Atoi(content)
					//} else if i == 7 {
					//	entry.intC, vecErr[i] = strconv.Atoi(content)
					//} else if i == 8 {
					//	entry.intG, vecErr[i] = strconv.Atoi(content)
					//} else if i == 9 {
					//	entry.intT, vecErr[i] = strconv.Atoi(content)
					//} else if i == 10 {
					//	entry.floatEntropy, vecErr[i] = strconv.ParseFloat(content, 32)
					//}
					//for _, err := range vecErr {
					//	if err != nil {
					//		fmt.Println(err)
					//		fmt.Println("\tError at line", strRowTrf, "\n\tError in File:", trfFileInfo.Name())
					//	}
					//}
					*/

					if i == 0 {
						entry.intPeriodSize = content
					} else if i == 1 {
						entry.floatCopyNumber = content
					} else if i == 2 {
						entry.intConsensusSize = content
					} else if i == 3 {
						entry.intPercentMatches = content
					} else if i == 4 {
						entry.intPercentIndels = content
					} else if i == 5 {
						entry.intScore = content
					} else if i == 6 {
						entry.intA = content
					} else if i == 7 {
						entry.intC = content
					} else if i == 8 {
						entry.intG = content
					} else if i == 9 {
						entry.intT = content
					} else if i == 10 {
						entry.floatEntropy = content
					}

					str = str[iEnd+9:]
				}

				mapEntry[entry.indiceRange] = entry
				newFrsInfo.vecEntry = append(newFrsInfo.vecEntry, entry)
			}
		}

		if err := scannerTrf.Err(); err != nil {
			fmt.Println(err)
			trfFile.Close()
			return nil
		}

		trfFile.Close()

		txtFileName := strings.Replace(trfFileInfo.Name(), ".html", ".txt.html", 1)
		txtFile, err := os.Open(dirTrf + "/" + txtFileName)
		if err != nil {
			fmt.Println(err)
			return nil
		}

		scannerTxt := bufio.NewScanner(txtFile)
		rowTxt := 0

		txtCurrentRange := indexRange{
			begin: -1,
			end:   -1,
		}
		consensusRow := -1
		for scannerTxt.Scan() {
			rowTxt++
			strRowTxt := scannerTxt.Text()

			// if find the line of indices
			iHeadRange := strings.Index(strRowTxt, "<A NAME=\"")
			if iHeadRange > -1 {
				// check if we didn't find the consensus the last time
				if consensusRow != -1 {
					fmt.Println("Warning: A line in .txt.html file has one missing consensus\n\t" + strRowTxt, "line NO.", rowTxt)
					fmt.Println("\tFile location:", txtFileName)
					panic(3)
				}

				iDash := strings.Index(strRowTxt, "--")
				iComma := strings.Index(strRowTxt, ",")
				if iHeadRange == -1 || iDash == -1 || iComma == -1 || iHeadRange > iDash || iDash > iComma {
					fmt.Println("Warning: A line in .txt.html file has range head but the format is strange\n\t" + strRowTxt, "line NO.", rowTxt)
					fmt.Println("\tFile location:", txtFileName)
					panic(3)
				}
				strBegin := strRowTxt[iHeadRange+9 : iDash]
				strEnd := strRowTxt[iDash+2 : iComma]
				intBegin, err1 := strconv.Atoi(strBegin)
				intEnd, err2 := strconv.Atoi(strEnd)
				if err1 != nil {
					fmt.Println(err1)
					panic(3)
				} else if err2 != nil {
					fmt.Println(err2)
					panic(3)
				}
				txtCurrentRange = indexRange{
					begin: intBegin,
					end:   intEnd,
				}
				continue
			}

			// if find the row of consensus
			iConsensus := strings.Index(strRowTxt, "Consensus pattern")
			if iConsensus > -1 {
				consensusRow = iConsensus
				continue
			}

			// if found consensus previously, see if range exists
			if consensusRow != -1 {
				if strRowTxt == "" {
					// reset both two variables
					txtCurrentRange.begin = -1; txtCurrentRange.end = -1
					consensusRow = -1
					continue
				}
				if txtCurrentRange.begin == -1 || txtCurrentRange.end == -1 {
					fmt.Println("Warning: A line in .txt.html file has a consensus without range info\n\t" + strRowTxt, "line NO.", rowTxt)
					fmt.Println("\tFile location:", txtFileName)
					panic(3)
				}
				mapEntry[txtCurrentRange].strConsensus += strRowTxt
			}
		}

		for _, entry := range mapEntry {
			if len(entry.strConsensus) == 0 {
				fmt.Println("Warning: A .txt.html file has an entry with no consensus\n\t", entry.indiceRange)
				fmt.Println("\tFile location:", txtFileName)
				panic(3)
			}
		}

		txtFile.Close()
	}

	return
}

type indexRange struct {
	// TR is in [begin,end]
	begin, end int
}

func isIn(a, b, c int) bool {
	if a >= b && a <= c {
		return true
	} else {
		return false
	}
}

func expandRange(r indexRange, scaffold *Scaffold) indexRange {
	// expand range to include ExpandRange Left flanking sequence and ExpandRange Right flanking sequence
	// and then make the begin to be 80N+1, and end to be 80M
	// and then make sure they are in range

	begin := r.begin
	end := r.end
	begin -= ExpandRange
	begin = ((begin - 1) / 80) * 80 + 1
	end += ExpandRange
	end = (((end - 1) / 80) + 1) * 80

	if begin <= 0 {
		begin = 1
	}
	if end > scaffold.intTotalBp {
		end = scaffold.intTotalBp
	}
	return indexRange{
		begin: begin,
		end:   end,
	}
}

func isOverlap(r1, r2 indexRange) bool {
	if r1.begin - 80 <= r2.begin && r2.begin <= r1.end + 80 {
		return true
	}
	if r1.begin - 80 <= r2.end && r2.end <= r1.end + 80 {
		return true
	}
	if r2.begin - 80 <= r1.begin && r1.begin <= r2.end + 80 {
		return true
	}
	if r2.begin - 80 <= r1.end && r1.end <= r2.end + 80 {
		return true
	}
	return false
}