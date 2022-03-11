package main

import (
	"bufio"
	"fmt"
	"os"
	"sort"
	"strconv"
	"strings"
)

type TR struct {
	index                                        int
	chr_name                                     string
	chr_ID                                       int
	pos_begin                                    float64
	pos_end                                      float64
	cp_harish                                    float64
	cp_candace                                   float64
	cp_ziheng                                    float64
	cp_str_harish, cp_str_candace, cp_str_ziheng string
	seq_harish, seq_candace, seq_ziheng          string
}

func main() {
	folderStr := "/Users/ziheng/Desktop/courses/CSE280A/project/repo/results"
	fileStr := folderStr + "/experiment_6_raw.csv"

	// open FNA file
	file, err := os.Open(fileStr)
	if err != nil {
		fmt.Println(err)
		return
	}
	defer file.Close()

	scanner := bufio.NewScanner(file)

	vecTR := []*TR{}

	rowNum := 0
	indexTR := 0
	for scanner.Scan() {
		rowNum++
		if rowNum == 1 {
			continue
		}
		newTR := &TR{}
		rowStr := scanner.Text()
		rowSlice := strings.Split(rowStr, ",")

		for i, column := range rowSlice {
			switch i {
			case 0:
				newTR.index = indexTR
				indexTR++
			case 1:
				newTR.chr_name = column
			case 2:
				newTR.chr_ID, _ = strconv.Atoi(column)
			case 3:
				newTR.pos_begin, _ = strconv.ParseFloat(column, 64)
			case 4:
				newTR.pos_end, _ = strconv.ParseFloat(column, 64)
			case 5:
				newTR.cp_harish, _ = strconv.ParseFloat(column, 64)
				newTR.cp_str_harish = column
			case 6:
				newTR.cp_candace, _ = strconv.ParseFloat(column, 64)
				newTR.cp_str_candace = column
			case 7:
				newTR.cp_ziheng, _ = strconv.ParseFloat(column, 64)
				newTR.cp_str_ziheng = column
			case 8:
				newTR.seq_harish = column
			case 9:
				newTR.seq_candace = column
			case 10:
				newTR.seq_ziheng = column
			}
		}

		vecTR = append(vecTR, newTR)
	}

	// 1. Draw Venn diagrams (Ziheng)
	var countOnlyZ, countOnlyH, countOnlyC, countZ_H, countH_C, countZ_C, countALL, countNone int

	var Q2_onlyTwo_same_copy int
	var Q5_same_three, Q5_same_two int

	for _, tr := range vecTR {
		if tr.cp_ziheng != 0 {
			if tr.cp_harish != 0 {
				if tr.cp_candace != 0 {
					countALL++
					if tr.seq_ziheng == tr.seq_harish && tr.seq_harish == tr.seq_candace {
						Q5_same_three++
					} else if tr.seq_ziheng == tr.seq_harish || tr.seq_harish == tr.seq_candace || tr.seq_candace == tr.seq_ziheng {
						Q5_same_two++
					}
				} else {
					countZ_H++
					if tr.cp_ziheng == tr.cp_harish {
						Q2_onlyTwo_same_copy++
					}
					if tr.seq_ziheng == tr.seq_harish {
						Q5_same_two++
					}
				}
			} else {
				if tr.cp_candace != 0 {
					countZ_C++
					if tr.cp_ziheng == tr.cp_candace {
						Q2_onlyTwo_same_copy++
					}
					if tr.seq_ziheng == tr.seq_candace {
						Q5_same_two++
					}
				} else {
					countOnlyZ++
				}
			}
		} else {
			if tr.cp_harish != 0 {
				if tr.cp_candace != 0 {
					countH_C++
					if tr.cp_candace == tr.cp_harish {
						Q2_onlyTwo_same_copy++
					}
					if tr.seq_candace == tr.seq_harish {
						Q5_same_two++
					}
				} else {
					countOnlyH++
				}
			} else {
				if tr.cp_candace != 0 {
					countOnlyC++
				} else {
					countNone++
				}
			}
		}
	}

	fmt.Printf("========Venn Diagram\n\tTotal TR Number:%d\n", len(vecTR))
	fmt.Printf("\tZiheng Only:%d\tZiheng_Harish:%d\tZiheng_Candace:%d\n", countOnlyZ, countZ_H, countZ_C)
	fmt.Printf("\tHarish Only:%d\tHarish_Candace:%d\n", countOnlyH, countH_C)
	fmt.Printf("\tCandace Only:%d\tAll Three:%d\t None:%d\n\n", countOnlyC, countALL, countNone)

	//How many of the consensus sequences have the same copy number in all three dolphins? (Ziheng)
	Q2_count_Same_Copy_All, Q2_count_Same_Copy_Two := 0, 0
	for _, tr := range vecTR {
		if tr.cp_ziheng == tr.cp_harish && tr.cp_ziheng == tr.cp_candace && tr.cp_ziheng != 0 {
			Q2_count_Same_Copy_All++
		} else if (tr.cp_harish == tr.cp_candace && tr.cp_harish != 0) || (tr.cp_harish == tr.cp_ziheng && tr.cp_harish != 0) || (tr.cp_ziheng == tr.cp_candace && tr.cp_ziheng != 0) {
			Q2_count_Same_Copy_Two++
		}
	}
	fmt.Printf("========How many TR have the same copy number\n")
	fmt.Printf("\tBetween three dolphins:%d\n\tBetween two dolphins:%d\n"+
		"\tOnly two have this TR and the copy num is same:%d\n", Q2_count_Same_Copy_All, Q2_count_Same_Copy_Two, Q2_onlyTwo_same_copy)

	//Which chromosomes have the most long (>100bp) TRs?  (Ziheng)
	type chr_count struct {
		chr_ID int
		count  int
	}
	mapChrID2countTR := make(map[int]int)
	for _, tr := range vecTR {
		if countTR, exist := mapChrID2countTR[tr.chr_ID]; exist {
			mapChrID2countTR[tr.chr_ID] = countTR + 1
		} else {
			mapChrID2countTR[tr.chr_ID] = 1
		}
	}
	vecChrCount := []chr_count{}
	for chr_ID, count := range mapChrID2countTR {
		vecChrCount = append(vecChrCount, chr_count{
			chr_ID: chr_ID,
			count:  count,
		})
	}
	sort.SliceStable(vecChrCount, func(i, j int) bool {
		return vecChrCount[i].count < vecChrCount[j].count
	})
	fmt.Printf("========Which chromosomes have most TR or least TR\n\t===Least\n")
	for i := 0; i < 3; i++ {
		fmt.Printf("\tNO.%d Chromosome %d has %d Tandeom Repeats\n", i+1, vecChrCount[i].chr_ID, vecChrCount[i].count)
	}
	fmt.Println("\t===Most")
	for i := len(vecChrCount) - 3; i < len(vecChrCount); i++ {
		fmt.Printf("\tNO.%d Chromosome %d has %d Tandeom Repeats\n", i+1, vecChrCount[i].chr_ID, vecChrCount[i].count)
	}

	sort.SliceStable(vecChrCount, func(i, j int) bool {
		return vecChrCount[i].chr_ID < vecChrCount[j].chr_ID
	})
	//for _, c := range vecChrCount {
	//	fmt.Printf("%d ", c.count)
	//}
	//fmt.Println()
	//for _, c := range vecChrCount {
	//	fmt.Printf("%d ", c.chr_ID)
	//}

	// Q4:How many have float copy number
	Q4_float_copy := 0
	Q4_int_copy := 0
	for _, tr := range vecTR {
		if tr.cp_ziheng != 0 {
			if strings.Contains(tr.cp_str_ziheng, ".0") {
				Q4_int_copy++
			} else {
				Q4_float_copy++
			}
		}
		if tr.cp_harish != 0 {
			if strings.Contains(tr.cp_str_harish, ".0") {
				Q4_int_copy++
			} else {
				Q4_float_copy++
			}
		}
		if tr.cp_candace != 0 {
			if strings.Contains(tr.cp_str_candace, ".0") {
				Q4_int_copy++
			} else {
				Q4_float_copy++
			}
		}
	}
	fmt.Printf("=======Among all copy numbers that are not 0, %d are integer, and %d are float number\n", Q4_int_copy, Q4_float_copy)

	// Q5: How many consensus are the same

	fmt.Printf("=======Among all existing consensus of the same TR, %d are all three exactly the same, and %d have two exactly the same\n", Q5_same_three, Q5_same_two)

	//Matrix of consensus sequences with different copy numbers? (Ziheng)

	// Head: #CHROM POS NA00001

	outFileStr := folderStr + "/snp.txt"
	outFile, err := os.Create(outFileStr)
	if err != nil {
		fmt.Println(err)
		panic(1)
	}
	defer outFile.Close()

	writer := bufio.NewWriter(outFile)
	writer.WriteString("#CHROM\tPOS\t\tDolphin_H\tDolphin_Z\tDolphin_C\n")
	sort.SliceStable(vecTR, func(i, j int) bool {
		if vecTR[i].chr_ID != vecTR[j].chr_ID {
			return vecTR[i].chr_ID < vecTR[j].chr_ID
		} else {
			return vecTR[i].pos_begin < vecTR[j].pos_begin
		}
	})
	for _, tr := range vecTR {

		str := strconv.Itoa(tr.chr_ID+1) + "\t"
		str += strconv.Itoa(int(tr.pos_begin)) + "\t"
		if tr.pos_begin < 10000000 {
			str += "\t"
		}
		str += tr.cp_str_harish + "\t\t"
		str += tr.cp_str_ziheng + "\t\t"
		str += tr.cp_str_candace
		str += "\n"

		str = strings.ReplaceAll(str, "0.0", "-")
		str = strings.ReplaceAll(str, ".0", "")

		writer.WriteString(str)
	}
	writer.Flush()
}
