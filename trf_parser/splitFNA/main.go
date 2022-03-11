package main

import (
	"bufio"
	"fmt"
	"os"
	"strconv"
	"strings"
)

func main() {
	folder1 := "/home/ziheng/courses/CSE280/data/harish/ncbi-genomes-2022-03-04/compress_harish"
	folder2 := "/home/ziheng/courses/CSE280/data/candace/ncbi-genomes-2022-03-04/compress_candace"
	folder3 := "/home/ziheng/courses/CSE280/data/ncbi-genomes-2022-03-01/compress_ziheng"
	filename := "/compressedFNA.fna"

	vecFolder := []string{folder1, folder2, folder3}



	for _, folder := range vecFolder {
		// open FNA file
		fna, err := os.Open(folder + filename)
		if err != nil {
			fmt.Println(err)
			return
		}
		defer fna.Close()

		scannerFNA := bufio.NewScanner(fna)
		rowFNA := 0



		indexOutFile := 1
		// create output file
		outFNA, err := os.Create(folder + "/compressedFNA_Num_" + strconv.Itoa(indexOutFile) + "_part.fna")
		if err != nil {
			fmt.Println(err)
			panic(1)
		}
		writerOut := bufio.NewWriter(outFNA)

		rowOut := 0
		for scannerFNA.Scan() {
			rowFNA++
			str := scannerFNA.Text() + "\n"
			rowOut++

			if strings.Contains(str, ">") {
				if rowOut * 81 > 850000 {
					// time to have a new file
					rowOut = 0
					writerOut.Flush()

					outFNA.Close()
					indexOutFile++
					outFNA, err = os.Create(folder + "/compressedFNA_Num_" + strconv.Itoa(indexOutFile) + "_part.fna")
					if err != nil {
						fmt.Println(err)
						panic(1)
					}
					writerOut = bufio.NewWriter(outFNA)
				}
			}
			writerOut.WriteString(str)
		}

		writerOut.Flush()
		outFNA.Close()
	}
}
