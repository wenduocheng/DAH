package main

import (
	"bufio"
	"fmt"
	"os"
	"strconv"
)

func GetUniqueKmerCounts(kmerCounts map[string]int) map[int]int {
	kmerUniqCounts := make(map[int]int)

	for _, count := range kmerCounts {
		kmerUniqCounts[count]++
	}

	return kmerUniqCounts
}

//
func GenerateUniqueKmerCountTextFile(uniqueKmerCounts map[int]int, kmerLength int) {

	f, err := os.Create(strconv.Itoa(kmerLength) + ".txt")
	if err != nil {
		fmt.Println(err)
	}
	w := bufio.NewWriter(f)
	// var list []int
	for key, val := range uniqueKmerCounts {
		for i := 0; i < val; i++ {
			w.WriteString(strconv.Itoa(key) + ",")
		}
	}
	w.Flush()

}
func arraytextfile(GenerateReadsPlot []int, name string) {
	f, err := os.Create(name + ".txt")
	if err != nil {
		fmt.Println(err)
	}
	w := bufio.NewWriter(f)
	for idx := range GenerateReadsPlot {
		w.WriteString(strconv.Itoa(GenerateReadsPlot[idx]) + ",")
	}
}
