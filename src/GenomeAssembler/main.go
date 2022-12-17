package main

import (
	"bufio"
	"fmt"
	"log"
	"os"
	"strconv"
)
//Sample input in command line 
//Input is Genome:
//./run genome default
//./run genome range 21 30

// }

func main() {
	//os.Args[1] is going to be if it is genome or reads file
	filetype := os.Args[1]
	if filetype != "genome" && filetype != "reads" {
		panic("Wrong filetype, try again.")
	}
	//os.Args[2] is going to be the file name
	filename := os.Args[2]

	//os.Args[3] is going to be the kmer choice
	kmerchoice := os.Args[3]
	if kmerchoice != "default" && kmerchoice != "range" {
		panic("Wrong kmerchoice, try again.")
	}

	//if kmerchoice is range: we have two more Args to represent min and max of the range
	var min, max int
	if kmerchoice == "range" {
		//min
		min1, err1 := strconv.Atoi(os.Args[4])
		if err1 != nil {
			panic(err1)
		}
		min = min1
		fmt.Println(min1)
		if min1 <= 0 {
			panic("kmer min less or equal to 0 wrong, try again.")
		}

		//max
		max1, err1 := strconv.Atoi(os.Args[5])
		if err1 != nil {
			panic(err1)
		}
		max = max1
		if max1 < min {
			panic("max is less than min, wrong sequence, try again")
		}
	}

	var reads []string
	if filetype == "genome" {
		genome := ReadSequence(filename)
		readLength := 100
		readCounts := 300000
		reads = GenerateReads(readLength, readCounts, genome)
		readsForPlot := GenerateReadsPlot(readLength, readCounts, viralSequence)
		DrawBarPlot(readsForPlot)
	} else {
		reads = ReadReads(filename)
	}

	//generate the distinct and unique kmer counts plot:
	// uniquecounts, distinctcounts := GenerateUniqueAndDistinctCount(reads, 15)

	//output the kmer counts to draw in python
	// fmt.Println("Now output the kmer counts.")

	// //unique counts file
	// file1, err := os.OpenFile("UniqueKmerCounts.txt", os.O_APPEND|os.O_CREATE|os.O_WRONLY, 0644)
	// if err != nil {
	// 	log.Fatalf("failed creating file: %s", err)
	// }

	// datawriter1 := bufio.NewWriter(file1)

	// for k, _ := range uniquecounts {
	// 	_, _ = datawriter1.WriteString(strconv.Itoa(k) + ":" + strconv.Itoa(uniquecounts[k]) + "\n")
	// }

	// datawriter1.Flush()
	// file1.Close()

	// file2, err := os.OpenFile("DistinctKmerCounts.txt", os.O_APPEND|os.O_CREATE|os.O_WRONLY, 0644)
	// if err != nil {
	// 	log.Fatalf("failed creating file: %s", err)
	// }

	// datawriter2 := bufio.NewWriter(file2)

	// for k, _ := range distinctcounts {
	// 	_, _ = datawriter2.WriteString(strconv.Itoa(k) + ":" + strconv.Itoa(distinctcounts[k]) + "\n")
	// }

	// datawriter2.Flush()
	// file2.Close()

	//kmer selection
	var kmer_length int
	var coverage int
	coverage = 10
	if kmerchoice == "default" {
		fmt.Println("The kmer choice is default, now select the adjusted optimal kmer length.")
		kmer_length = OptimalKmerSize(reads, coverage)
	} else {
		fmt.Println("The kmer choice is range, the range is from", min, "to", max)
		kmer_length = OptimalKmerSizeWithRange(min, max, reads)
	}

	//print out the optiaml kmerlength
	fmt.Println("The optimal kmerlegnth is:")
	fmt.Println(kmer_length)

	//contigs assembly
	fmt.Println("Now perform the de novo assembly.")
	contigs := DenovoAssembler(reads, kmer_length)

	fmt.Println("De novo assembly was finished!")

	//output the generated contigs
	fmt.Println("Now output the contigs.")
	file3, err := os.OpenFile("contigs.fasta", os.O_APPEND|os.O_CREATE|os.O_WRONLY, 0644)

	if err != nil {
		log.Fatalf("failed creating file: %s", err)
	}

	datawriter3 := bufio.NewWriter(file3)

	for i, contig := range contigs {
		_, _ = datawriter3.WriteString(">contigs" + strconv.Itoa(i) + " length: " + strconv.Itoa(len(contig)) + "\n")
		_, _ = datawriter3.WriteString(contig + "\n")
	}

	datawriter3.Flush()
	file3.Close()
	fmt.Println("Contigs are written to the file")

	//get the N50 length
	n50 := N50(contigs)
	fmt.Print("The N50 is: ")
	fmt.Println(n50)

	fmt.Println("Done")

}

