package main

import (
	"bufio"
	"fmt"
	"math/rand"
	"os"
	"strconv"
	"time"
)

// GenerateReads randomly select positions of the genome and generate reads
// Wenduo
func GenerateReads(readLength int, readCounts int, sequence string) []string {
	n := len(sequence)
	reads := make([]string, 0)
	for i := 0; i < readCounts; i++ {
		start := rand.Intn(n)
		if start < n-readLength {
			read := sequence[start : start+readLength]
			reads = append(reads, read)
		}
	}
	reads = Noisify(reads, 300, 30)
	f, err := os.Create("GeneratedReads.fasta")
	if err != nil {
		fmt.Println(err)
	}
	w := bufio.NewWriter(f)
	// var list []int
	for i := range reads {

		w.WriteString(">reads" + strconv.Itoa(i) + "\n" + reads[i] + "\n")

	}
	w.Flush()

	return reads

}

// Generate Read Coverage Plot
// Wenduo
func GenerateReadsPlot(readLength int, readCounts int, sequence string) []int {
	n := len(sequence)
	reads := make([]int, n)
	for i := 0; i < readCounts; i++ {
		start := rand.Intn(n)
		if start < n-readLength {
			// read := sequence[start : start+readLength]
			// reads = append(reads, read)
			for j := start; j < start+readLength; j++ {
				reads[j]++
			}
		}
	}
	return reads
}

// GenerateSequence randomly genertes a sequence composed of A/T/C/G
// Wenduo
func GenerateSequence(length int) string {
	if length < 0 {
		panic("Error: Negative integer is given. Please give a nonnegative integer.")
	}
	nucleotides := []string{"A", "T", "C", "G"}
	sequence := ""
	for i := 0; i < length; i++ {
		index := rand.Intn(4)
		sequence += nucleotides[index]
	}
	return sequence
}

// GenerateReadsNaive will generate a list of reads.
// Input: an integer corresponding to read length, an integer corresponding to the number of copies of each read to be generated, sequence string
// Output: a list of strings corresponding to generated reads
// Wenduo
func GenerateReadsNaive(readLength int, numberOfCopies int, sequence string) []string {
	if readLength > len(sequence) {
		panic("Error: The read length has to be less than the sequence length.")
	}

	reads := make([]string, 0)
	for i := 0; i <= len(sequence)-readLength; i++ {
		reads = append(reads, sequence[i:i+readLength])
	}

	reads_copy := make([]string, len(reads))
	for j := 0; j < len(reads); j++ {
		reads_copy[j] = reads[j]
	}

	for k := 0; k < numberOfCopies-1; k++ {
		reads = append(reads, reads_copy...)
	}

	return reads
}

// GenerateReadsNaive will generate random reads by given read length. It will generte a list of strings with given readlength and different num of reads
// input: int,int
// tianyue
//
//output:reads []string
func GenerateReadsNaive0(readlength int, differentReads int) []string {
	rand.Seed(time.Now().UnixNano())
	var temp []string
	var LetterRunes = []rune("ATCG")
	Sequence := make([]rune, readlength)

	for j := 0; j < differentReads; j++ {
		for i := range Sequence {
			Sequence[i] = LetterRunes[rand.Intn(len(LetterRunes))]

		}
		//sequence will transfer rune type into string
		sequence := string(Sequence)
		temp = append(temp, sequence)
	}
	var ListSequence []string
	ListSequence = CopyReads(temp, 10)

	return ListSequence
}

// Copy each reads by given number of times
// input []string as sequence
// output []string reads with repeats
// tianyue
func CopyReads(ListSequence []string, times int) []string {
	var newReads []string
	for i := 0; i < times; i++ {
		newReads = append(newReads, ListSequence...)
	}
	return newReads
}

// Generate read seqs with the num pf each distinct read is normaly distrubuted
// input int int
// output []string as read seqs
// tianyue
func GenerateReadsNorm0(readlength, differentreads int) []string {
	var temp []string
	var LetterRunes = []rune("ATCG")
	// var NormalDist []int
	Sequence := make([]rune, readlength)

	for j := 0; j < differentreads; j++ {
		for i := range Sequence {
			Sequence[i] = LetterRunes[rand.Intn(len(LetterRunes))]

		}
		//sequence will transfer rune type into string
		sequence := string(Sequence)
		temp = append(temp, sequence)
	}

	//Generate a list number with Random Distribution
	var b []int
	for i := 0; i < len(temp); i++ {
		b = append(b, int(rand.NormFloat64()*5+10))
	}
	var SeqReads []string
	for j := 0; j < len(temp); j++ {
		for i := 0; i < len(b); i++ {

			SeqReads = append(SeqReads, temp[j])
		}
	}
	return SeqReads
}
func GenerateReadsNorm(readLength int, sequence string) []string {
	if readLength > len(sequence) {
		panic("Error: The read length has to be less than the sequence length.")
	}

	var temp_reads []string
	for i := 0; i <= len(sequence)-readLength; i++ {
		temp_reads = append(temp_reads, sequence[i:i+readLength])
	}

	//Generate a list number with Random Distribution
	var b []int
	for i := 0; i < len(temp_reads); i++ {
		b = append(b, int(rand.NormFloat64()*1+10))
	}

	var SeqReads []string
	for j := 0; j < len(temp_reads); j++ {
		for i := 0; i < b[j]; i++ {

			SeqReads = append(SeqReads, temp_reads[j])

		}
	}

	return SeqReads
}

// RandomMutate randomly mutates a slice of kmers
// It is possible that even after mutation, genome is still the same.
// Wenduo
func RandomMutate(kmers []string, numMutations int) []string {
	nucleotides := []string{"A", "T", "C", "G"}
	mutatedKmers := make([]string, len(kmers))
	// Copy kmers
	for index, val := range kmers {
		mutatedKmers[index] = val
	}
	for i := 0; i < numMutations; i++ {
		kmerIndex := rand.Intn(len(kmers))        // which kmer to be mutated
		positionIndex := rand.Intn(len(kmers[0])) // which position of a kmer to be mutated
		nucleotideIndex := rand.Intn(4)           // randomly choose A/T/C/G
		mutatedKmers[kmerIndex] = kmers[kmerIndex][0:positionIndex] + nucleotides[nucleotideIndex] + kmers[kmerIndex][positionIndex+1:]
	}
	return mutatedKmers
}

// Shuffle shuffles a slice of kmers
// Wenduo
func Shuffle(kmers []string) []string {
	rand.Seed(time.Now().UnixNano())
	rand.Shuffle(len(kmers), func(i, j int) { kmers[i], kmers[j] = kmers[j], kmers[i] })
	return kmers
}

// Delete deletes numDeletions kmers from a slice of kmers
// Wenduo
func Delete(kmers []string, numDeletions int, shuffle bool) []string {
	if shuffle {
		kmers = Shuffle(kmers)
	}
	for i := 0; i < numDeletions; i++ {
		kmers = kmers[:(len(kmers) - 1)]
	}
	return kmers
}

// Noisify introduces noise to a set of kmers by shuffling, mutating and deleting them
// Wenduo
func Noisify(kmers []string, numMutations, numDeletions int) []string {
	deletedkmers := Delete(kmers, numDeletions, true)
	mutatedKmers := RandomMutate(deletedkmers, numMutations)
	return mutatedKmers
}
