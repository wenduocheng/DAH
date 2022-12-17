package main

import (
	"math"
)

// N50 takes the list of contigs and return the shortest contig length needed to cover 50% genome,
// describing the completeness of genome assembly
// Input: the list of contigs
// Output: The length of the shortest contig cover 50% genome
// Lilin
func N50(contigs []string) int { //contigs
	var total_length int
	contigs = Sort(contigs)
	for _, contig := range contigs {
		total_length += len(contig)
	}

	N50_length := (float64(total_length) / 2.0)
	var currLength int
	var N50_string string
	for _, contig := range contigs {
		currLength += len(contig)
		if float64(currLength) >= N50_length {
			N50_string = contig
			break
		}
	}

	return len(N50_string)
}

// Sort sorts the contigs according to their length from longest to shortest
// Input: the list of contigues
// Output: the list of the sorted contigs from longest to shortest
// Lilin
func Sort(contigs []string) []string {
	lenHash := make(map[int][]string)
	for _, contig := range contigs {
		length := len(contig)
		_, Exists := lenHash[length]
		if Exists {
			lenHash[length] = append(lenHash[length], contig)
		} else {
			lenHash[length] = []string{contig}
		}
	}
	var lengthList []int
	for key, _ := range lenHash {
		lengthList = append(lengthList, key)
	}
	lengthList = DivideConquer(lengthList)
	var sortList []string
	for _, num := range lengthList {
		stringList := lenHash[num]
		sortList = append(sortList, stringList...)
	}
	return sortList
}

// DivideConquer sorts the list of int from largest to smallest using divide&conquer algorithm
// Intput: a list of integer
// Output: a sortest list of integer from largest to smallest
// Lilin
func DivideConquer(numlist []int) []int {
	if len(numlist) == 1 || len(numlist) == 0 {
		return numlist
	}
	i := int(len(numlist) / 2)
	right := DivideConquer(numlist[:i])
	left := DivideConquer(numlist[i:])
	var combineList []int
	for len(right) != 0 && len(left) != 0 {
		if right[0] >= left[0] {
			combineList = append(combineList, right[0])
			right = right[1:]
		} else if right[0] < left[0] {
			combineList = append(combineList, left[0])
			left = left[1:]
		}
	}
	if len(right) > 0 {
		combineList = append(combineList, right...)
	}
	if len(left) > 0 {
		combineList = append(combineList, left...)
	}
	return combineList
}

// DistinctKmerCount counts the number of distinct kmer in the string
// Input: a map of kmer and frequency
// Output: the number of distinct kmer
// Lilin
func DistinctKmerCount(kmerCounts map[string]int) int {

	count := len(kmerCounts)

	return count
}

// KmerSizeAdjustment counts the kmer size according to the readlength and kmer and genome coverage
// Input: the readlength, kmer coverage and genome coverage
// Output: the estimate optimal kmer size
// Lilin
func KmerSizeAdjustment(readlength int, kmer_coverage, genome_coverage float64) int {

	var kmer_size int
	kmer_size = 1 + readlength - int((kmer_coverage*float64(readlength))/genome_coverage)

	return kmer_size
}

// GenomeCoverage counts the genome coverage for the current kmer set
// Input: the estimate genome length, the readlength and the reads number
// Output: the genome coverage
// Lilin
func GenomeCoverage(genome_length, readlength, reads_number int) float64 {
	var genome_coverage float64

	//(A total number of reads * read length)/ (Estimated genome size)
	genome_coverage = float64(reads_number) * float64(readlength) / float64(genome_length)

	return genome_coverage
}

// KmerCoverage counts the kmer coverage for the current kmer set
// Input: the read length, kmer size and genome coverage
// Output: the kmer coverage
// Lilin
func KmerCoverage(readlength, kmersize int, genome_coverage float64) float64 {

	var kmer_coverage float64
	base_coverage := genome_coverage
	//𝐶𝑘=𝐶⋅(𝑅−𝐾+1)/𝑅
	kmer_coverage = base_coverage * float64(readlength-kmersize+1) / float64(readlength)

	return kmer_coverage
}

// KmerSizeSet takes the readlength and generate the set of interger for range over to pick the optimal kmer size form this set
// Input: the readlength
// Output: the generated kmer size set
// Lilin
func KmerSizeSet(readlength int) []int {
	var size_set []int
	//find the prome number less or equal to the readlength
	size_set = FindPrime(readlength)
	return size_set
}

// FindPrime takes in the value and returns the list of prime number smaller than that value
// in this case, although 1 is not prime number, we also take it into consideration
// Input: the value
// Output: the int list contains the prime number smaller than the value
// Lilin
func FindPrime(n int) []int {
	var prime_list []int
	//special case even 1 is not prime
	if n == 1 {
		prime_list = append(prime_list, 1)
		return prime_list
	}
	if n == 2 {
		prime_list = append(prime_list, 1)
		prime_list = append(prime_list, 2)
		return prime_list
	}

	var number_list []int
	for i := 0; i <= n; i++ {
		number_list = append(number_list, i)
	}

	var if_prime_list []bool
	if_prime_list = make([]bool, n+1)
	//1 and 2 are prime
	if_prime_list[0] = false
	if_prime_list[1] = true

	for i := 2; i <= n; i++ {
		if_prime_list[i] = true
	}
	sqrtn := int(math.Sqrt(float64(n))) + 1
	//cross out multiple of the prime numebr
	for i := 2; i < sqrtn; i++ {
		if if_prime_list[i] {
			for j := i * i; j <= n; j += i {
				if_prime_list[j] = false
			}
		}
	}

	// if the nth position of the if_prime_list is true,
	// the number at nth position of number_list is prime
	for n := range if_prime_list {
		if if_prime_list[n] {
			prime_list = append(prime_list, number_list[n])
		}
	}

	return prime_list
}

// OptimalKmerSize gives the optiaml kmersize for the later analysis
// Input: genome, list of reads
// Output: optimal kmer size
// Lilin
func OptimalKmerSize(reads []string, coverage int) int {

	var optimalk int

	readlength := len(reads[0])
	//method1: find the k gives most number of distinct kmers
	k_size_set := KmerSizeSet(readlength)
	max_distinct_kmer_count := 0
	var optimalk1 int
	for i := range k_size_set {
		k := k_size_set[i]

		var kmerCounts map[string]int

		kmerCounts = KmerHashFromReads(k, reads)

		distinct_kmer_count := DistinctKmerCount(kmerCounts)
		if distinct_kmer_count > max_distinct_kmer_count {
			max_distinct_kmer_count = distinct_kmer_count
			optimalk1 = k
		}
	}

	//methods improve result of first method from genome and kmer coverage
	var kmer_coverage float64
	var genome_coverage float64
	var readsnumber int
	var genome_length int

	genome_length = GenomeSizeEstimate(reads, coverage)
	readsnumber = len(reads)
	genome_coverage = GenomeCoverage(genome_length, readlength, readsnumber)
	kmer_coverage = KmerCoverage(readlength, optimalk1, genome_coverage)

	optimalk = KmerSizeAdjustment(readlength, kmer_coverage, genome_coverage)

	//if the kmer is too small, the memory would be large, and mismatch could also be large
	if readlength >= 10 && optimalk < 4 {
		optimalk = 4
	}
	//if too large accuracy problem
	if optimalk > 25 {
		optimalk = 25
	}
	return optimalk

}

// GenomeSizeEstimate estimates the size of the genome
// Input: the reads set, the coverage rate
// Output: the estimate size of genome
//Lilin
func GenomeSizeEstimate(reads []string, numberOfCopies int) int {
	var estimate_size int
	//our reads is 10X
	estimate_size = len(reads) / numberOfCopies
	return estimate_size
}

//OptimalKmerSizeWithRange finds the best kmer length within a range
//Input: the range of the kmerlength the min value and the max value, the slice of reads
//Output: the kmer length within this range of best fit
//Lilin
func OptimalKmerSizeWithRange(kmermin, kmermax int, reads []string) int {
	var k_size_set []int

	for i := kmermin; i <= kmermax; i++ {
		k_size_set = append(k_size_set, i)
	}

	max_distinct_kmer_count := 0
	var optimalk1 int
	for i := range k_size_set {
		k := k_size_set[i]

		var kmerCounts map[string]int

		kmerCounts = KmerHashFromReads(k, reads)

		distinct_kmer_count := DistinctKmerCount(kmerCounts)
		if distinct_kmer_count > max_distinct_kmer_count {
			max_distinct_kmer_count = distinct_kmer_count
			optimalk1 = k
		}
	}
	return optimalk1
}

// GenerateUniqueAndDistinctCount generate the unique counts and distinct counts map
// Input: The reads and number n represents the maximum kmer for the map
// Output: The unique and distinct kmers counts map, key: kmer length, value:kmer counts
// Lilin
func GenerateUniqueAndDistinctCount(genome string, n int) (map[int]int, map[int]int) {
	uniqueCount := make(map[int]int)
	distinctCount := make(map[int]int)
	for k := 1; k <= n; k++ {
		kmerCounts := KmerHash(genome, k)
		countu := UniqueCount(kmerCounts)
		countd := DistinctCount(kmerCounts)

		uniqueCount[k] = countu
		distinctCount[k] = countd
	}
	return uniqueCount, distinctCount
}

// UniqueCount finds the unique kmer counts
// Intput: kmerCounts map, key: kmers, value: counts of this kmer
// Output: the number of unique kmers
// Lilin
func UniqueCount(kmerCounts map[string]int) int {
	var count int

	for k, _ := range kmerCounts {
		if kmerCounts[k] <= 10 {
			count++
		}
	}
	return count
}

// DistinctCount finds the distinct kmer counts
// Input: kmerCounts map, key: kmers, value: counts of this kmer
// Output: the number of distinct kmers
// Lilin
func DistinctCount(kmerCounts map[string]int) int {
	var count int
	count = len(kmerCounts)
	return count
}

// KmerHash hash the kmers to the map and record the number of time that kmer appears in the genome
// Input: the reads and the kmer length
// Output: the Kmer hash map, key: kmers value: kmer counts
// Lilin
func KmerHash(genome string, k int) map[string]int {
	kmerCounts := make(map[string]int)

	for i := 0; i <= len(genome)-k; i++ {
		kmer := genome[i : i+k]
		_, exists := kmerCounts[kmer]
		if exists {
			kmerCounts[kmer]++
		} else {
			kmerCounts[kmer] = 1
		}
	}

	return kmerCounts
}
