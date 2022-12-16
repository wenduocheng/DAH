package main

// GenerateUniqueAndDistinctCount generate the unique counts and distinct counts map
// Input: The reads and number n represents the maximum kmer for the map
// Output: The unique and distinct kmers counts map, key: kmer length, value:kmer counts
// Lilin
func GenerateUniqueAndDistinctCount(reads []string, n int) (map[int]int, map[int]int) {
	var uniqueCount map[int]int
	var distinctCount map[int]int
	for k := 1; k <= n; k++ {
		kmerCounts := KmerHash(reads, k)
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
		if kmerCounts[k] == 1 {
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
func KmerHash(reads []string, k int) map[string]int {
	kmerCounts := make(map[string]int)
	for _, read := range reads {
		for i := 0; i <= len(read)-k; i++ {
			kmer := read[i : i+k]
			_, exists := kmerCounts[kmer]
			if exists {
				kmerCounts[kmer]++
			} else {
				kmerCounts[kmer] = 1
			}
		}
	}
	return kmerCounts
}
