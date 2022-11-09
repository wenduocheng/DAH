package main

import (
	"fmt"
	"math/rand"
	"testing"
)

// Wenduo
func TestGenerateSequence(t *testing.T) {
	rand.Seed(1.0)
	outcome := GenerateSequence(10)
	answer := "TGGGTCTAAA"
	if outcome != answer {
		t.Errorf("Error! The answer is %v, and your code gives %v", answer, outcome)
	} else {
		fmt.Println("Correct!")
	}

	outcome2 := GenerateSequence(0)
	answer2 := ""
	if outcome2 != answer2 {
		t.Errorf("Error! The answer is %v, and your code gives %v", answer, outcome)
	} else {
		fmt.Println("Correct!")
	}
}

// Wenduo
func TestGetKmers(t *testing.T) {
	sequence := "ATCGTTACCGT"
	kmerLength := 4
	kmers := GetKmers(sequence, kmerLength)
	answer := []string{"ATCG", "TCGT", "CGTT", "GTTA", "TTAC", "TACC", "ACCG", "CCGT"}

	if !SameStringSlices(kmers, answer) {
		t.Errorf("Error! The answer is %v, and your code gives %v", answer, kmers)
	} else {
		fmt.Println("Correct!")
	}

}

// Wenduo
func TestRandomMutate(t *testing.T) {
	rand.Seed(1.0)
	kmers := []string{"ATCG", "TCGT", "CGTT", "GTTA", "TTAC", "TACC", "ACCG", "CCGT"}
	outcome := RandomMutate(kmers, 3)
	answer := []string{"ATCG", "ACGT", "CGTT", "GCTA", "TTAC", "TACC", "ACCG", "CCGT"}

	if !SameStringSlices(outcome, answer) {
		t.Errorf("Error! The answer is %v, and your code gives %v", answer, outcome)
	} else {
		fmt.Println("Correct!")
	}
}

// Wenduo
func TestShuffle(t *testing.T) {
	originalKmers := []string{"ATCG", "TCGT", "CGTT", "GTTA", "TTAC", "TACC", "ACCG", "CCGT", "CCCC", "ATTC"}
	kmers := []string{"ATCG", "TCGT", "CGTT", "GTTA", "TTAC", "TACC", "ACCG", "CCGT", "CCCC", "ATTC"}
	outcome := Shuffle(kmers)
	if SameStringSlices(outcome, originalKmers) {
		t.Errorf("Error! The Shuffle function is not working.")
	} else {
		fmt.Println("Correct!")
	}
}

// Wenduo
func TestDelete(t *testing.T) {
	kmers := []string{"ATCG", "TCGT", "CGTT", "GTTA", "TTAC", "TACC", "ACCG", "CCGT", "CCCC", "ATTC"}
	outcome := Delete(kmers, 2, false)
	answer := []string{"ATCG", "TCGT", "CGTT", "GTTA", "TTAC", "TACC", "ACCG", "CCGT"}
	if !SameStringSlices(outcome, answer) {
		t.Errorf("Error! The answer is %v, and your code gives %v", answer, outcome)
	} else {
		fmt.Println("Correct!")
	}
}

// Wenduo
func TestGetUniqueKmerCounts(t *testing.T) {
	kmers := []string{"ATCG", "ATCG", "CGTT", "GTTA", "ATCG", "GTTA", "AAAA", "CCCC"}
	kmerCounts := KmerHash2(kmers)
	outcome := GetUniqueKmerCounts(kmerCounts)

	answer := make(map[int]int)
	answer[3] = 1
	answer[2] = 1
	answer[1] = 3

	if len(outcome) != 3 {
		t.Errorf("Error!")
	} else {
		fmt.Println("Correct!")
	}

	if outcome[3] != 1 || outcome[2] != 1 || outcome[1] != 3 {
		t.Errorf("Error! The answer is %v, and your code gives %v", answer, outcome)
	} else {
		fmt.Println("Correct!")
	}
}

// Wenduo
func TestKmerCountSort(t *testing.T) {
	kmers := []string{"ATCG", "ATCG", "CGTT", "GTTA", "ATCG", "GTTA", "AAAA", "CCCC"}
	kmerCounts := KmerHash2(kmers)
	kmerUniqCounts := GetUniqueKmerCounts(kmerCounts)
	outcome := KmerCountSort(kmerUniqCounts)
	answer := []int{1, 2, 3}

	if !SameIntegerSlices(answer, outcome) {
		t.Errorf("Error! The answer is %v, and your code gives %v", answer, outcome)
	} else {
		fmt.Println("Correct!")
	}
}
