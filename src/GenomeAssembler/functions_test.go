package main

import (
	"fmt"
	"math/rand"
	"testing"
)

// Wenduo
func TestKmerHashFromReads(t *testing.T) {
	kmers := []string{"ATCG", "TCGT", "CGTT", "GTTA", "TTAC", "TACC", "ACCG", "CCGT", "CCCC", "ATTC"}
	kmerCounts := KmerHashFromReads(3, kmers)

	if len(kmerCounts) != 11 {
		t.Errorf("Error! KmerHashFromReads miss out some kmers!")
	} else if kmerCounts["ACC"] != 2 {
		t.Errorf("Error! KmerHashFromReads gives the wrong counts!")
	} else if kmerCounts["ATC"] != 1 {
		t.Errorf("Error! KmerHashFromReads gives the wrong counts!")
	} else if kmerCounts["TCG"] != 2 {
		t.Errorf("Error! KmerHashFromReads gives the wrong counts!")
	} else if kmerCounts["TTC"] != 1 {
		t.Errorf("Error! KmerHashFromReads gives the wrong counts!")
	} else {
		fmt.Println("KmerHashFromReads passes all tests!!")
	}
}

// Wenduo
func TestGetKmerComposition(t *testing.T) {
	kmers := []string{"ATCG", "TCGT", "CGTT", "GTTA", "TTAC", "TACC", "ACCG", "CCGT", "CCCC", "ATTC"}
	kmerCounts := KmerHashFromReads(3, kmers)
	kmerComposition := GetKmerComposition(kmerCounts)
	fmt.Println(kmerComposition)

	if len(kmerCounts) != 11 {
		t.Errorf("Error! GetKmerComposition miss out some kmers!")
	} else if !Contain(kmerComposition, "TCG") {
		t.Errorf("Error! KmerHashFromReads should have contained TCG!")
	} else if !Contain(kmerComposition, "TAC") {
		t.Errorf("Error! KmerHashFromReads should have contained TAC!")
	} else if !Contain(kmerComposition, "ATT") {
		t.Errorf("Error! KmerHashFromReads should have contained ATT!")
	} else if !Contain(kmerComposition, "TTC") {
		t.Errorf("Error! KmerHashFromReads should have contained TTC!")
	} else if !Contain(kmerComposition, "CGT") {
		t.Errorf("Error! KmerHashFromReads should have contained CGT!")
	} else if !Contain(kmerComposition, "TTA") {
		t.Errorf("Error! KmerHashFromReads should have contained TTA!")
	} else {
		fmt.Println("KmerHashFromReads passes all tests!!")
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
		fmt.Println("RandomMutate passes all tests!!")
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
		fmt.Println("Shuffle passes all tests!!")
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
		fmt.Println("Delete passes all tests!")
	}
}

// Wenduo
func TestGetUniqueKmerCounts(t *testing.T) {
	kmers := []string{"ATCG", "ATCG", "CGTT", "GTTA", "ATCG", "GTTA", "AAAA", "CCCC"}
	kmerCounts := KmerHashFromReads(4, kmers)
	outcome := GetUniqueKmerCounts(kmerCounts)

	answer := make(map[int]int)
	answer[3] = 1
	answer[2] = 1
	answer[1] = 3

	if len(outcome) != 3 {
		t.Errorf("Error!")
	}
	// } else {
	// 	fmt.Println("Correct!")
	// }

	if outcome[3] != 1 || outcome[2] != 1 || outcome[1] != 3 {
		t.Errorf("Error! The answer is %v, and your code gives %v", answer, outcome)
	} else {
		fmt.Println("GetUniqueKmerCounts passes all tests!")
	}
}

// Wenduo
func TestKmerCountSort(t *testing.T) {
	kmers := []string{"ATCG", "ATCG", "CGTT", "GTTA", "ATCG", "GTTA", "AAAA", "CCCC"}
	kmerCounts := KmerHashFromReads(4, kmers)
	kmerUniqCounts := GetUniqueKmerCounts(kmerCounts)
	outcome := KmerCountSort(kmerUniqCounts)
	answer := []int{1, 2, 3}

	if !SameIntegerSlices(answer, outcome) {
		t.Errorf("Error! The answer is %v, and your code gives %v", answer, outcome)
	} else {
		fmt.Println("KmerCountSort passes all tests!")
	}
}

// Wenduo
func TestGenerateSequence(t *testing.T) {
	rand.Seed(1.0)
	outcome := GenerateSequence(10)
	answer := "TGGGTCTAAA"
	if outcome != answer {
		t.Errorf("Error! The answer is %v, and your code gives %v", answer, outcome)
	} else {
		fmt.Println("GenerateSequence passes first test!")
	}

	outcome2 := GenerateSequence(0)
	answer2 := ""
	if outcome2 != answer2 {
		t.Errorf("Error! The answer is %v, and your code gives %v", answer, outcome)
	} else {
		fmt.Println("GenerateSequence passes second test!")
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
		fmt.Println("GetKmers pass all tests!")
	}

}

// Lilin
func TestDivideConquer(t *testing.T) {

	//base case
	test1 := []int{1, 3, 2}
	result1 := DivideConquer(test1)
	expect1 := []int{3, 2, 1}
	for i, _ := range result1 {
		if result1[i] != expect1[i] {
			t.Errorf("DivideConquer base case1 fails")
		}
	}
	test2 := []int{6}
	result2 := DivideConquer(test2)
	expect2 := []int{6}
	for i, _ := range result2 {
		if result2[i] != expect2[i] {
			t.Errorf("DivideConquer base case2 fails")
		}
	}
	test3 := []int{2, 3, 4, 5, 6}
	result3 := DivideConquer(test3)
	expect3 := []int{6, 5, 4, 3, 2}
	for i, _ := range result3 {
		if result3[i] != expect3[i] {
			t.Errorf("DivideConquer base case3 fails")
		}
	}

	//duplicate case
	test4 := []int{1, 5, 7, 7, 2, 2, 3, 3}
	result4 := DivideConquer(test4)
	expect4 := []int{7, 7, 5, 3, 3, 2, 2, 1}
	for i, _ := range result4 {
		if result4[i] != expect4[i] {
			t.Errorf("DivideConquer duplicate case fails")
		}
	}

	fmt.Println("DivideConquer passes all tests!")
}

// Lilin
func TestSort(t *testing.T) {
	//base case
	test1 := []string{"A", "AC", "ACT"}
	result1 := Sort(test1)
	expect1 := []string{"ACT", "AC", "A"}
	for i, _ := range result1 {
		if len(result1[i]) != len(expect1[i]) {
			t.Errorf("Sort base case fails")
		}
	}
	//duplicate case
	test2 := []string{"AGGT", "ACCT", "AGGG", "CTGCGC", "AA", "GGG"}
	result2 := Sort(test2)
	expect2 := []string{"CTGCGC", "AGGT", "ACCT", "AGGG", "GGG", "AA"}
	for i, _ := range result2 {
		if len(result2[i]) != len(expect2[i]) {
			t.Errorf("Sort duplicate case fails")
		}
	}
	fmt.Println("Sort passes all tests!")
}

// Lilin
func TestN50(t *testing.T) {
	test1 := []string{"AGG"}
	result1 := N50(test1)
	expect1 := 3
	if result1 != expect1 {
		t.Errorf("N50 fails at list of length 1")
	}

	test2 := []string{"ATT", "AGCT", "AGGGG", "AGCCG", "ATTT"}
	result2 := N50(test2)
	expect2 := 4
	if result2 != expect2 {
		t.Errorf("N50 fails at duplicate length")
	}

	fmt.Println("N50 passes all tests!")
}

// Lilin
func TestEulerianPath(t *testing.T) {
	//base case
	var text1 []string
	text1 = append(text1, "AAGATT")
	graph1 := DeBruijnGraph(4, text1)
	path1 := EulerianPath(graph1)
	genome1 := ReconstructStringFromGenomePath(path1[0], 4)
	if genome1 != text1[0] {
		t.Errorf("EulerianPath fails at base case")
	}
	//genome length = kmer length
	var text2 []string
	text2 = append(text2, "AGCT")
	graph2 := DeBruijnGraph(4, text2)
	path2 := EulerianPath(graph2)
	genome2 := ReconstructStringFromGenomePath(path2[0], 4)
	if genome2 != text2[0] {
		t.Errorf("EulerianPath fails at kmer length = genome length")
	}

	fmt.Println("EulerianPath passes all tests!")

}

// Lilin
func TestFindPrime(t *testing.T) {
	a := 20
	a_prime_list := []int{1, 2, 3, 5, 7, 11, 13, 17, 19}
	a_findprime := FindPrime(a)
	for i, _ := range a_prime_list {
		if a_prime_list[i] != a_findprime[i] {
			t.Errorf("FindPrime function is wrong!")
		}
	}

	b := 1
	b_prime_list := []int{1}
	b_findprime := FindPrime(b)
	for i, _ := range b_prime_list {
		if b_prime_list[i] != b_findprime[i] {
			t.Errorf("FindPrime function is wrong!")
		}
	}
	fmt.Println("FindPrime passes all tests!")
}

// Lilin
func TestDistinctKmerCount(t *testing.T) {

	var text1 []string
	text1 = append(text1, "AAGAAGAAA")
	expectedDisCount1 := 4
	kmerCounts1 := KmerHashFromReads(3, text1)
	realDisCount1 := DistinctKmerCount(kmerCounts1)
	if expectedDisCount1 != realDisCount1 {
		t.Errorf("DistinctKmerCount is wrong!")
	}

	var text2 []string
	text2 = append(text2, "AGTCGCTA")
	expectedDisCount2 := 5
	kmerCounts2 := KmerHashFromReads(4, text2)
	realDisCount2 := DistinctKmerCount(kmerCounts2)
	if expectedDisCount2 != realDisCount2 {
		t.Errorf("DistinctKmerCount is wrong!")
	}
	fmt.Println("DistinctKmerCount passes all tests!")
}

// Lilin
func TestGenomeCoverage(t *testing.T) {
	genome := "AGTCAGTCAGTCAGTCAGTC" //20
	readlength := 5
	reads_number := 10
	genome_length := len(genome)
	genome_coverage := GenomeCoverage(genome_length, readlength, reads_number)

	if genome_coverage != 2.5 {
		t.Errorf("GenomeCoverage is wrong!")
	}

	fmt.Println("GenomeCoverage passes all tests!")
}

//Lilin
func TestKmerHash(t *testing.T) {
	var reads1 []string
	reads1 = append(reads1, "AGTC")
	reads1 = append(reads1, "GGCT")

	result1 := KmerHash(reads1, 3)
	expect1 := map[string]int{"AGT": 1, "GCT": 1, "GGC": 1, "GTC": 1}

	if result1["AGT"] != expect1["AGT"] || result1["GCT"] != expect1["GCT"] {
		t.Errorf("KmerHash fails at base case!")
	}

	var reads2 []string
	reads2 = append(reads2, "GGGG")
	reads2 = append(reads2, "GGCT")

	result2 := KmerHash(reads2, 3)
	expect2 := map[string]int{"GGG": 2, "GGC": 1, "GTC": 1}

	if result2["GGG"] != expect2["GGG"] || result2["GGC"] != expect2["GGC"] {
		t.Errorf("KmerHash fails at repeat case!")
	}

	fmt.Println("KmerHash passes all tests")
}
//tianyue
func TestGetInputForHistogram(t *testing.T) {
	input := map[int]int{
		1: 3,
		2: 2,
		3: 1,
	}
	output := GetInputForHistogram(input)

	answer := []int{1, 1, 1, 2, 2, 3}
	if len(output) != len(answer) {
		t.Errorf("Error! Length of output doesn't match with answer!")
	} else {
		for i := 0; i < len(answer); i++ {
			if output[i] != answer[i] {
				t.Errorf("Error! Output doesn't match with answer!")
			}
		}

	}
	fmt.Println("GetInputForHistogram passes all tests!")
}

//tianyue
func TestCopyReads(t *testing.T) {
	reads := []string{
		"ATCGC", "TTCGA", "ACTTA",
	}
	times := 3

	output := CopyReads(reads, times)
	answer := []string{
		"ATCGC", "TTCGA", "ACTTA", "ATCGC", "TTCGA", "ACTTA", "ATCGC", "TTCGA", "ACTTA",
	}
	if len(output) != len(answer) {
		t.Errorf("Error! Length of output is wrong!")
	} else {
		for i := 0; i < len(output); i++ {
			if output[i] != answer[i] {
				t.Errorf("Error! Output is wrong!")
			}
		}
	}
	fmt.Println("CopyReads passes all tests!")
}
//tianyue
func TestGetRidOfN(t *testing.T) {
	Input := "ACGGANA"
	output := GetRidOfN(Input)
	answer := "ACGGAA"
	if output != answer {
		t.Errorf("Error! Output1 is wrong!")
	}
	Input2 := "AGCNNNNN"
	output2 := GetRidOfN(Input2)
	answer2 := "AGC"
	if output2 != answer2 {
		t.Errorf("Error! Output2 is wrong!")
	} else {
		fmt.Println("GetRidOfN passes all tests!")
	}

}
