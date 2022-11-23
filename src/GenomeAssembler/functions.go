package main

import (
	"math/rand"
	"sort"
	"time"
)

// String Composition Problem: Generate the k-mer composition of a string.
// Input: An integer k and a string Text.
// Output: Compositionk(Text), where the k-mers are arranged in lexicographic order.
// Wenduo
func StringComposition(k int, text string) []string {
	//count how many times each kmer occurred in this sequence
	kmerCounts := KmerHash(k, text)

	kmerComposition := make([]string, 0, len(kmerCounts))
	for kmer := range kmerCounts {
		kmerComposition = append(kmerComposition, kmer)
	}

	return kmerComposition
}

// KmerHash
// Input: An integer k and a string Text.
// Output: a map of kmer and frequency
// Wenduo
func KmerHash(k int, text string) map[string]int {
	//count how many times each kmer occurred in this sequence
	kmerCounts := make(map[string]int)
	for j := 0; j <= len(text)-k; j++ {
		kmer := text[j : j+k]
		_, exists := kmerCounts[kmer]
		if exists {
			kmerCounts[kmer]++
		} else {
			kmerCounts[kmer] = 1
		}
	}

	return kmerCounts
}

func KmerHash2(kmers []string) map[string]int {
	//count how many times each kmer occurred in this sequence
	kmerCounts := make(map[string]int)
	for _, kmer := range kmers {
		_, exists := kmerCounts[kmer]
		if exists {
			kmerCounts[kmer]++
		} else {
			kmerCounts[kmer] = 1
		}
	}

	return kmerCounts
}

// String Spelled by a Genome Path Problem. Reconstruct a string from its genome path.
// Input: A sequence path of k-mers Pattern1, … ,Patternn such that the last k - 1 symbols of Patterni are equal to the first k-1 symbols of Patterni+1 for 1 ≤ i ≤ n-1.
// Output: A string Text of length k+n-1 such that the i-th k-mer in Text is equal to Patterni (for 1 ≤ i ≤ n).
// Wenduo
func ReconstructStringFromGenomePath(sequencePath []string) string {
	genome := sequencePath[0]
	k := len(sequencePath[0])
	for i := 1; i < len(sequencePath); i++ {
		genome = genome + sequencePath[i][k-1:k]
	}
	return genome
}

// Overlap Graph Problem
// Construct the overlap graph of a collection of k-mers.
// Input: A collection Patterns of k-mers.
// Output: The overlap graph Overlap(Patterns), in the form of an adjacency matix/list.
// Wenduo
func OverlapGraph(kmersCollection []string) [][]int {
	if len(kmersCollection) <= 1 {
		panic("Erros: too few kmers")
	}
	adjacencyMatrix := make([][]int, len(kmersCollection))
	for i := range kmersCollection {
		adjacencyMatrix[i] = make([]int, len(kmersCollection))
		for j := range kmersCollection {
			adjacencyMatrix[i][j] = 0
			if i != j && Suffix(kmersCollection[i]) == Prefix(kmersCollection[j]) {
				adjacencyMatrix[i][j] = 1
			}
		}
	}
	return adjacencyMatrix
}

// Suffix
// Given: A string respresenting a kmer
// Return: A string corresponding to the suffix of the kmer
// Wenduo
func Suffix(kmer string) string {
	if len(kmer) <= 1 {
		panic("Error: kmer size too small")
	}
	return kmer[1:len(kmer)]
}

// Prefix
// Given: A string respresenting a kmer
// Return: A string corresponding to the prefix of the kmer
// Wenduo
func Prefix(kmer string) string {
	if len(kmer) <= 1 {
		panic("Error: kmer size too small")
	}
	return kmer[0:(len(kmer) - 1)]
}

// De Bruijn Graph from a String Problem
// Construct the de Bruijn graph of a string.
// Input: An integer k and a string Text.
// Output: DeBruijn_k(Tejxt), in the form of an adjacency matrix/list.
// Wenduo; Lilin
func DeBruijnGraph(k int, sequence string) map[string][]string {
	kmerComposition := StringComposition(k, sequence)

	adjacencyList := make(map[string][]string)
	for _, kmer := range kmerComposition {
		node := Prefix(kmer)

		adjacencyList[node] = append(adjacencyList[node], Suffix(kmer))
		// _, exists := adjacencyList[node]
		// if exists {
		// 	adjacencyList[node] = append(adjacencyList[node], Suffix(kmer))
		// } else {
		// 	adjacencyList[node] = append(adjacencyList[node], Suffix(kmer))
		// }
	}
	return adjacencyList
}

// De Bruijn Graph from a String Problem
// Construct the de Bruijn graph of a string.
// Input: An integer k and a string Text.
// Output: a graph
// Wenduo; Lilin
func DeBruijnGraph2(k int, text string) Graph {
	var dbgraph Graph
	dbnodes := make(map[string]*Node)
	dbedges := make(map[string]*Edge)

	kmerComposition := StringComposition(k, text) // list of kmers
	kmers := KmerHash(k, text)                    // map: kmer and frequency

	for _, kmer := range kmerComposition {
		// add prefix of kmer to nodes
		_, prefixExists := dbnodes[Prefix(kmer)]

		if prefixExists {
			dbnodes[Prefix(kmer)].outDegree++
		} else {
			var newFromNode Node
			newFromNode.label = Prefix(kmer)
			newFromNode.outDegree++
			dbnodes[Prefix(kmer)] = &newFromNode

		}
		// add suffix of kmer to nodes
		_, suffixExists := dbnodes[Suffix(kmer)]
		if suffixExists {
			dbnodes[Suffix(kmer)].inDegree++
		} else {
			var newToNode Node
			newToNode.label = Suffix(kmer)
			newToNode.inDegree++
			dbnodes[Suffix(kmer)] = &newToNode
		}
		// add the edge (prefix->suffix) to edges
		var newEdge Edge
		newEdge.from = dbnodes[Prefix(kmer)]
		newEdge.to = dbnodes[Suffix(kmer)]
		newEdge.weight = kmers[kmer]
		newEdge.label = kmer

		dbedges[kmer] = &newEdge

		// add sufix node as a child to prefix node
		dbnodes[Prefix(kmer)].children = append(dbnodes[Prefix(kmer)].children, dbnodes[Suffix(kmer)])

	}
	dbgraph.nodes = dbnodes
	dbgraph.edges = dbedges
	dbgraph.root = dbgraph.nodes[Prefix(kmerComposition[0])]
	return dbgraph
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

// GenerateSequence returns all the kmers of a string given a kmerLength
// Wenduo
func GetKmers(sequence string, kmerLength int) []string {
	kmers := make([]string, 0)
	for j := 0; j <= len(sequence)-kmerLength; j++ {
		kmer := sequence[j:(j + kmerLength)]
		kmers = append(kmers, kmer)
	}
	return kmers
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

// GetUniqueKmerCounts
// Input: map of kmer to frequency
// Output: a list of unique kmer counts (and the number of kmers with this count) present in the kmerCounts map.
// Wenduo
func GetUniqueKmerCounts(kmerCounts map[string]int) map[int]int {
	kmerUniqCounts := make(map[int]int)

	for _, count := range kmerCounts {
		kmerUniqCounts[count]++
	}

	return kmerUniqCounts
}

// KmerCountSort
// Input: kmerUniqCounts
// Output: a sorted list of kmer count values
// to generate distribution plot
// Wenduo
func KmerCountSort(kmerUniqCounts map[int]int) []int {
	sortedCounts := make([]int, 0)
	for count := range kmerUniqCounts {
		sortedCounts = append(sortedCounts, count)
	}
	sort.Ints(sortedCounts)
	return sortedCounts
}

// SameStringSlices returns true if two slices are the same
// Wendu
func SameStringSlices(s1, s2 []string) bool {
	if len(s1) != len(s2) {
		return false
	}
	for i := range s1 {
		if s1[i] != s2[i] {
			return false
		}
	}
	return true
}

// SameIntegerSlices returns true if two slices are the same
// Wendu
func SameIntegerSlices(s1, s2 []int) bool {
	if len(s1) != len(s2) {
		return false
	}
	for i := range s1 {
		if s1[i] != s2[i] {
			return false
		}
	}
	return true
}

// If there are two connected nodes in the graph without a divergence, merge the two nodes.
// ChainMerging
// Input: Graph
// Output: Graph
// Wenduo
func (dbGraph Graph) ChainMerging() Graph {
	var newGraph Graph
	newNodes := make(map[string]*Node)
	newEdges := make(map[string]*Edge)

	visited := make(map[string]bool)
	for _, node := range dbGraph.nodes {
		visited[node.label] = false
	}

	node := dbGraph.root // choose a starting node
	kmerLength := len(node.label) + 1
	toBeMerged := make([]*Node, 0)
	// make a stack
	stack := make([]*Node, 0)
	stack = append(stack, node)
	// run a DFS
	for len(stack) != 0 {
		// pop an element out of stack

		node = stack[len(stack)-1]
		stack = stack[:(len(stack) - 1)]

		visited[node.label] = true
		if len(node.children) <= 1 {
			toBeMerged = append(toBeMerged, node)
		} else { // have two or more children
			// add node to newGraph.nodes
			newNodes[node.label] = node
			// add edge to newGraph.edges
			for _, child := range node.children {
				kmer := node.label + child.label[kmerLength-2:kmerLength-1] // kmer is node.label+last character of child.label
				newEdges[kmer] = dbGraph.edges[kmer]
			}
			// Merge nodes in toBeMerged and add to graph
			if len(toBeMerged) >= 1 {
				var mergedNode Node
				newLabel := toBeMerged[0].label
				for i := 1; i < len(toBeMerged); i++ {
					newLabel += toBeMerged[i].label[kmerLength-2 : kmerLength-1]
				}
				mergedNode.label = newLabel
				mergedNode.children = toBeMerged[len(toBeMerged)-1].children
				mergedNode.inDegree = toBeMerged[0].inDegree
				mergedNode.outDegree = toBeMerged[len(toBeMerged)-1].outDegree
				newNodes[mergedNode.label] = &mergedNode
			}

			toBeMerged = make([]*Node, 0)

		}
		for _, nextNode := range node.children {
			if !visited[nextNode.label] {
				stack = append(stack, nextNode)
			}
		}

	}
	newGraph.nodes = newNodes
	newGraph.edges = newEdges
	
	for key, _ := range newGraph.nodes {
		newGraph.root = newGraph.nodes[key]
		break
	}
	
	return newGraph
}

//EulerianPath find the Eulerian path for the De brujin graph
//Input: the graph object, representing the built de brujin graph
//Output: the string list represent the eulerian path for kemrs, if no Eulerian path return empty string list
//Lilin
func EulerianPath(graph Graph) []string {
	//find the start node
	var start *Node
	//count the total number of edges
	var numEdge int

	var numOddVertices int
	for _, n := range graph.nodes {
		if n.inDegree < n.outDegree {
			start = n
		}
		if (n.inDegree+n.outDegree)%2 != 0 {
			numOddVertices += 1
		}
		numEdge += n.inDegree
	}
	//If a graph has more than two vertices of odd degree, no eulerian path
	if numOddVertices > 2 {
		var emptyString []string
		//return empty string list
		return emptyString
	}

	var genome []string
	genome = append(genome, start.label)

	var currentnode *Node
	currentnode = start

	for i := 0; i < numEdge; i++ {
		prefix := currentnode.children[len(currentnode.children)-1].label
		currentnode.children = currentnode.children[:len(currentnode.children)-1]
		currentnode = graph.nodes[prefix]
		genome = append(genome, currentnode.label)
	}

	return genome
}

//N50 takes the list of contigs and return the shortest contig length needed to cover 50% genome,
// describing the completeness of genome assembly
//Input: the list of contigs
//Output: The length of the shortest contig cover 50% genome
//Lilin
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

//Sort sorts the contigs according to their length from longest to shortest
//Input: the list of contigues
//Output: the list of the sorted contigs from longest to shortest
//Lilin
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

//DivideConquer sorts the list of int from largest to smallest using divide&conquer algorithm
//Intput: a list of integer
//Output: a sortest list of integer from largest to smallest
//Lilin
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
//DistinctKmerCount counts the number of distinct kmer in the string
//Input: a map of kmer and frequency
//Output: the number of distinct kmer
func DistinctKmerCount(kmerCounts map[string]int) int {

	count := len(kmerCounts)

	return count
}

//KmerSizeAdjustment counts the kmer size according to the readlength and kmer and genome coverage
//Input: the readlength, kmer coverage and genome coverage
//Output: the estimate optimal kmer size
func KmerSizeAdjustment(readlength int, kmer_coverage, genome_coverage float64) int {

	var kmer_size int
	kmer_size = 1 + readlength - int((kmer_coverage*float64(readlength))/genome_coverage)

	return kmer_size
}

//GenomeCoverage counts the genome coverage for the current kmer set
//Input: the current genome, the kmer set
//Output: the genome coverage
func GenomeCoverage(genome string, readlength, reads_number int) float64 {
	var genome_coverage float64

	genome_size := len(genome)

	//(A total number of reads * read length)/ (Estimated genome size)
	genome_coverage = float64(reads_number) * float64(readlength) / float64(genome_size)

	return genome_coverage
}

//KmerCoverage counts the kmer coverage for the current kmer set
//Input: the read length, kmer size and genome coverage
//Output: the kmer coverage
func KmerCoverage(readlength, kmersize int, genome_coverage float64) float64 {

	var kmer_coverage float64
	base_coverage := genome_coverage
	//𝐶𝑘=𝐶⋅(𝑅−𝐾+1)/𝑅
	kmer_coverage = base_coverage * float64(readlength-kmersize+1) / float64(readlength)

	return kmer_coverage
}

//KmerSizeSet takes the readlength and generate the set of interger for range over to pick the optimal kmer size form this set
//Input: the readlength
//Output: the generated kmer size set
func KmerSizeSet(readlength int) []int {
	var size_set []int
	//find the prome number less or equal to the readlength
	size_set = FindPrime(readlength)
	return size_set
}

//FindPrime takes in the value and returns the list of prime number smaller than that value
// in this case, although 1 is not prime number, we also take it into consideration
//Input: the value
//Output: the int list contains the prime number smaller than the value
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

//OptimalKmerSize gives the optiaml kmersize for the later analysis
//Input: genome, list of reads
//Output: optimal kmer size
func OptimalKmerSize(genome string, reads []string) int {

	var optimalk int

	readlength := len(reads[0])
	//method1: find the k gives most number of distinct kmers
	k_size_set := KmerSizeSet(readlength)
	max_distinct_kmer_count := 0
	var optimalk1 int
	for i := range k_size_set {
		k := k_size_set[i]
		kmerCounts := KmerHash(k, genome)
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

	readsnumber = len(reads)
	genome_coverage = GenomeCoverage(genome, readlength, readsnumber)
	kmer_coverage = KmerCoverage(readlength, optimalk1, genome_coverage)

	optimalk = KmerSizeAdjustment(readlength, kmer_coverage, genome_coverage)

	return optimalk

}

