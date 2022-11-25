package main

import (
	"sort"
)

// Nomenclature
// kmer: a single kmer
// kmerLength: length of kmer
// kmers: a list of kmers
// read: a single read
// kmerCounts: a [string]int map of kmers; key is kmer and value is frequency
// dbGraph: a de Bruijn graph

// DenovoAssembler is the highest level function of our assembler
// Input: a list of strings corresponding to sequencing reads; an integer corresponding to kmer size (to be modified)
// Output: a list of strings corresponding to output contigs
// Wenduo; Lilin; Tianyue
func DenovoAssembler(reads []string, kmerLength int) [][]string {
	// First step: Determine kmer size

	// Second step: Hash the reads
	// kmerCounts := KmerHashFromReads(kmerLength, reads)
	// Generate a Kmer Distribution Plot
	// uniqueKmerCounts := GetUniqueKmerCounts(kmerCounts)
	// sortedUniqCounts := KmerCountSort(uniqueKmerCounts)
	// DrawKmerScatter(sortedUniqCounts)

	// Third step: Construct the de Bruijn graph
	dbGraph := DeBruijnGraph(kmerLength, reads)

	// Fourth step: Simplify the de Bruijn graph
	mergedGraph := dbGraph.ChainMerging()

	// Fifth step: Output the contigs
	contigsPath := EulerianPath(mergedGraph)
	contigs := AssembleContigs(contigsPath)

	return contigs
}

// kmerHashFromReads hashes a list of reads
// Input: an integer corresponding to kmerLength, a list of string corresponding to sequencing reads
// Output: a map[string] int; the key is kmer(string) and the value is frequency(int)
// Wenduo
func KmerHashFromReads(kmerLength int, reads []string) map[string]int {
	//count how many times each kmer occurred in a list of reads
	kmerCounts := make(map[string]int)

	for _, read := range reads {
		for j := 0; j <= len(read)-kmerLength; j++ {
			kmer := read[j : j+kmerLength]
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

// GetKmerComposition gets the kmer composition from the kmerCounts map
// Adapted from the String Composition Problem: Generate the k-mer composition of a string.
// Input: a map corresponding to kmerCounts.
// Output: a list of strings corresponding to kmers
// Wenduo
func GetKmerComposition(kmerCounts map[string]int) []string {
	kmerComposition := make([]string, 0, len(kmerCounts))
	for kmer := range kmerCounts {
		kmerComposition = append(kmerComposition, kmer)
	}

	return kmerComposition
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
// Input: An integer k and a list of reads.
// Output: a graph
// Wenduo; Lilin
func DeBruijnGraph(kmerLength int, reads []string) Graph {
	var dbGraph Graph
	dbnodes := make(map[string]*Node)
	dbedges := make(map[string]*Edge)

	kmerCounts := KmerHashFromReads(kmerLength, reads)

	kmerComposition := GetKmerComposition(kmerCounts)

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
		newEdge.weight = kmerCounts[kmer]
		newEdge.label = kmer

		dbedges[kmer] = &newEdge

		// add sufix node as a child to prefix node
		dbnodes[Prefix(kmer)].children = append(dbnodes[Prefix(kmer)].children, dbnodes[Suffix(kmer)])

	}
	dbGraph.nodes = dbnodes
	dbGraph.edges = dbedges
	
	lowestInDegree := len(dbGraph.nodes)
	var rootNode *Node
	for key, _ := range dbGraph.nodes {
		if dbGraph.nodes[key].inDegree < lowestInDegree {
			lowestInDegree = dbGraph.nodes[key].inDegree
			rootNode = dbGraph.nodes[key]
		}
	}
	dbGraph.root = rootNode
	
	return dbGraph
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

	// Merge nodes in toBeMerged after running the DFS
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

	newGraph.nodes = newNodes
	newGraph.edges = newEdges

	for key, _ := range newGraph.nodes {
		newGraph.root = newGraph.nodes[key]
		break
	}

	return newGraph
}

// EulerianPath find the Eulerian path for the De brujin graph
// Input: the graph object, representing the built de brujin graph
// Output: the string list represent the eulerian path for kemrs, if no Eulerian path return empty string list
// Lilin
func EulerianPath(graph Graph) [][]string {
	//find the start node
	var start []*Node
	for _, n := range graph.nodes {
		if n.inDegree < n.outDegree {
			start = append(start, n)
		}

	}

	var contigs [][]string
	for i, _ := range start {
		startNode := start[i]
		currentnode := startNode
		contigs[i] = append(contigs[i], startNode.label)
		//if start node has no more children nodes,continue
		for len(startNode.children) != 0 {
			prefix := currentnode.children[len(currentnode.children)-1].label
			currentnode.children = currentnode.children[:len(currentnode.children)-1]
			currentnode = graph.nodes[prefix]
			contigs[i] = append(contigs[i], currentnode.label)
		}
	}
	return contigs
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

// SameStringSlices returns true if two slices are the same
// This function is for test
// Wenduo
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
// This function is for test
// Wenduo
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

//AssembleContigs takes in path of different contigs and assemble them
//Input: kmer path of different contigs
//Output: the list of assembled kmers
//Lilin
func AssembleContigs(contigs [][]string) []string {
	var contigs_list []string
	num_contigs := len(contigs)
	contigs_list = make([]string, num_contigs)
	for i, _ := range contigs {
		contigs_list[i] = ReconstructStringFromGenomePath(contigs[i])
	}
	return contigs_list
}
