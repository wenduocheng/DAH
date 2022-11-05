package main

import "fmt"

func main() {
	k := 5
	text := "CAATCCAAC"
	fmt.Println(StringComposition(k, text))
	fmt.Println(KmerHash(k, text))

	sequencePath := make([]string, 5)
	sequencePath[0] = "ACCGA"
	sequencePath[1] = "CCGAA"
	sequencePath[2] = "CGAAG"
	sequencePath[3] = "GAAGC"
	sequencePath[4] = "AAGCT"
	fmt.Println(ReconstructStringFromGenomePath(sequencePath))

	fmt.Println(OverlapGraph(sequencePath))

	text2 := "AAGATTCTCTACAAGA"
	fmt.Println("DeBrujinGraph was created")
	fmt.Println(DeBruijnGraph(4, text2))

	graph := DeBruijnGraph2(4, text2)
	fmt.Println("DeBrujinGraph2 was created")
	for i := range graph.nodes {
		fmt.Println(graph.nodes[i].label, ":")
		for j := range graph.nodes[i].children {
			fmt.Println(graph.nodes[i].children[j].label)
		}
	}
	for i := range graph.edges {
		fmt.Println(graph.edges[i].label, graph.edges[i].from.label, graph.edges[i].to.label, graph.edges[i].weight)
	}

}

// String Composition Problem: Generate the k-mer composition of a string.
// Input: An integer k and a string Text.
// Output: Compositionk(Text), where the k-mers are arranged in lexicographic order.
func StringComposition(k int, text string) []string {
	//count how many times each kmer occurred in this sequence
	kmers := KmerHash(k, text)

	kmerComposition := make([]string, 0, len(kmers))
	for kmer := range kmers {
		kmerComposition = append(kmerComposition, kmer)
	}

	return kmerComposition
}

// KmerHash
// Input: An integer k and a string Text.
// Output: a map of kmer and frequency
func KmerHash(k int, text string) map[string]int {
	//count how many times each kmer occurred in this sequence
	kmers := make(map[string]int)
	for j := 0; j <= len(text)-k; j++ {
		kmer := text[j : j+k]
		_, exists := kmers[kmer]
		if exists {
			kmers[kmer]++
		} else {
			kmers[kmer] = 1
		}
	}

	return kmers
}

// String Spelled by a Genome Path Problem. Reconstruct a string from its genome path.
// Input: A sequence path of k-mers Pattern1, … ,Patternn such that the last k - 1 symbols of Patterni are equal to the first k-1 symbols of Patterni+1 for 1 ≤ i ≤ n-1.
// Output: A string Text of length k+n-1 such that the i-th k-mer in Text is equal to Patterni (for 1 ≤ i ≤ n).
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
func Suffix(kmer string) string {
	if len(kmer) <= 1 {
		panic("Error: kmer size too small")
	}
	return kmer[1:len(kmer)]
}

// Prefix
// Given: A string respresenting a kmer
// Return: A string corresponding to the prefix of the kmer
func Prefix(kmer string) string {
	if len(kmer) <= 1 {
		panic("Error: kmer size too small")
	}
	return kmer[0:(len(kmer) - 1)]
}

type Graph struct {
	nodes map[string]*Node //k-1 mer
	edges map[string]*Edge //kmer

}

type Node struct {
	label     string
	inDegree  int
	outDegree int
	children  []*Node
}

type Edge struct {
	label  string
	from   *Node
	to     *Node
	weight int
}

// De Bruijn Graph from a String Problem
// Construct the de Bruijn graph of a string.
// Input: An integer k and a string Text.
// Output: DeBruijn_k(Tejxt), in the form of an adjacency matrix/list.
func DeBruijnGraph(k int, text string) map[string][]string {
	kmerComposition := StringComposition(k, text)

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
	return dbgraph
}

