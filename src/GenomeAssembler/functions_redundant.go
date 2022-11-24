package main

// KmerHash
// Input: An integer k and a string Text, and already exists kmer map.
// Output: a map of kmer and frequency
// Wenduo Lilin
func KmerHash0(k int, text string, kmerCounts map[string]int) map[string]int {
	//count how many times each kmer occurred in this sequence
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

func KmerHash1(kmers []string) map[string]int {
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

// String Composition Problem: Generate the k-mer composition of a string.
// Input: An integer k and a string Text.
// Output: Compositionk(Text), where the k-mers are arranged in lexicographic order.
// Wenduo
func StringComposition0(k int, text string) []string {
	//count how many times each kmer occurred in this sequence
	kmers := make(map[string]int)
	kmerCounts := KmerHash0(k, text, kmers)

	kmerComposition := make([]string, 0, len(kmerCounts))
	for kmer := range kmerCounts {
		kmerComposition = append(kmerComposition, kmer)
	}

	return kmerComposition
}

// De Bruijn Graph from a String Problem
// Construct the de Bruijn graph of a string.
// Input: An integer k and a string Text.
// Output: DeBruijn_k(Tejxt), in the form of an adjacency matrix/list.
// Wenduo; Lilin
func DeBruijnGraph0(k int, sequence string) map[string][]string {
	kmerComposition := StringComposition0(k, sequence)

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
func DeBruijnGraph1(k int, texts []string) Graph {
	var dbgraph Graph
	dbnodes := make(map[string]*Node)
	dbedges := make(map[string]*Edge)

	var kmerComposition []string
	for _, text := range texts {
		kmer_list := StringComposition0(k, text)
		kmerComposition = append(kmerComposition, kmer_list...)
	}

	kmers := make(map[string]int)
	for _, text := range texts {
		kmers = KmerHash0(k, text, kmers)
	}

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

// EulerianPath find the Eulerian path for the De brujin graph
// Input: the graph object, representing the built de brujin graph
// Output: the string list represent the eulerian path for kemrs, if no Eulerian path return empty string list
// Lilin
func EulerianPath0(graph Graph) [][]string {
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
