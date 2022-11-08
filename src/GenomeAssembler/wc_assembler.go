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
// Wenduo
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
// Wenduo
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

func EulerianCycle(graph Graph) []string {
	//find the start node
	var start *Node
	//count the total number of edges
	var numEdge int
	for _, n := range graph.nodes {
		if n.inDegree < n.outDegree {
			start = n
		}
		numEdge += n.inDegree
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

func N50(reads []string) int {
	var total_length int
	reads = Sort(reads)
	for _, read := range reads {
		total_length += len(read)
	}

	N50_length := int(total_length / 2)
	var currLength int
	var N50_string string
	for _, read := range reads {
		currLength += len(read)
		if currLength >= N50_length {
			N50_string = read
			break
		}
	}

	return len(N50_string)
}
func Sort(reads []string) []string {
	lenHash := make(map[int][]string)
	for _, read := range reads {
		length := len(read)
		_, Exists := lenHash[length]
		if Exists {
			lenHash[length] = append(lenHash[length], read)
		} else {
			lenHash[length] = []string{read}
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

// apply divide and conquer
func DivideConquer(numlist []int) []int {
	if len(numlist) == 1 || len(numlist) == 0 {
		return numlist
	}
	i := int(len(numlist) / 2)
	right := DivideConquer(numlist[:i])
	left := DivideConquer(numlist[i:])
	var combineList []int
	for len(right) != 0 && len(left) != 0 {
		if right[0] <= left[0] {
			combineList = append(combineList, right[0])
			right = right[1:]
		}
		if right[0] > left[0] {
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
	dbgraph.root = dbgraph.nodes[kmerComposition[0]]
	return dbgraph
}


// GenerateSequence randomly genertes a sequence composed of A/T/C/G
// Wenduo
func GenerateSequence(length int) string {
	nucleotides := []string{"A", "T", "C", "G"}
	sequence := ""
	for i := 0; i < length; i++ {
		index := rand.Intn(4)
		sequence += nucleotides[index]
	}
	return sequence
}

// RandomMutate randomly mutates a slice of kmers
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
func kmerCountSort(kmerUniqCounts map[int]int) []int {
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

	visited := make(map[string]bool)
	for _, node := range dbGraph.nodes {
		visited[node.label] = false
	}
	node := dbGraph.root // choose a starting node
	kmerLength := len(node.label) + 1
	toBeMerged := make([]*Node, 0)
	stack := make([]*Node, 0)
	stack = append(stack, node)
	// run a DFS
	for len(stack) != 0 {
		node = stack[len(stack)-1]
		stack = stack[:(len(stack) - 1)]
		if len(node.children) <= 1 {
			toBeMerged = append(toBeMerged, node)
		} else { // have two or more children
			// add node to newGraph.nodes
			newGraph.nodes[node.label] = node
			// add edge to newGraph.edges
			for _, child := range node.children {
				kmer := node.label + child.label[kmerLength-1:kmerLength] // kmer is node.label+last character of child.label
				newGraph.edges[kmer] = dbGraph.edges[kmer]
			}
			// Merge nodes in toBeMerged and add to graph
			// Merge(toBeMerged) // merge the nodes into one
			// // add to new graph
			var mergedNode Node
			newLabel := toBeMerged[0].label
			for i := 1; i < len(toBeMerged); i++ {
				newLabel += toBeMerged[i].label[kmerLength-1 : kmerLength]
			}
			mergedNode.label = newLabel
			mergedNode.children = toBeMerged[len(toBeMerged)].children
			mergedNode.inDegree = toBeMerged[0].inDegree
			mergedNode.outDegree = toBeMerged[len(toBeMerged)].outDegree
			newGraph.nodes[mergedNode.label] = &mergedNode

			toBeMerged = make([]*Node, 0)

		}
		for _, nextNode := range node.children {
			if !visited[nextNode.label] {
				stack = append(stack, nextNode)
			}
		}
	}
	return newGraph
}
