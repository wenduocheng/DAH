package main

import (
	"bufio"
	"fmt"
	"os"
	"sort"
	"strconv"
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
func DenovoAssembler(reads []string, kmerLength int) []string {
	// First step: Determine kmer size

	// Second step: Hash the reads
	//kmerCounts := KmerHashFromReads(kmerLength, reads)
	// Generate a Kmer Frequency Distribution Plot
	// uniqueKmerCounts := GetUniqueKmerCounts(kmerCounts)
	// sortedUniqCounts := KmerCountSort(uniqueKmerCounts)

	// DrawHistogram(sortedUniqCounts)

	// Third step: Construct the de Bruijn graph
	dbGraph := DeBruijnGraph(kmerLength, reads)
	// SaveGraphToGFA(dbGraph, "deBruijnGraph")
	// Fourth step: Simplify the de Bruijn graph
	mergedGraph := dbGraph.ChainMerging()
	// SaveGraphToGFA(mergedGraph, "merged")

	tipclipedGraph := mergedGraph.TipClip(kmerLength)
	// SaveGraphToGFA(tipclipedGraph, "tipclip")

	// Fifth step: Output the contigs
	contigsPath := EulerianPath(tipclipedGraph)
	contigs := AssembleContigs(contigsPath, kmerLength)

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
		// add prefix node as a parent to suffix node
		dbnodes[Suffix(kmer)].parents = append(dbnodes[Suffix(kmer)].parents, dbnodes[Prefix(kmer)])
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
// Input: a string representing a sequence and an integer representing the kmer length
// Output: a list of kmers
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

		// add the node to toBeMerged if it satifies the conditions
		if len(node.children) == 1 && !visited[node.children[0].label] && (len(node.parents) == 1 || node.parents == nil) {
			toBeMerged = append(toBeMerged, node)
		} else if node.children == nil {
			// if the node has no children, clear out the toBeMerged
			toBeMerged = append(toBeMerged, node)
			var mergedNode Node
			newLabel := toBeMerged[0].label
			for i := 1; i < len(toBeMerged); i++ {
				newLabel += toBeMerged[i].label[kmerLength-2 : kmerLength-1]
			}
			// Create the merged node
			mergedNode.label = newLabel
			mergedNode.inDegree = toBeMerged[0].inDegree
			mergedNode.outDegree = toBeMerged[len(toBeMerged)-1].outDegree
			newNodes[mergedNode.label] = &mergedNode

			for _, parent := range toBeMerged[0].parents {
				mergedNode.parents = append(mergedNode.parents, newNodes[parent.label])
				newNodes[parent.label].children = append(newNodes[parent.label].children, &mergedNode)
				// Create the edge
				var newEdge Edge
				newEdge.label = parent.label + mergedNode.label[kmerLength-2:len(mergedNode.label)]
				newEdge.from = parent
				newEdge.to = &mergedNode
				newEdge.weight = dbGraph.edges[parent.label+toBeMerged[0].label[kmerLength-2:kmerLength-1]].weight
				newEdges[newEdge.label] = &newEdge
			}
			// Reset the toBeMerged
			toBeMerged = make([]*Node, 0)

		} else if len(node.children) == 1 && visited[node.children[0].label] {
			// if a node has a child that is itself
			// this will generate a self-pointing ring
			if node.label == node.children[0].label {
				var newEdge Edge
				newEdge.label = node.label + node.children[0].label[kmerLength-2:kmerLength-1]
				newEdge.from = newNodes[node.label]
				newEdge.to = newNodes[node.label]
				newEdge.weight = dbGraph.edges[newEdge.label].weight
				newEdges[newEdge.label] = &newEdge

				newNodes[node.label].children = append(newNodes[node.label].children, newNodes[node.label])
				newNodes[node.label].parents = append(newNodes[node.label].parents, newNodes[node.label])
			} else {
				toBeMerged = append(toBeMerged, node)
			}
			if len(toBeMerged) >= 1 {
				var mergedNode Node
				newLabel := toBeMerged[0].label
				for i := 1; i < len(toBeMerged); i++ {
					newLabel += toBeMerged[i].label[kmerLength-2 : kmerLength-1]
				}
				mergedNode.label = newLabel
				mergedNode.inDegree = toBeMerged[0].inDegree
				mergedNode.outDegree = toBeMerged[len(toBeMerged)-1].outDegree
				newNodes[mergedNode.label] = &mergedNode

				for _, parent := range toBeMerged[0].parents {
					mergedNode.parents = append(mergedNode.parents, newNodes[parent.label])
					newNodes[parent.label].children = append(newNodes[parent.label].children, &mergedNode)
					var newEdge Edge
					newEdge.label = parent.label + mergedNode.label[kmerLength-2:len(mergedNode.label)]
					newEdge.from = parent
					newEdge.to = &mergedNode
					newEdge.weight = dbGraph.edges[parent.label+toBeMerged[0].label[kmerLength-2:kmerLength-1]].weight
					newEdges[newEdge.label] = &newEdge
				}

				mergedNode.children = append(mergedNode.children, newNodes[node.children[0].label])
				newNodes[node.children[0].label].parents = append(newNodes[node.children[0].label].parents, &mergedNode)
				var newEdge Edge
				newEdge.label = mergedNode.label + node.children[0].label[kmerLength-2:kmerLength-1]
				newEdge.from = &mergedNode
				newEdge.to = newNodes[node.children[0].label]
				newEdge.weight = dbGraph.edges[toBeMerged[len(toBeMerged)-1].label+node.children[0].label[kmerLength-2:kmerLength-1]].weight
				newEdges[newEdge.label] = &newEdge
			}
			toBeMerged = make([]*Node, 0)
		} else {
			// if a node has two or more children
			// add node to newGraph.nodes
			// newNodes[node.label] = node

			var newNode Node
			newNode.label = node.label
			newNode.inDegree = node.inDegree
			newNode.outDegree = node.outDegree
			newNodes[newNode.label] = &newNode

			for _, child := range node.children {
				// if a node has a child that is the node itself
				if child.label == node.label {
					var newEdge Edge
					newEdge.label = node.label + node.children[0].label[kmerLength-2:kmerLength-1]
					newEdge.from = newNodes[node.label]
					newEdge.to = newNodes[node.label]
					newEdge.weight = dbGraph.edges[newEdge.label].weight
					newEdges[newEdge.label] = &newEdge
					// newNodes[node.label].children = append(newNodes[node.label].children, newNodes[node.label])
					// newNodes[node.label].parents = append(newNodes[node.label].parents, newNodes[node.label])
				}
			}

			// check if the parent nodes of the node are in the newNodes
			for _, parent := range node.parents {
				_, exists := newNodes[parent.label]
				if exists {
					newNodes[parent.label].children = append(newNodes[parent.label].children, newNodes[node.label])
					newNodes[node.label].parents = append(newNodes[node.label].parents, newNodes[parent.label])
					var newEdge Edge
					newEdge.label = newNodes[parent.label].label + newNodes[node.label].label[kmerLength-2:kmerLength-1]
					newEdge.from = newNodes[parent.label]
					newEdge.to = newNodes[node.label]
					newEdge.weight = dbGraph.edges[newEdge.label].weight
					newEdges[newEdge.label] = &newEdge
				}
			}

			// check if the children nodes of the node are in the newNodes
			for _, child := range node.children {
				_, exists := newNodes[child.label]
				if exists && !newNodes[child.label].InChildren(newNodes[node.label]) {
					newNodes[node.label].children = append(newNodes[node.label].children, newNodes[child.label])
					newNodes[child.label].parents = append(newNodes[child.label].parents, newNodes[node.label])
					var newEdge Edge
					newEdge.label = newNodes[node.label].label + newNodes[child.label].label[kmerLength-2:kmerLength-1]
					newEdge.from = newNodes[node.label]
					newEdge.to = newNodes[child.label]
					newEdge.weight = dbGraph.edges[newEdge.label].weight
					newEdges[newEdge.label] = &newEdge
				}
			}

			// Merge the nodes in toBeMerged into one node and add it to graph
			if len(toBeMerged) >= 1 {
				var mergedNode Node
				newLabel := toBeMerged[0].label
				// Get the label of the merged node and the weight sum
				for i := 1; i < len(toBeMerged); i++ {
					newLabel += toBeMerged[i].label[kmerLength-2 : kmerLength-1]
				}
				mergedNode.label = newLabel
				// mergedNode.children = toBeMerged[len(toBeMerged)-1].children
				mergedNode.inDegree = toBeMerged[0].inDegree
				mergedNode.outDegree = toBeMerged[len(toBeMerged)-1].outDegree
				newNodes[mergedNode.label] = &mergedNode

				// mergedNode parents
				for _, parent := range toBeMerged[0].parents {
					mergedNode.parents = append(mergedNode.parents, newNodes[parent.label])
					newNodes[parent.label].children = append(newNodes[parent.label].children, &mergedNode)

					var newEdge Edge
					newEdge.label = parent.label + mergedNode.label[kmerLength-2:len(mergedNode.label)]
					newEdge.from = parent
					newEdge.to = &mergedNode
					newEdge.weight = dbGraph.edges[parent.label+toBeMerged[0].label[kmerLength-2:kmerLength-1]].weight
					newEdges[newEdge.label] = &newEdge
				}

				// Add the mergedNode as the parent node of current node
				newNode.parents = append(newNode.parents, &mergedNode)
				// Add the currentNode as the children node of merged node
				mergedNode.children = append(mergedNode.children, &newNode)

				// Add the edge between the current node and the merged node
				var newEdge Edge
				newEdge.label = mergedNode.label + node.label[kmerLength-2:kmerLength-1]
				newEdge.from = &mergedNode
				newEdge.to = node
				newEdge.weight = dbGraph.edges[toBeMerged[len(toBeMerged)-1].label+node.label[kmerLength-2:kmerLength-1]].weight
				newEdges[newEdge.label] = &newEdge
			}
			toBeMerged = make([]*Node, 0)
		}
		visited[node.label] = true

		// Add unvisited children nodes to the stack
		for _, nextNode := range node.children {
			if !visited[nextNode.label] {
				stack = append(stack, nextNode)
			}
		}

	}

	newGraph.nodes = newNodes
	newGraph.edges = newEdges

	// Set the root
	for key, _ := range newGraph.nodes {
		newGraph.root = newGraph.nodes[key]
		break
	}

	return newGraph
}

// InChildren is a method of *Node to check if node is in parent.children
// Wenduo
func (node *Node) InChildren(parent *Node) bool {
	for _, child := range parent.children {
		if child.label == node.label {
			return true
		}
	}
	return false
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
	var contigs_path [][]string
	//special case:
	//when results of the ChainMergeing is only one complete genome no indegree
	//and outdegree of the node
	if len(start) == 0 {
		first := make([]string, 1)
		first[0] = graph.root.label
		contigs_path = append(contigs_path, first)
		return contigs_path
	}

	for i, _ := range start {
		startNode := start[i]
		currentnode := startNode
		//index := len(contigs_path)
		first := make([]string, 1)
		contigs_path = append(contigs_path, first)
		index := len(contigs_path) - 1
		contigs_path[index][0] = startNode.label
		//if start node has no more children nodes,continue
		for len(startNode.children) != 0 {
			prefix := currentnode.children[len(currentnode.children)-1].label
			currentnode.children = currentnode.children[:len(currentnode.children)-1]
			currentnode = graph.nodes[prefix]
			if currentnode == nil {
				break
			}
			contigs_path[index] = append(contigs_path[index], currentnode.label)
			startNode = currentnode

		}
	}

	//to check if finish all edges if not continue
	a := true
	for a {
		var start_node *Node
		for _, n := range graph.nodes {
			a = false
			if len(n.children) > 0 {
				start_node = n
				a = true
				break
			}
		}
		if a == false {
			break
		}
		currentnode := start_node
		//index := len(contigs_path)
		first := make([]string, 1)
		contigs_path = append(contigs_path, first)
		index := len(contigs_path) - 1
		contigs_path[index][0] = start_node.label
		//if start node has no more children nodes,continue
		for len(start_node.children) != 0 {
			prefix := currentnode.children[len(currentnode.children)-1].label
			currentnode.children = currentnode.children[:len(currentnode.children)-1]
			currentnode = graph.nodes[prefix]
			if currentnode == nil {
				break
			}
			contigs_path[index] = append(contigs_path[index], currentnode.label)
			start_node = currentnode

		}
	}

	return contigs_path
}

// String Spelled by a Genome Path Problem. Reconstruct a string from its genome path.
// Input: A sequence path of k-mers Pattern1, … ,Patternn such that the last k - 1 symbols of Patterni are equal to the first k-1 symbols of Patterni+1 for 1 ≤ i ≤ n-1.Kmer_length.
// Output: A string Text of length k+n-1 such that the i-th k-mer in Text is equal to Patterni (for 1 ≤ i ≤ n).
// Wenduo; Lilin
func ReconstructStringFromGenomePath(sequencePath []string, kmer_length int) string {
	genome := sequencePath[0]
	start_point := kmer_length - 2

	for i := 1; i < len(sequencePath); i++ {
		k := len(sequencePath[i])
		genome = genome + sequencePath[i][start_point:k]
	}
	return genome
}

// Contain returns true if an element is in the slice
// This function is for test
// Wenduo
func Contain(kmerComposition []string, kmer string) bool {
	for _, val := range kmerComposition {
		if val == kmer {
			return true
		}
	}
	return false
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

// AssembleContigs takes in path of different contigs and assemble them
// Input: kmer path of different contigs, kmerlength
// Output: the list of assembled kmers
// Lilin
func AssembleContigs(contigs [][]string, kmer_length int) []string {
	var contigs_list []string
	num_contigs := len(contigs)
	contigs_list = make([]string, num_contigs)
	for i, _ := range contigs {
		contigs_list[i] = ReconstructStringFromGenomePath(contigs[i], kmer_length)
	}
	return contigs_list
}

// CalculateGCConent determines the GC content of a sequence string
// Wenduo
func CalculateGCConent(sequence string) int {
	count := 0
	for i := 0; i < len(sequence); i++ {
		if sequence[i:i+1] == "c" || sequence[i:i+1] == "g" || sequence[i:i+1] == "C" || sequence[i:i+1] == "G" {
			count++
		}
	}
	return count
}

// //Tip Clip will iterate through the de Bruijn Graph and find the dead end part that is
// //less 2kmer length
// input: de Bruijn Graph after chain merging
// output: a new de Bruijn Graph after tip clipped
// tianyue
func (graph Graph) TipClip(kmerLength int) Graph {

	newGraph := graph.CopyGraph()
	var currentNode *Node
	var pastNode *Node
	//starting from the root of the graph
	currentNode = newGraph.root
	//fmt.Println(currentNode, newGraph.root)
	// fmt.Println("old root", graph.root)

	//iterate the nodes and find the node with more than one child
	//keey track of how many times each node will be visit
	visited := make(map[string]int)
	//divergenodes will keey track of the pointes where diverge happened
	var divergenodes []*Node
	// divergenodes = append(divergenodes, currentNode)
	// map path is used to track the starting diverge point and the end point of the current bunch. a string is matching to a list of string means that a diverge node can have more than 1 child
	path := make(map[string][]string)
	//path weight used to calculate the weight of this divergent path, and the key is storing the start node lable + end node label
	// var pathweight map[string]int
	//iterate the graph, and check the nodes where key represents k-1 mer
	for currentNode != nil {
		//first update visited, which is a map to keep track the nodes that are visiting and in the end, it will outoput the
		_, isExit := visited[currentNode.label]
		if !isExit {
			visited[currentNode.label] = 1
		} else {
			visited[currentNode.label] += 1
		}
		//second, going through the paths in the graph,
		if len(currentNode.children) > 1 {
			// 			//fmt.Println("currentNode.children >1")
			//this means that the current divergence node is downstream of another divergenodes, and since I need to iterate the path map later, i also want to indicate the linear relation between divergennodes. If i want to check if a given divergenode is the very last/tip one, just iterate the whole map with the key as the last string attached to the "currentKey"
			// fmt.Println("current Node", currentNode, "first", currentNode.children)
			divergenodes = append(divergenodes, currentNode)
			path[currentNode.label] = append(path[currentNode.label], currentNode.label)
			currentNode.inDegree -= 1
			pastNode = currentNode
			currentNode = currentNode.children[0]
			//update the path of the path

			// path[divergenodes[len(divergenodes)-1].label] = append(path[divergenodes[len(divergenodes)-1].label], currentNode.label)

			newGraph.nodes[pastNode.label].children = pastNode.children[1:]
			// pastNode.outDegree--
			newGraph.nodes[pastNode.label].outDegree -= 1
		} else if len(currentNode.children) == 1 {
			// fmt.Println("currentNode.child =1")
			// fmt.Println("current Node", currentNode, "first", currentNode.children)
			// pastNode.outDegree = 0
			// fmt.Println(currentNode)
			currentNode.inDegree -= 1
			currentNode.outDegree -= 1
			if len(divergenodes) != 0 {
				path[divergenodes[len(divergenodes)-1].label] = append(path[divergenodes[len(divergenodes)-1].label], currentNode.label)
				divergenodes = divergenodes[:len(divergenodes)-1]
			}
			pastNode = currentNode
			currentNode = currentNode.children[0]
			newGraph.nodes[pastNode.label].children = nil

		} else {
			pastNode = currentNode
			pastNode.children = nil

			if len(divergenodes) != 0 {
				currentNode = divergenodes[len(divergenodes)-1]

			} else {
				currentNode = nil
			}
		}
	}
	// fmt.Println("path", path)
	// fmt.Println("divergeNODE", divergenodes)
	//generate a list that stores the lable of nodes that need to be deleted
	//not considering the weight of the tips
	deleteEdge := make(map[string][]string)
	for _, val := range path {
		var pathsequence string
		for i := 0; i < len(val); i++ {

			if i == 0 {
				pathsequence = val[i]
			} else {
				pathsequence += val[i][2:]
			}
		}
		//if the sequence in this tip is less than 2 kmer length, and the last node lable is not recorded as another divergence node, which means that this tip is indeeded not integrated into the graph
		if len(pathsequence) <= kmerLength+1 && path[val[len(val)-1]] == nil {
			//the last one need to be add specially
			deleteEdge[val[0]] = append(deleteEdge[val[0]], path[val[0]][1:]...)
		}
	}
	//print out the edges that need to be delet if needed
	// fmt.Println("edges need to be deleted", deleteEdge)
	//final step, check the graph (nodes and edges) and delete the nodes in the map deleteEdge and the edges that is connecting to it
	//this for loop is to check nodes
	for key, val := range deleteEdge {
		// 	//if there are only one child or not child attached to this key
		if len(graph.nodes[key].children) <= 1 {
			panic("should not include this node")
		} else {
			for i := 0; i < len(val); i++ {
				graph.nodes[val[i]] = nil
				delete(graph.nodes, val[i])
				// fmt.Println("this node should be deleted", graph.nodes[val[i]])
			}
			//find the kmer of the edge
			for i := 0; i < len(val)-1; i++ {
				kmer := val[i][:1] + val[i+1][1:]
				graph.edges[kmer] = nil
				delete(graph.edges, kmer)
				// fmt.Println("this edge should be deleted", graph.edges[kmer], kmer)
			}
		}
	}
	return graph
}

// //Copy graph will take a Graph and make a copy of it
// //input: De bruijn Graph
// //output: a copy of a de Bruijn graph
func (graph Graph) CopyGraph() Graph {
	// fmt.Println("copy graph old graph root", graph.root)
	var newGraph Graph
	newNodes := make(map[string]*Node, len(graph.nodes))
	newEdges := make(map[string]*Edge, len(graph.edges))
	var newRoot Node
	newGraph.edges = newEdges
	// fmt.Println("check if map edges have the same memory address", "old:", &graph.edges, "new:", &newGraph.edges, "and their len, old:", len(graph.edges), "new:", len(newGraph.edges))
	newGraph.nodes = newNodes
	// fmt.Println("check if map nodes have the same memory address", "old:", &graph.nodes, "new:", &newGraph.nodes, "and their len, old:", len(graph.nodes), "new:", len(newGraph.nodes))
	newGraph.root = &newRoot

	//copying all the value of root into the new copy of root
	newRoot.label = graph.root.label
	var newChildren []*Node
	newGraph.root.children = newChildren
	// var newParents []*Node
	newRoot.inDegree = graph.root.inDegree
	newRoot.outDegree = graph.root.outDegree
	for key, _ := range graph.nodes {
		// fmt.Println("start to iterate graph.nodes")
		//fmt.Println(key, val)
		_, isExist := newGraph.nodes[key]
		// fmt.Println("get to here", isExist)
		if isExist {
			// newGraph.nodes[key].outDegree+=1
			//If exist, there must be something wrong since I am doing a direct copy, so there indegree/outdegree should be the most updatedvales, thus if this condition is met, it measn there must be something wrong

			// fmt.Println("finished")
		} else {
			//deal with the list of children at the very last, first make sure to make all the nodes of the new graph, and update the easy value as int/string
			// fmt.Println("else")
			var newNode Node
			var pointertoNewNode *Node
			pointertoNewNode = &newNode
			newNode.label = key
			// fmt.Println(key, "key")
			newNode.inDegree = graph.nodes[key].inDegree
			newNode.outDegree = graph.nodes[key].outDegree
			newGraph.nodes[newNode.label] = pointertoNewNode
			//for this node, copy the node relation
		}
	}
	// fmt.Print("finished ranging nodes")

	if len(newGraph.nodes) != len(graph.nodes) {
		panic("Error: the len of map nodes are not the same length")
	}
	for i := 0; i < len(graph.root.children); i++ {
		newGraph.root.children = append(newGraph.root.children, newGraph.nodes[graph.root.children[i].label])
	}

	for key, val := range graph.nodes {
		var newChild []*Node
		for i := 0; i < len(val.children); i++ {
			newChild = append(newChild, newGraph.nodes[graph.nodes[key].children[i].label])
			newGraph.nodes[key].children = newChild
			// fmt.Println("append new child", newGraph.nodes[key].label, newGraph.nodes[key].children, "old child", graph.nodes[key].children)
		}
	}
	for key, val := range graph.nodes {
		var newparents []*Node
		for i := 0; i < len(val.parents); i++ {

			newparents = append(newparents, newGraph.nodes[graph.nodes[key].parents[i].label])
			newGraph.nodes[key].parents = newparents
			// fmt.Println("append new parents", newGraph.nodes[key].parents, "old [parents]", graph.nodes[key].parents)
		}
	}
	// fmt.Println("old graoh", graph.nodes, "new graphj", newGraph.nodes)
	for key, val := range graph.edges {
		// fmt.Println("edges", key, val)
		var newEdge Edge
		var pointertoNewEdge *Edge
		pointertoNewEdge = &newEdge
		newEdge.label = val.label
		newEdge.weight = val.weight
		// fmt.Println("weight", newEdge.weight)

		var newFrom *Node = newGraph.nodes[val.from.label]
		newEdge.from = newFrom
		// fmt.Println("from,", newEdge.from)
		newEdge.to = newGraph.nodes[val.to.label]
		newEdges[key] = pointertoNewEdge
		// fmt.Println("checking each edges and make sure they have the same attributes but the memory address is not the same", "check the label of the old edge:", key, "new label:", newEdge.label, "old memory address", &val, "new memory address:", &newEdge, "old edge weight:", val.weight, "new weight:", newEdge.weight, "old from/to memory address", &val.from, &val.to, "new from /to  memory address", &newEdge.from, &newEdge.to, "old from/to label", val.from.label, val.to.label, "new from/to label:", newEdge.from.label, newEdge.to.label, "old edge memory address", &val, "new edge memory address", &newEdge)

	}
	if len(newGraph.edges) != len(graph.edges) {
		panic("Error: the len of map edges are not the same length")

	}
	// fmt.Println("finished edges, finished all")
	// // fmt.Println("graph vs new graph", graph.nodes)
	// fmt.Println("newGraph", newGraph.nodes)
	// fmt.Println("new root", newGraph.root, "old root", graph.root)
	return newGraph
}
//Generate a txt file based on uniqueKmerCounts
//Input: Unique kmer counts, kmer length
//output: txt file 
//tianyue
func GenerateUniqueKmerCountTextFile(uniqueKmerCounts map[int]int, kmerLength int) {

	f, err := os.Create(strconv.Itoa(kmerLength) + ".txt")
	if err != nil {
		fmt.Println(err)
	}
	w := bufio.NewWriter(f)
	// var list []int
	for key, val := range uniqueKmerCounts {
		for i := 0; i < val; i++ {
			w.WriteString(strconv.Itoa(key) + ",")
		}
	}
	w.Flush()

}
//Generate a text file for read coverage plotting
//Input: read plot 
//output: txt file
//tianyue
func arraytextfile(GenerateReadsPlot []int, name string) {
	f, err := os.Create(name + ".txt")
	if err != nil {
		fmt.Println(err)
	}
	w := bufio.NewWriter(f)
	for idx := range GenerateReadsPlot {
		w.WriteString(strconv.Itoa(GenerateReadsPlot[idx]) + ",")
	}
}
//input
func GetInputForHistogram(uniqueKmerCounts map[int]int) []int {
	var result []int
	for key, val := range uniqueKmerCounts {
		for i := 0; i < val; i++ {
			result = append(result, key)
		}
	}
	return result
}
