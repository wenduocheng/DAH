package main

import (
	"fmt"
)

//3.3
func PathToGenome(Kmers []string) string {
	var genome string
	k := len(Kmers[0])
	for i := 0; i < len(Kmers); i++ {
		if i == 0 {
			genome += Kmers[i]
		} else {
			genome += Kmers[i][k-1 : k]
		}
	}
	return genome
}

//3.4
//construct Debrujin Graph from text
func Debrujin(k int, text string) map[string][]string {
	DebjMap := make(map[string][]string)

	for num := 0; num < len(text)-k; num++ {
		key := text[num : num+k-1]
		value := text[num+1 : num+k]

		_, isPresent := DebjMap[key]
		if !isPresent {
			DebjMap[key] = []string{}
		}
		DebjMap[key] = append(DebjMap[key], value)
	}

	return DebjMap
}

//3.5
//construct Debrujin Graph from set of kmers
//Input: a set of kmers
//Output: Debrujin Graph (map[AGG:[GGG] CAG:[AGG AGG] GAG:[AGG] GGA:[GAG] GGG:[GGG GGA]])
func KmerDeBruijin(kmers []string) map[string][]string {
	KmerDebjMap := make(map[string][]string)
	k := len(kmers[0])

	for num := 0; num < len(kmers); num++ {

		kmer := kmers[num]
		key := kmer[0 : k-1]
		value := kmer[1:k]
		_, isPresent := KmerDebjMap[key]
		if !isPresent {
			KmerDebjMap[key] = []string{}
		}
		KmerDebjMap[key] = append(KmerDebjMap[key], value)
	}

	return KmerDebjMap
}

//construct OOD of kmer for debrujin graph
//node object
//name is the prefix of the kmer
//child list is the suffix of the k-1 mer
//in is the incoming edges of the node
//out is the outcoming edges of the node
type node struct {
	name      string
	childlist []string
	in        int
	out       int
}

//graph object forms a map the
//key of the map is the prefix of the kmers
//the value of the map[key] is the pointer to the node
type graph struct {
	nodelist map[string]*node
}

//take a set of kmers and return a graph
//input: as set of kmers
//output: graph obeject
func TKmerDeBruijin(kmers []string) graph {
	var currentnode *node
	k := len(kmers[0])
	var degraph graph
	degraph.nodelist = make(map[string]*node)
	for num := 0; num < len(kmers); num++ {

		kmer := kmers[num]
		key := kmer[0 : k-1]
		value := kmer[1:k]

		_, isPresentkey := degraph.nodelist[key]
		if !isPresentkey {
			var newNodenum1 node
			newNodenum1.name = key
			degraph.nodelist[key] = &newNodenum1
		}
		currentnode = degraph.nodelist[key]
		currentnode.childlist = append(currentnode.childlist, value)
		currentnode.out += 1

		//add income edge to value node
		_, isPresentvalue := degraph.nodelist[value]
		if !isPresentvalue {
			var newNodenum2 node
			newNodenum2.name = value
			degraph.nodelist[value] = &newNodenum2
		}
		currentnode = degraph.nodelist[value]
		currentnode.in += 1
	}

	return degraph
}
func TEulerianCycle(Graph graph) []string {
	//find start point
	var start *node
	var numEdge int
	for _, n := range Graph.nodelist {
		if n.in < n.out {
			start = n
		}
		numEdge += n.in
	}
	var genome []string
	genome = append(genome, start.name)
	currentnode := start
	for i := 0; i < numEdge; i++ {
		prefix := currentnode.childlist[len(currentnode.childlist)-1]
		currentnode.childlist = currentnode.childlist[:len(currentnode.childlist)-1]
		currentnode = Graph.nodelist[prefix]
		genome = append(genome, currentnode.name)
	}
	return genome
}
func main() {
	//a := 4
	b := []string{
		"CTTA",
		"ACCA",
		"TACC",
		"GGCT",
		"GCTT",
		"TTAC"}
	graph := TKmerDeBruijin(b)
	path := TEulerianCycle(graph)
	// for k, _ := range result.nodelist {
	// 	fmt.Println(result.nodelist[k].name)
	// 	fmt.Println(result.nodelist[k].childlist)
	// 	fmt.Println(result.nodelist[k].in)
	// 	fmt.Println(result.nodelist[k].out)
	// }
	genome := PathToGenome(path)
	fmt.Println(genome)

}
