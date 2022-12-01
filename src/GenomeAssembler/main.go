package main

import "fmt"

func main() {

		viralSequence := "TTGTAGATCTGTTCTCTAAACGAACTTTAAAATCTGTGTGGCTGTCACTCGGCTGCATGCTTAGTGCACTCACGCAGTATAATTAATAACTAATTACTGTCGTTGACAGGACACGAGTAACTCGTCTATCTTCTGCAGGCTGCTTACGGTTTCGTCCGTGTTGCAGCCGATCATCAGCACATCTAGGTTTTGTCCGGGTGTGACCGAAAGGTAAGATGGAGAGCCTTGTCCCTGGTTTCAACGAGAAAACACACGTCCAACTCAGTTTGCCTGTTTTACAGGTTCGCGACGTGCTCGTACGTGGCTTTGGAGACTCCGTGGAGGAGGTCTTATCAGAGGCACGTCAACATCTTAAAGATGGCACTTGTGGCTTAGTAGAAGTTGAAAAAGGCGTTTTGCCTCAACTTGAACAGCCCTATG"
	readLength := 150
	readCounts := 3000
	reads := GenerateReads(readLength, readCounts, viralSequence)
	// fmt.Println(reads)

	kmerLength := 8

	// graph := DeBruijnGraph(kmerLength, reads)
	// fmt.Println("DeBrujinGraph was created")
	// for i := range graph.nodes {
	// 	fmt.Println("Node: ", graph.nodes[i].label, "Indegree: ", graph.nodes[i].inDegree, "Outdegree: ", graph.nodes[i].outDegree, "Children:")
	// 	for j := range graph.nodes[i].children {
	// 		fmt.Println(graph.nodes[i].children[j].label)
	// 	}
	// 	fmt.Println("Parents: ")
	// 	for j := range graph.nodes[i].parents {
	// 		fmt.Println(graph.nodes[i].parents[j].label)
	// 	}
	// }
	// for i := range graph.edges {
	// 	fmt.Println(graph.edges[i].label, graph.edges[i].from.label, graph.edges[i].to.label, graph.edges[i].weight)
	// }
	//
	// fmt.Println(graph.root)

	// mergedGraph := graph.ChainMerging()
	// fmt.Println("Chain merging was performed!")
	// for i := range mergedGraph.nodes {
	// 	fmt.Println("Node: ", mergedGraph.nodes[i].label, "Indegree: ", mergedGraph.nodes[i].inDegree, "Outdegree: ", mergedGraph.nodes[i].outDegree, "Children:")
	// 	for j := range mergedGraph.nodes[i].children {
	// 		fmt.Println(mergedGraph.nodes[i].children[j].label)
	// 	}
	// 	fmt.Println("Parents: ")
	// 	for j := range mergedGraph.nodes[i].parents {
	// 		fmt.Println(mergedGraph.nodes[i].parents[j].label)
	// 	}
	// }
	// for i := range mergedGraph.edges {
	// 	fmt.Println(mergedGraph.edges[i].label, mergedGraph.edges[i].from.label, mergedGraph.edges[i].to.label, mergedGraph.edges[i].weight)
	// }

	// fmt.Println(mergedGraph.root)

	contigs := DenovoAssembler(reads, kmerLength)

	fmt.Println("De novo assembly was finished!")

	fmt.Println(contigs)

}
