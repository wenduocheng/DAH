package main

import "fmt"

func main() {

	pseudoSequence := GenerateSequence(20)
	fmt.Println(pseudoSequence)
	readLength := 10
	numberOfCopies := 3
	reads := GenerateReadsNaive(readLength, numberOfCopies, pseudoSequence)
	fmt.Println(reads)

	kmerLength := 5

	graph := DeBruijnGraph(kmerLength, reads)
	fmt.Println("DeBrujinGraph was created")
	for i := range graph.nodes {
		fmt.Println(graph.nodes[i].label, ":")
		for j := range graph.nodes[i].children {
			fmt.Println(graph.nodes[i].children[j].label)
		}
	}
	for i := range graph.edges {
		fmt.Println(graph.edges[i].label, graph.edges[i].from.label, graph.edges[i].to.label, graph.edges[i].weight)
	}

	fmt.Println(graph.root)

	contigs := DenovoAssembler(reads, kmerLength)

	fmt.Println("De novo assembly was finished!")

	fmt.Println(contigs)

}
