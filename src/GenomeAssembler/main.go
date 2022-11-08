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
