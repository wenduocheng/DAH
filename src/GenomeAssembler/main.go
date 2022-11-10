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
	//Example of graphing
	//Because previous variable is conflit with package, so I change it to another one
	graph_plot := graph
	graph_plot := graph.New(graph.StringHash, graph.Directed())

	for key := range graph_plot.nodes {
		_ = g.AddVertex(graph_plot.nodes[key].label)
		fmt.Println(graph_plot.nodes[key].label)
		//cannot label edges as integer so please ignore these codes
		// _ = g.AddEdge(dbj.edges[key].from.label, dbj.edges[key].to.label)
		// fmt.Println(dbj.edges[key].from.label, dbj.edges[key].to.label)

	}
	for key := range graph_plot.edges {
		_ = g.AddEdge(graph_plot.edges[key].from.label, graph_plot.edges[key].to.label)
		fmt.Println(graph_plot.edges[key].from.label, graph_plot.edges[key].to.label)

	}
	//The file is generated as a DOT file
	file, _ := os.Create("./mygraph.gv")
	_ = draw.DOT(graph_plot, file)

}
