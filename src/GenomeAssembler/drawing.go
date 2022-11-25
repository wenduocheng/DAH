package main

import (
	"bufio"
	"fmt"
	"os"
	"strconv"

	"github.com/dominikbraun/graph"
	"github.com/dominikbraun/graph/draw"
)

// func main() {
// 	// fmt.Println("Hello, 世界")
// 	var sudocounts = []int{233, 125, 98, 70, 23, 4, 1}

// 	DrawKmerScatter(sudocounts)
// }

// Draw Scatter plot
// Input: UniqueKmerCount []int
// Output: Scatter plot
//tianyue

// Draw graph based on given De Bruijn Graph
// Input: De Bruijn Graph (datatype as nodes map[string]*Node //k-1 mer; edges map[string]*Edge)
// Output: De Bruijn Graph where vertices are k-1 mers and edges are k-mers
//tianyue
func (deBruijn Graph) DrawDeBruijnGraph() {
	deBruijn_plot := graph.New(graph.StringHash, graph.Directed())

	for key := range deBruijn.nodes {
		_ = deBruijn_plot.AddVertex(deBruijn.nodes[key].label)
		fmt.Println(deBruijn.nodes[key].label)
		// _ = g.AddEdge(dbj.edges[key].from.label, dbj.edges[key].to.label)
		// fmt.Println(dbj.edges[key].from.label, dbj.edges[key].to.label)

	}
	for key := range deBruijn.edges {
		_ = deBruijn_plot.AddEdge(deBruijn.edges[key].from.label, deBruijn.edges[key].to.label)
		fmt.Println(deBruijn.edges[key].from.label, deBruijn.edges[key].to.label)

	}
	file, _ := os.Create("./mygraph.gv")
	_ = draw.DOT(deBruijn_plot, file)

}

// In our graph construct, we have Node as a map where keys represent the k-1 mers where we are connecting to the next one, and generate a gfa file that can be read by bandage
// Input: de Bruijn Graph
// output: GFA file
// tianyue
func SaveGraphToGFA(graph Graph) {
	fmt.Println("SaveGraphToGFA")
	//create a reference map
	reference := make(map[string]int)
	counter := 0
	for key := range graph.nodes {
		counter += 1
		reference[key] = counter
	}

	f, err := os.Create("output.gfa")
	if err != nil {
		fmt.Println(err)
	}
	w := bufio.NewWriter(f)
	w.WriteString("H" + "\t" + "\n")

	//write the segment part
	for key := range graph.nodes {
		idx := reference[key]

		w.WriteString("S" + "\t" + strconv.Itoa(idx) + "\t" + key + "\n")

		// w.WriteString(key + "\n")
	}
	//write the linkage part
	// counter_2 :=0
	for _, val := range graph.edges {
		from := reference[val.from.label]
		to := reference[val.to.label]
		w.WriteString("L" + "\t" + strconv.Itoa(from) + "\t" + "+" + "\t" + strconv.Itoa(to) + "\t" + "+" + "\t" + "5M" + "\n")

	}
	w.Flush()

}
