package main

import (
	"fmt"
	"os"

	"github.com/dominikbraun/graph"
	"github.com/dominikbraun/graph/draw"
	"gonum.org/v1/plot"
	"gonum.org/v1/plot/plotter"
)

func main() {
	fmt.Println("Hello, 世界")
	var sudocounts = map[int]int{
		1: 3,
		2: 5,
		3: 8,
		4: 21,
		5: 23,
		6: 4,
		7: 8,
	}
	DrawKmerHistogram(sudocounts)
}

// Draw Histogrm; To draw histogram, make sure that you alrea
// Input: KmerCount
// Output: Histogram

func DrawKmerHistogram(counts map[int]int) {
	//due to the package is using interface, I need to check more about this

	// XYs implements the XYer interface.
	// 	type XYer interface {
	// 		// Len returns the number of x, y pairs.
	// 		Len() int

	// 		// XY returns an x, y pair.
	// 		XY(int) (x, y float64)
	// 	}
	// 	type XY struct{ X, Y float64 }

	// 	type XYs []XY

	// // XY is an x and y value.

	// 	Plot_val := make(XYs,len(counts)+1)
	// 	Plot_val_1 :=make(plotter.XYer,len(counts)+1)

	// 	for i:=1; i<len(counts)+1; i++ {
	// 		Plot_val[i].X,Plot_val[i].Y = float64(i), float64(counts[i])
	// 		// fmt.Println(key, val)
	// 	}
	// 	var Plot_val_1 XYer{
	// 		Plot_val,

	// 	}

	// for key, val := range counts {
	// 	Plot_val =

	// 		fmt.Println(Plot_val.X, Plot_val.Y)
	// }

	Plot_val := make(plotter.Values, len(counts))
	for key, val := range counts {
		Plot_val[key] = float64(val)
	}

	fmt.Println("finished data copying")
	plot := plot.New()

	plot.Title.Text = "K-mer distribution"
	plot.Y.Label.Text = "number of times Kmer seen"
	plot.X.Label.Text = "Number of Kmers"
	hist, err2 := plotter.NewHist(Plot_val, 16)
	if err2 != nil {
		panic("Filure generating histogram")
	}
	plot.Add(hist)
	plot.Save(200, 200, "histogram_trial.png")

}

// Draw graph based on given De Bruijn Graph
// Input: De Bruijn Graph (datatype as nodes map[string]*Node //k-1 mer; edges map[string]*Edge)
// Output: De Bruijn Graph where vertices are k-1 mers and edges are k-mers
func (deBruijn Graph) DrawDeBruijnGraph() {
	deBruijn_plot := graph.New(graph.StringHash, graph.Directed())

	for key := range deBruijn_plot.nodes {
		_ = deBruijn_plot.AddVertex(deBruijn_plot.nodes[key].label)
		fmt.Println(deBruijn_plot.nodes[key].label)
		// _ = g.AddEdge(dbj.edges[key].from.label, dbj.edges[key].to.label)
		// fmt.Println(dbj.edges[key].from.label, dbj.edges[key].to.label)

	}
	for key := range deBruijn_plot.edges {
		_ = deBruijn.AddEdge(deBruijn_plot.edges[key].from.label, deBruijn_plot.edges[key].to.label)
		fmt.Println(deBruijn_plot.edges[key].from.label, deBruijn_plot.edges[key].to.label)

	}
	file, _ := os.Create("./mygraph.gv")
	_ = draw.DOT(deBruijn_plot, file)

}
