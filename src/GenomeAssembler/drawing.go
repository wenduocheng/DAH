package main

import (
	"fmt"
	"image/color"
	"log"
	"os"

	"github.com/dominikbraun/graph"
	"github.com/dominikbraun/graph/draw"
	"gonum.org/v1/plot"
	"gonum.org/v1/plot/plotter"
	"gonum.org/v1/plot/vg"
)

func main() {
	// fmt.Println("Hello, 世界")
	var sudocounts = []int{233, 125, 98, 70, 23, 4, 1}

	DrawKmerScatter(sudocounts)
}

// Draw Scatter plot
// Input: UniqueKmerCount []int
// Output: Scatter plot
//tianyue

func DrawKmerScatter(UniqueCounts []int) {
	//due to the package is using interface, I need to check more about this

	type XY struct{ X, Y float64 }

	type XYs []XY

	//generate data in XY format
	pts := make(plotter.XYs, len(UniqueCounts))
	for i := 0; i < len(UniqueCounts); i++ {
		pts[i].X = float64(i)
		pts[i].Y = float64(UniqueCounts[i])

	}
	//generate data
	// for i := 0; i < len(UniqueCounts); i++ {
	// 	Plot_val = append(Plot_val, float64(UniqueCounts[i]))
	// 	fmt.Println(Plot_val, UniqueCounts[i])
	// }
	LineData := pts
	fmt.Println("finished data copying")
	plot := plot.New()

	plot.Title.Text = "K-mer distribution"
	plot.Y.Label.Text = "number of times Kmer seen"
	plot.X.Label.Text = "Number of Kmers"
	sca, err2 := plotter.NewScatter(pts)
	if err2 != nil {
		panic("Filure generating histogram")
	}
	sca.GlyphStyle.Color = color.RGBA{R: 255, B: 128, A: 255}
	sca.GlyphStyle.Radius = vg.Points(1)
	line, err := plotter.NewLine(LineData)
	if err != nil {
		log.Panic(err)
	}
	line.LineStyle.Width = vg.Points(1)
	line.LineStyle.Dashes = []vg.Length{vg.Points(5), vg.Points(5)}
	line.LineStyle.Color = color.RGBA{B: 255, A: 255}
	plot.Add(sca, line)
	plot.Save(200, 200, "scatter_trial.png")

}

// Draw graph based on given De Bruijn Graph
// Input: De Bruijn Graph (datatype as nodes map[string]*Node //k-1 mer; edges map[string]*Edge)
// Output: De Bruijn Graph where vertices are k-1 mers and edges are k-mers
//tianyue
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
