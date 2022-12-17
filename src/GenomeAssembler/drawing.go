package main

import (
	"bufio"
	"fmt"
	"image/color"
	"os"
	"strconv"

	"github.com/dominikbraun/graph"
	"github.com/dominikbraun/graph/draw"
	"gonum.org/v1/plot"
	"gonum.org/v1/plot/plotter"
	"gonum.org/v1/plot/vg"
)


// Draw histogram
// Wenduo
func DrawHistogram(counts []int) {
	//make data
	var values plotter.Values

	for i := 0; i < len(counts); i++ {
		values = append(values, float64(counts[i]))
		// fmt.Println(float64(counts[i]))
	}

	//boxPlot(values)
	//barPlot(values[:4])
	histPlot(values)
}

// Wenduo
func histPlot(values plotter.Values) {
	p := plot.New()

	p.Title.Text = "Kmer frequency distribution"

	hist, err := plotter.NewHist(values, 18)
	if err != nil {
		panic(err)
	}
	p.Add(hist)

	if err := p.Save(3*vg.Inch, 3*vg.Inch, "kmerFrequencyDistribution.png"); err != nil {
		panic(err)
	}
}

// Draw bar plot
// Wenduo
func DrawBarPlot(counts []int) {
	//make data
	var values plotter.Values

	for i := 0; i < len(counts); i++ {
		values = append(values, float64(counts[i]))
		// fmt.Println(float64(counts[i]))
	}

	//boxPlot(values)
	//barPlot(values[:4])
	barPlot(values)
}

// Wenduo
func barPlot(values plotter.Values) {
	p := plot.New()

	p.Title.Text = "Coverage count plot"

	bar, err := plotter.NewBarChart(values, 1)
	if err != nil {
		panic(err)
	}
	p.Add(bar)

	if err := p.Save(3*vg.Inch, 3*vg.Inch, "CoverageCountPlot.png"); err != nil {
		panic(err)
	}
}



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
func SaveGraphToGFA(graph Graph, prefix string) {
	fmt.Println("SaveGraphToGFA")
	//create a reference map
	reference := make(map[string]int)
	counter := 0
	for key := range graph.nodes {
		counter += 1
		reference[key] = counter
	}

	f, err := os.Create(prefix + "output.gfa")
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
		w.WriteString("L" + "\t" + strconv.Itoa(from) + "\t" + "-" + "\t" + strconv.Itoa(to) + "\t" + "-" + "\t" + "5M" + "\n")

	}
	w.Flush()

}
//Draw a scatter plot based on uniqueKmerCounts
//Input: uniqueKmerCounts
//output: a png image
//tianyue
func DrawKmerScatter(uniqueKmerCounts map[int]int, title, xlabel, ylabel string) {
	//due to the package is using interface, I need to check more about this

	type XY struct{ X, Y float64 }

	type XYs []XY

	//generate data in XY format
	pts := make(plotter.XYs, len(uniqueKmerCounts))
	var x []int
	var y []int

	for key, val := range uniqueKmerCounts {
		x = append(x, val)
		y = append(y, key)

	}
	for i := 0; i < len(pts); i++ {
		pts[i].Y = float64(x[i])
		pts[i].X = float64(y[i])
	}

	fmt.Println(pts, "pts")

	//generate data
	// for i := 0; i < len(UniqueCounts); i++ {
	// 	Plot_val = append(Plot_val, float64(UniqueCounts[i]))
	// 	fmt.Println(Plot_val, UniqueCounts[i])
	// }

	fmt.Println("finished data copying")
	plot := plot.New()

	plot.Title.Text = title
	plot.Y.Label.Text = ylabel
	plot.X.Label.Text = xlabel
	sca, err2 := plotter.NewScatter(pts)
	if err2 != nil {
		panic("Filure generating histogram")
	}
	sca.GlyphStyle.Color = color.RGBA{R: 255, B: 128, A: 255}
	sca.GlyphStyle.Radius = vg.Points(1)
	// line, err := plotter.NewLine(LineData)
	// if err != nil {
	// 	log.Panic(err)
	// }
	// line.LineStyle.Width = vg.Points(1)
	// line.LineStyle.Dashes = []vg.Length{vg.Points(5), vg.Points(5)}
	// line.LineStyle.Color = color.RGBA{B: 255, A: 255}
	plot.Add(sca)
	plot.Save(200, 200, "scatter_trial.png")

}
