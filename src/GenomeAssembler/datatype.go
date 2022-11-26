package main

type Graph struct {
	root  *Node
	nodes map[string]*Node //k-1 mer
	edges map[string]*Edge //kmer
}

type Node struct {
	label     string
	inDegree  int
	outDegree int
	children  []*Node
	parents   []*Node
}

type Edge struct {
	label  string
	from   *Node
	to     *Node
	weight int
}
