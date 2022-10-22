package main

//add small sample data for instruction about what this function is taking and what it is outputing

import (
	"fmt"
	"io/ioutil"
	"os"
	"strconv"
	"strings"
)

func main() {
	//those are test cases, please change/ignore
	inputdir := "C:/Users/tiant/go/src/GenomeAssembler/input/"
	var inputfile []os.FileInfo
	inputfile = ReadFilesFromDirectory(inputdir)
	fmt.Println(inputfile)
	fmt.Print(ReadTwoIntandTwoPairString(inputdir, inputfile[2]))

}

// by given directory, load ALL the files under that specific directory
// input: directory of the file that you want to take into a function
// for windows: "C:/Users/tiant/go/src/go_parctice/" and make sure to keep the last /
func ReadFilesFromDirectory(directory string) []os.FileInfo {
	dirContents, err := ioutil.ReadDir(directory)
	if err != nil {
		panic("Error reading directory: " + directory)
	}

	return dirContents
}

// For a given text file with integer and string, output parameters for integer and string
// eg: 50
// "GATCTTAGTTTGACGCTTTGAAGAAAAGGCAG"
func ReadIntegerandStringFromFile(directory string, input_file os.FileInfo) (int, string) {
	//now, consult the associated output file.
	fileName := input_file.Name() //grab file name

	//now, read out the file
	fileContents, err := ioutil.ReadFile(directory + fileName)
	if err != nil {
		panic(err)
	}

	//trim out extra space and store as a slice of strings, each containing one line.
	inputLines := strings.Split(strings.TrimSpace(strings.Replace(string(fileContents), "\r\n", "\n", -1)), "\n")

	//parse the float
	integer, err := strconv.Atoi(inputLines[0])

	if err != nil {
		panic(err)
	}
	var text string
	text = strings.Join(inputLines[1:], text)

	return integer, text
}

// Read a text file with slice of strings patterns (ie k-mer patterns) and output []string
// eg input:GTGAAGACAAAGTTAGGGTTCGCGA
// TGAAGACAAAGTTAGGGTTCGCGAT
// GAAGACAAAGTTAGGGTTCGCGATA
// AAGACAAAGTTAGGGTTCGCGATAC
// AGACAAAGTTAGGGTTCGCGATACA
func ReadStringSliceFromFile(directory string, input_file os.FileInfo) []string {
	fileName := input_file.Name() //grab file name

	//now, read out the file
	fileContents, err := ioutil.ReadFile(directory + fileName)
	if err != nil {
		panic(err)
	}

	//trim out extra space and store as a slice of strings, each containing one line.
	inputLines := strings.Split(strings.TrimSpace(strings.Replace(string(fileContents), "\r\n", "\n", -1)), "\n")

	return inputLines
}

// for adjancency list, I am not sure what kind of function need to write so I will leave this for now
// For input with 2 different integers and a []STRING, specifically for Integers k and d followed by a collection of paired k-mers PairedReads.
// eg input: 30 100
// GACAAGCGTATTTCCTTTAGAGGACGTTAT|CTTGGATCGCGTTGGGATCTTGGATCGCGA
// AATTGGGATCTTGGATCGCGTTTCCGTCTG|TCTTGGATCGCGTTGGTAACCCTGGCTCAA
// TGAAAGGCAACACCAAACGGATACCCGTAT|TCGTGTGATCGCCTTTTTAGGGGCGGTTAT
func ReadTwoIntandTwoPairString(directory string, input_file os.FileInfo) (int, int, []string) {
	//now, consult the associated output file.
	fileName := input_file.Name() //grab file name

	//now, read out the file
	fileContents, err := ioutil.ReadFile(directory + fileName)
	if err != nil {
		panic(err)
	}

	//trim out extra space and store as a slice of strings, each containing one line.
	inputLines := strings.Split(strings.TrimSpace(strings.Replace(string(fileContents), "\r\n", "\n", -1)), "\n")
	fmt.Println(inputLines[0])
	intinput := strings.Split(inputLines[0], " ")
	fmt.Println(intinput)

	//parse the float
	integer1, err := strconv.Atoi(intinput[0])

	if err != nil {
		panic("int1 error")
	}
	integer2, err := strconv.Atoi(intinput[1])

	if err != nil {
		panic("int2 error")
	}
	var text []string
	text = inputLines[1:]

	return integer1, integer2, text

}

//map as input (adj matrix) ->
