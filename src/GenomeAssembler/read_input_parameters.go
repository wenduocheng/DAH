package main

//add small sample data for instruction about what this function is taking and what it is outputing

import (
	"fmt"
	"io/ioutil"
	"os"
	"strconv"
	"strings"
)

// func main() {
// 	//those are test cases, please change/ignore
// 	inputdir := "C:/Users/tiant/go/src/GenomeAssembler/input/"
// 	var inputfile []os.FileInfo
// 	inputfile = ReadFilesFromDirectory(inputdir)
// 	fmt.Println(inputfile)
// 	fmt.Print(ReadTwoIntandTwoPairString(inputdir, inputfile[2]))

// }

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

//read and input the sequencing infomation from a text file in the given directory
//input: genome sequence file directory; or if in the same folder, just input the name of the file
//output: a string of a sequence
func ReadSequence(path string) string {
	fileContents, err := ioutil.ReadFile(path)
	if err != nil {
		panic(err)
	}
	inputLines := strings.Split(strings.TrimSpace(strings.Replace(string(fileContents), "\r\n", "\n", -1)), "\n")
	inputLines = inputLines[1:]
	var text string
	text = strings.Join(inputLines, text)
	text = GetRidOfN(text)
	var result string
	for i := 0; i < len(text); i++ {
		if string(text[i]) == "A" || string(text[i]) == "a" || string(text[i]) == "C" || string(text[i]) == "c" || string(text[i]) == "T" || string(text[i]) == "t" || string(text[i]) == "G" || string(text[i]) == "g" || string(text[i]) == " " {
			result += string(text[i])
		}
	}
	// fmt.Println(result, "GE")
	return result
}

//ReadReads takes a fasta file and keep the sequencing reads into a slice of strings
//input: fasta file
//output: []strings of illumina sequencing short reads
//tianyue

func ReadReads(path string) []string {
	//get and read the file
	fileContents, err := ioutil.ReadFile(path)
	if err != nil {
		panic(err)
	}
	//trim lines
	inputLines := strings.Split(strings.TrimSpace(strings.Replace(string(fileContents), "\r\n", "\n", -1)), "\n")
	var text []string
	text = inputLines
	var reads []string
	//only keep the information that we need, but do consider if we still want to keep N in our test sets
	for i := 0; i < len(text); i++ {
		if text[i][0:1] != ">" {
			reads = append(reads, text[i])

		}

	}

	return reads
}

//GetRidOf N will go through all the bases in the sequence and remove allthe part with N
//input : stringh
//output: string with no N inside
func GetRidOfN(sequence string) string {
	var result string
	for i := 0; i < len(sequence); i++ {
		if string(sequence[i]) != "N" {
			result += string(sequence[i])
			// fmt.Println(string(sequence[i]))
		}
	}

	return result
}
