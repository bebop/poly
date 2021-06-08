package fasta

import (
	"bufio"
	"bytes"
	"compress/gzip"
	"io"
	"io/ioutil"
	"os"
	"strings"
)

/******************************************************************************
Apr 25, 2021

 Parser begins here

Many thanks to Jordan Campbell (https://github.com/0x106) for building the first
 parser for Poly and thanks to Tim Stiles (https://github.com/TimothyStiles)
for helping complete that PR. This work expands on the previous work by allowing
for concurrent  parsing and giving Poly a specific  parser subpackage,
as well as few bug fixes.

 is a very simple file format for working with DNA, RNA, or protein sequences.
It was first released in 1985 and is still widely used in bioinformatics.

https://en.wikipedia.org/wiki/_format

One interesting use of the concurrent  parser is working with the Uniprot
 dump files, which are far too large to fit into RAM. This parser is able
to easily handle those files by doing computation actively while the data dump
is getting parsed.

https://www.uniprot.org/downloads

I have removed the  Parsers from the io.go file and moved them into this
subpackage.

Hack the Planet,

Keoni

******************************************************************************/

// Fasta is a struct representing a single Fasta file element with a Name and its corresponding Sequence.
type Fasta struct {
	Name     string `json:"name"`
	Sequence string `json:"sequence"`
}

// Parse parses a given Fasta file into an array of Fasta structs. Internally, it uses ParseFastaConcurrent.
func Parse(r io.Reader) []Fasta {
	fastas := make(chan Fasta, 1000) // A buffer is used so that the functions runs as it is appending to outputFastas
	go ParseConcurrent(r, fastas)

	var outputFastas []Fasta
	for fasta := range fastas {
		outputFastas = append(outputFastas, fasta)
	}
	return outputFastas
}

// ParseConcurrent concurrently parses a given Fasta file in an io.Reader into a channel of Fasta structs.
func ParseConcurrent(r io.Reader, sequences chan<- Fasta) {
	// Initialize necessary variables
	var sequenceLines []string
	var name string
	start := true

	// Start the scanner
	scanner := bufio.NewScanner(r)
	for scanner.Scan() {
		line := scanner.Text()
		switch {
		// if there's nothing on this line skip this iteration of the loop
		case len(line) == 0:
			continue
		// if it's a comment skip this line
		case line[0:1] == ";":
			continue
		// start of a fasta line
		case line[0:1] != ">":
			sequenceLines = append(sequenceLines, line)
		// Process normal new lines
		case line[0:1] == ">" && !start:
			sequence := strings.Join(sequenceLines, "")
			newFasta := Fasta{
				Name:     name,
				Sequence: sequence}
			// Reset sequence lines
			sequenceLines = []string{}
			// New name
			name = line[1:]
			sequences <- newFasta
		// Process first line of file
		case line[0:1] == ">" && start:
			name = line[1:]
			start = false
		}
	}
	// Add final sequence in file to channel
	sequence := strings.Join(sequenceLines, "")
	newFasta := Fasta{
		Name:     name,
		Sequence: sequence}
	sequences <- newFasta
	close(sequences)
}

/******************************************************************************

Start of  Read functions

******************************************************************************/

// ReadGzConcurrent concurrently reads a gzipped Fasta file into a Fasta channel.
func ReadGzConcurrent(path string, sequences chan<- Fasta) {
	file, _ := os.Open(path)
	r, _ := gzip.NewReader(file)
	go ParseConcurrent(r, sequences)
}

// ReadConcurrent concurrently reads a flat Fasta file into a Fasta channel.
func ReadConcurrent(path string, sequences chan<- Fasta) {
	file, _ := os.Open(path)
	go ParseConcurrent(file, sequences)
}

// ReadGz reads a gzipped  file into an array of Fasta structs.
func ReadGz(path string) []Fasta {
	file, _ := os.Open(path)
	r, _ := gzip.NewReader(file)
	fastas := Parse(r)
	return fastas
}

// Read reads a  file into an array of Fasta structs
func Read(path string) []Fasta {
	file, _ := os.Open(path)
	fastas := Parse(file)
	return fastas
}

/******************************************************************************

Start of  Write functions

******************************************************************************/

// Build writes a Fasta struct to a  string.
func Build(fastas []Fasta) []byte {
	var fastaString bytes.Buffer
	for _, fasta := range fastas {
		fastaString.WriteString(">")
		fastaString.WriteString(fasta.Name)
		fastaString.WriteString("\n")
		fastaString.WriteString(fasta.Sequence)
		fastaString.WriteString("\n")
	}
	return fastaString.Bytes()
}

// Write writes a  string to a file.
func Write(fastas []Fasta, path string) {
	_ = ioutil.WriteFile(path, Build(fastas), 0644)
}
