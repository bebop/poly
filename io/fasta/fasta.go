/*
Package fasta contains fasta parsers and writers.

Fasta is a flat text file format developed in 1985 to store nucleotide and
amino acid sequences. It is extremely simple and well supported across many
languages. However, this simplicity means that annotation of genetic objects
is not supported.

This package provides a parser and writer for working with Fasta formatted
genetic sequences.
*/
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

Fasta Parser begins here

Many thanks to Jordan Campbell (https://github.com/0x106) for building the first
parser for Poly and thanks to Tim Stiles (https://github.com/TimothyStiles)
for helping complete that PR. This work expands on the previous work by allowing
for concurrent  parsing and giving Poly a specific  parser subpackage,
as well as few bug fixes.

Fasta is a very simple file format for working with DNA, RNA, or protein sequences.
It was first released in 1985 and is still widely used in bioinformatics.

https://en.wikipedia.org/wiki/_format

One interesting use of the concurrent  parser is working with the Uniprot
fasta dump files, which are far too large to fit into RAM. This parser is able
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
func Parse(r io.Reader) ([]Fasta, error) {
	fastas := make(chan Fasta, 1000) // A buffer is used so that the functions runs as it is appending to outputFastas
	go ParseConcurrent(r, fastas)

	var outputFastas []Fasta
	for fasta := range fastas {
		outputFastas = append(outputFastas, fasta)
	}
	return outputFastas, nil
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
	file, _ := os.Open(path) // these errors need to be handled/logged

	reader, _ := gzip.NewReader(file)
	defer reader.Close()
	go ParseConcurrent(reader, sequences)
}

// ReadConcurrent concurrently reads a flat Fasta file into a Fasta channel.
func ReadConcurrent(path string, sequences chan<- Fasta) {
	file, _ := os.Open(path) // these errors need to be handled/logged
	go ParseConcurrent(file, sequences)
}

// ReadGz reads a gzipped  file into an array of Fasta structs.
func ReadGz(path string) ([]Fasta, error) {
	file, err := os.Open(path)
	if err != nil {
		return nil, err
	}
	reader, err := gzip.NewReader(file)
	if err != nil {
		return nil, err
	}
	fastas, err := Parse(reader)
	if err != nil {
		return nil, err
	}
	return fastas, nil
}

// Read reads a  file into an array of Fasta structs
func Read(path string) ([]Fasta, error) {
	file, err := os.Open(path)
	if err != nil {
		return nil, err
	}

	fastas, err := Parse(file)
	if err != nil {
		return nil, err
	}
	return fastas, nil
}

/******************************************************************************

Start of  Write functions

******************************************************************************/

// Build writes a Fasta struct to a  string.
func Build(fastas []Fasta) ([]byte, error) {
	var fastaString bytes.Buffer
	for _, fasta := range fastas {
		fastaString.WriteString(">")
		fastaString.WriteString(fasta.Name)
		fastaString.WriteString("\n")
		fastaString.WriteString(fasta.Sequence)
		fastaString.WriteString("\n")
	}
	return fastaString.Bytes(), nil
}

// Write writes a fasta array to a file.
func Write(fastas []Fasta, path string) error {
	fastaBytes, err := Build(fastas)
	if err != nil {
		return err
	}
	err = ioutil.WriteFile(path, fastaBytes, 0644)
	return err
}
