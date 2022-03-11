/*
Package fastq contains fastq parsers and writers.

Fastq is a flat text file format to store nucleotide and amino acid sequences
along with their sequence read quality. It is extremely simple and well
supported across many languages. However, this simplicity means that annotation
of genetic objects is not supported.

This package provides a parser and writer for working with Fastq formatted
genetic sequences.
*/
package fastq

import (
	"bufio"
	"io"
	"strings"
)

// Fastq is a struct representing a single Fastq file element with a Name, its corresponding Sequence, and sequence Quality scores.
type Fastq struct {
	Name     string `json:"name"`
	Sequence string `json:"sequence"`
	Quality  string `json:"quality"`
}

// ParseConcurrent concurrently parses a given Fasta file in an io.Reader into a channel of Fasta structs.
func ParseConcurrent(r io.Reader, sequences chan<- Fastq) {
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
		// start of a fastq title line
		case line[0:1] != "@":
			sequenceLines = append(sequenceLines, line)
		// Process normal new lines
		case line[0:1] == "@" && !start:
			sequence := strings.Join(sequenceLines, "")
			newFasta := Fastq{
				Name:     name,
				Sequence: sequence}
			// Reset sequence lines
			sequenceLines = []string{}
			// New name
			name = line[1:]
			sequences <- newFasta
		// Process first line of file
		case line[0:1] == "@" && start:
			name = line[1:]
			start = false
		}
	}
	// Add final sequence in file to channel
	sequence := strings.Join(sequenceLines, "")
	newFasta := Fastq{
		Name:     name,
		Sequence: sequence}
	sequences <- newFasta
	close(sequences)
}
