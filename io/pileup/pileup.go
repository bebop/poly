/*
Package pileup contains pileup parsers and writers.
*/
package pileup

import (
	"bufio"
	"fmt"
	"io"
	"strconv"
	"strings"
	"unicode"
)

// https://en.wikipedia.org/wiki/Pileup_format

type Pileup struct {
	Sequence       string   `json:"sequence"`
	Position       uint     `json:"position"`
	ReferenceBase  string   `json:"reference_base"`
	ReadCount      uint     `json:"read_count"`
	ReadResults    []string `json:"read_results"`
	Quality        string   `json:"quality"`
	Starts         uint     `json:"starts"`
	StartQualities []string `json:"start_qualities"`
	Ends           uint     `json:"ends"`
}

type Parser struct {
	reader bufio.Reader
	line   uint
}

func NewParser(r io.Reader) *Parser {
	const maxLineSize = 2 * 32 * 1024
	return &Parser{
		reader: *bufio.NewReaderSize(r, maxLineSize),
	}
}

func (parser *Parser) ParseNext() (Pileup, error) {
	if _, err := parser.reader.Peek(1); err != nil {
		// Early return on error. Probably will be EOF.
		return Pileup{}, err
	}
	// Parse out a single line
	lineBytes, err := parser.reader.ReadSlice('\n')
	parser.line++
	line := string(lineBytes)
	line = line[:len(line)-1] // Exclude newline delimiter.

	// Check that there are 6 values, as defined by the pileup format
	values := strings.Split(line, "\t")
	if len(values) != 6 {
		return Pileup{}, fmt.Errorf("Error on line %d: Got %d values, expected 6.", parser.line, len(values))
	}

	// Convert Position and ReadCount to integers
	positionInteger, err := strconv.Atoi(values[1])
	if err != nil {
		return Pileup{}, fmt.Errorf("Error on line %d. Got error: %w", parser.line, err)
	}
	readCountInteger, err := strconv.Atoi(values[3])
	if err != nil {
		return Pileup{}, fmt.Errorf("Error on line %d. Got error: %w", parser.line, err)
	}

	// Parse ReadResults
	readResults := make([]string, readCountInteger)
	var starts uint
	var ends uint
	var skip int
	var readCount uint
	resultsString := values[4]
	for resultIndex := range resultsString {
		if skip != 0 {
			skip = skip - 1
			continue
		}
		resultRune := resultsString[resultIndex]
		switch resultRune {
		case '^':
			starts = starts + 1
			skip = skip + 2
			readResults[readCount] = resultsString[resultIndex : resultIndex+2]
		case '$':
			ends = ends + 1
			// This applies to the last read segement
			readResults[readCount] = readResults[readCount] + "$"
		case '.', ',', '*', 'A', 'T', 'G', 'C', 'N', 'a', 't', 'g', 'c', 'n':
			readResults[readCount] = string(resultRune)
		case '-', '+':
			// formatted in `+4ATGC` format. We need to know the number of jumps
			// because you can have +10AAAAAAAAAA
			var numberOfJumps string
			for numberIndex := range resultsString[resultIndex:] {
				runeToCheck := resultsString[resultIndex+numberIndex+1]
				if unicode.IsDigit(rune(runeToCheck)) {
					numberOfJumps = numberOfJumps + string(runeToCheck)
					continue
				}
				break
			}
			regularExpressionInt, err := strconv.Atoi(numberOfJumps)
			if err != nil {
				return Pileup{}, fmt.Errorf("Error on line %d within read results. Got error: %w", parser.line, err)
			}
			readResults[readCount] = resultsString[resultIndex : resultIndex+regularExpressionInt+2]
			//fmt.Println(readResults[readCount])
			skip = skip + regularExpressionInt + len(numberOfJumps) // The 1 makes sure to include the regularExpressionInt in readResult string
		}
		fmt.Printf("ReadCount %d: %s\n", readCount, readResults[readCount])
		readCount = readCount + 1
	}

	return Pileup{Sequence: values[0], Position: uint(positionInteger), ReferenceBase: values[2], ReadCount: uint(readCountInteger), ReadResults: readResults, Quality: values[5]}, nil
}
