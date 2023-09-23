/*
Package pileup contains pileup parsers and writers.

The pileup format is a text-based bioinformatics format to summarize aligned
reads against a reference sequence. In comparison to simply getting a consensus
sequence from sequencing data, pileup files can contain more context about the
mutations in a sequencing run, which is especially useful when analyzing
plasmid sequencing data from Nanopore sequencing runs.

Pileup files are basically tsv files with 6 columns: Sequence Identifier, Position,
Reference Base, Read Count, Read Results, and Quality. An example from
wikipedia (https://en.wikipedia.org/wiki/Pileup_format) is shown below:

	```
	seq1 	272 	T 	24 	,.$.....,,.,.,...,,,.,..^+. 	<<<+;<<<<<<<<<<<=<;<;7<&
	seq1 	273 	T 	23 	,.....,,.,.,...,,,.,..A 	<<<;<<<<<<<<<3<=<<<;<<+
	seq1 	274 	T 	23 	,.$....,,.,.,...,,,.,... 	7<7;<;<<<<<<<<<=<;<;<<6
	seq1 	275 	A 	23 	,$....,,.,.,...,,,.,...^l. 	<+;9*<<<<<<<<<=<<:;<<<<
	seq1 	276 	G 	22 	...T,,.,.,...,,,.,.... 	33;+<<7=7<<7<&<<1;<<6<
	seq1 	277 	T 	22 	....,,.,.,.C.,,,.,..G. 	+7<;<<<<<<<&<=<<:;<<&<
	seq1 	278 	G 	23 	....,,.,.,...,,,.,....^k. 	%38*<<;<7<<7<=<<<;<<<<<
	seq1 	279 	C 	23 	A..T,,.,.,...,,,.,..... 	75&<<<<<<<<<=<<<9<<:<<<
	```

	1. Sequence Identifier: The sequence identifier of the reference sequence
	2. Position: Position of row in the reference sequence (indexed at 1)
	3. Reference Base: Base pair in reference sequence
	4. Read Count: Number of aligned reads to this particular base pair
	5. Read Results: The resultant alignments
	6. Quality: Phred quality scores associated with each base

This package provides a parser and writer for working with pileup files.
*/
package pileup

import (
	"bufio"
	"bytes"
	"errors"
	"fmt"
	"io"
	"math"
	"os"
	"strconv"
	"strings"
	"unicode"
)

// https://en.wikipedia.org/wiki/Pileup_format

// Pileup struct is a single position in a pileup file. Pileup files "pile"
// a bunch of separate bam/sam alignments into something more readable at a per
// base pair level, so are only useful as a grouping.
type Pileup struct {
	Sequence      string   `json:"sequence"`
	Position      uint     `json:"position"`
	ReferenceBase string   `json:"reference_base"`
	ReadCount     uint     `json:"read_count"`
	ReadResults   []string `json:"read_results"`
	Quality       string   `json:"quality"`
}

// Parse parses a given Pileup file into an array of Pileup structs.
func Parse(r io.Reader) ([]Pileup, error) {
	// 32kB is a magic number often used by the Go stdlib for parsing. We multiply it by two.
	const maxLineSize = 2 * 32 * 1024
	parser := NewParser(r, maxLineSize)
	return parser.ParseAll()
}

// Parser is a pileup parser.
type Parser struct {
	reader bufio.Reader
	line   uint
}

// NewParser creates a parser from an io.Reader for pileup data.
func NewParser(r io.Reader, maxLineSize int) *Parser {
	return &Parser{
		reader: *bufio.NewReaderSize(r, maxLineSize),
	}
}

// ParseAll parses all sequences in underlying reader only returning non-EOF errors.
// It returns all valid pileup sequences up to error if encountered.
func (parser *Parser) ParseAll() ([]Pileup, error) {
	return parser.ParseN(math.MaxInt)
}

// ParseN parses up to maxRows pileup sequences from the Parser's underlying reader.
// ParseN does not return EOF if encountered.
// If an non-EOF error is encountered it returns it and all correctly parsed sequences up to then.
func (parser *Parser) ParseN(maxRows int) (pileups []Pileup, err error) {
	for counter := 0; counter < maxRows; counter++ {
		pileup, err := parser.ParseNext()
		if err != nil {
			if errors.Is(err, io.EOF) {
				err = nil // EOF not treated as parsing error.
			}
			return pileups, err
		}
		pileups = append(pileups, pileup)
	}
	return pileups, nil
}

// ParseNext parses the next pileup row in a pileup file.
// ParseNext returns an EOF if encountered.
func (parser *Parser) ParseNext() (Pileup, error) {
	if _, err := parser.reader.Peek(1); err != nil {
		// Early return on error. Probably will be EOF.
		return Pileup{}, err
	}
	// Parse out a single line
	lineBytes, err := parser.reader.ReadSlice('\n')
	if err != nil {
		return Pileup{}, err
	}
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
	var readResults []string
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
			readResults = append(readResults, resultsString[resultIndex:resultIndex+3])
		case '$':
			ends = ends + 1
			// This applies to the last read segement
			readResults[len(readResults)-1] = readResults[len(readResults)-1] + "$"
		case '.', ',', '*', 'A', 'T', 'G', 'C', 'N', 'a', 't', 'g', 'c', 'n':
			readResults = append(readResults, string(resultRune))
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
			regularExpressionInt, _ := strconv.Atoi(numberOfJumps) // Because of the above check, this will never err
			readResult := resultsString[resultIndex : resultIndex+regularExpressionInt+2]
			for _, letter := range readResult {
				switch letter {
				case '0', '1', '2', '3', '4', '5', '6', '7', '8', '9', 'A', 'T', 'G', 'C', 'N', 'a', 't', 'g', 'c', 'n', '-', '+':
					continue
				default:
					return Pileup{}, fmt.Errorf("Rune within +,- not found on line %d. Got %c: only runes allowed are: [0 1 2 3 4 5 6 7 8 9 A T G C N a t g c n - +]", parser.line, letter)
				}
			}
			readResults = append(readResults, readResult)
			skip = skip + regularExpressionInt + len(numberOfJumps) // The 1 makes sure to include the regularExpressionInt in readResult string
		default:
			return Pileup{}, fmt.Errorf("Rune not found on line %d. Got %c: only runes allowed are: [^ $ . , * A T G C N a t g c n - +]", parser.line, resultRune)
		}
		readCount = readCount + 1
	}

	return Pileup{Sequence: values[0], Position: uint(positionInteger), ReferenceBase: values[2], ReadCount: uint(readCountInteger), ReadResults: readResults, Quality: values[5]}, nil
}

// Reset discards all data in buffer and resets state.
func (parser *Parser) Reset(r io.Reader) {
	parser.reader.Reset(r)
	parser.line = 0
}

/******************************************************************************

Start of  Read functions

******************************************************************************/

// Read reads a  file into an array of Pileup structs
func Read(path string) ([]Pileup, error) {
	file, err := os.Open(path)
	if err != nil {
		return nil, err
	}
	defer file.Close()
	return Parse(file)
}

/******************************************************************************

Start of  Write functions

******************************************************************************/

// WritePileups writes a pileup array to an io.Writer
func WritePileups(pileups []Pileup, w io.Writer) error {
	for _, pileup := range pileups {
		_, err := w.Write([]byte(
			strings.Join(
				[]string{
					pileup.Sequence,
					strconv.FormatUint(uint64(pileup.Position), 10),
					pileup.ReferenceBase,
					strconv.FormatUint(uint64(pileup.ReadCount), 10),
					strings.Join(pileup.ReadResults, ""),
					pileup.Quality,
				},
				"\t") + "\n"))
		if err != nil {
			return err
		}
	}
	return nil
}

// Write writes a pileup array to a file
func Write(pileups []Pileup, path string) error {
	var b bytes.Buffer
	if err := WritePileups(pileups, &b); err != nil {
		return err
	}
	return os.WriteFile(path, b.Bytes(), 0644)
}
