/*
Package fastq contains fastq parsers and writers.

Fastq is a flat text file format developed in ~2000 to store nucleotide
sequencing data. While similar to fastq, fastq has a few differences. First,
the sequence identifier begins with @ instead of >, and includes quality
values for a sequence.

This package provides a parser and writer for working with Fastq formatted
sequencing data.
*/
package fastq

import (
	"bufio"
	"errors"
	"fmt"
	"io"
	"strings"
)

/******************************************************************************
March 22, 2023

Fastq Parser begins here

I basically stole everything from the fasta parser, and added a few bits
for parsing out additional data from fastq nanopore files. Mwhahaha, stealing
code!

Keoni

******************************************************************************/

// Read is a struct representing a single Fastq read element with an Identifier, its corresponding sequence, its quality score, and any optional pieces of data.
type Read struct {
	Identifier string            `json:"identifier"`
	Optionals  map[string]string `json:"optionals"` // Nanopore, for example, carries along data like: `read=13956 ch=53 start_time=2020-11-11T01:49:01Z`
	Sequence   string            `json:"sequence"`
	Quality    string            `json:"quality"`
}

// Header is a blank struct, needed for compatibility with bio parsers. It contains nothing.
type Header struct{}

// WriteTo is a blank function, needed for compatibility with bio parsers. It doesn't do anything.
func (header *Header) WriteTo(w io.Writer) (int64, error) {
	return 0, nil
}

// Parser is a flexible parser that provides ample
// control over reading fastq-formatted sequences.
// It is initialized with NewParser.
type Parser struct {
	// reader keeps state of current reader.
	reader bufio.Reader
	line   uint
	atEOF  bool
}

// Header returns nil,nil.
func (parser *Parser) Header() (*Header, error) {
	return &Header{}, nil
}

// NewParser returns a Parser that uses r as the source
// from which to parse fastq formatted sequences.
func NewParser(r io.Reader, maxLineSize int) *Parser {
	return &Parser{
		reader: *bufio.NewReaderSize(r, maxLineSize),
	}
}

// Next reads next fastq genome in underlying reader and returns the result
// and the amount of bytes read during the call.
// Next only returns an error if it:
//   - Attempts to read and fails to find a valid fastq sequence.
//   - Returns reader's EOF if called after reader has been exhausted.
//   - If a EOF is encountered immediately after a sequence with no newline ending.
//     In this case the Read up to that point is returned with an EOF error.
//
// It is worth noting the amount of bytes read are always right up to before
// the next fastq starts which means this function can effectively be used
// to index where fastqs start in a file or string.
//
// Next is simplified for fastq files from fasta files. Unlike fasta
// files, fastq always have 4 lines following each other - not variable with
// a line limit of 80 like fasta files have. So instead of a for loop, you
// can just parse 4 lines at once.
func (parser *Parser) Next() (*Read, error) {
	if parser.atEOF {
		return &Read{}, io.EOF
	}
	// Initialization of parser state variables.
	var (
		// Parser looks for a line starting with '@'
		// that contains the next fastq sequence identifier.
		lookingForIdentifier   = true
		seqIdentifier, quality string
		optionals              map[string]string
		sequence, line         []byte
		err                    error
	)

	// parse identifier
	line, err = parser.reader.ReadSlice('\n')
	parser.line++
	if err != nil {
		return &Read{}, err
	}

	line = line[:len(line)-1] // Exclude newline delimiter.
	if string(line)[0] == '@' {
		lookingForIdentifier = false
	}
	lineSplits := strings.Split(string(line), " ")
	seqIdentifier = lineSplits[0][1:]
	optionals = make(map[string]string)
	for _, optionalDatum := range lineSplits[1:] {
		optionalSplits := strings.Split(optionalDatum, "=")
		optionalKey := optionalSplits[0]
		optionalValue := optionalSplits[1]
		optionals[optionalKey] = optionalValue
	}

	// parse sequence
	line, err = parser.reader.ReadSlice('\n')
	parser.line++
	if err != nil {
		return &Read{}, err
	}
	if len(line) <= 1 { // newline delimiter - actually checking for empty line
		return &Read{}, fmt.Errorf("empty fastq sequence for %q,  got to line %d: %w", seqIdentifier, parser.line, err)
	}
	sequence = line[:len(line)-1] // Exclude newline delimiter.
	var newSequence []byte
	for _, char := range sequence {
		newSequence = append(newSequence, char)
		if !strings.ContainsRune("ATGCN", rune(char)) {
			return &Read{}, errors.New("Only letters ATGCN are allowed for DNA/RNA in fastq file. Got letter: " + string(char))
		}
	}

	// skip +
	_, err = parser.reader.ReadSlice('\n')
	parser.line++
	if err != nil {
		return &Read{}, err
	}

	// parse quality
	line, err = parser.reader.ReadSlice('\n')
	parser.line++
	if err != nil {
		// If the line is EOF, just continue and finish off returning the new read.
		if err != io.EOF {
			return &Read{}, nil
		}
		parser.atEOF = true
	}
	if len(line) <= 1 { // newline delimiter - actually checking for empty line
		return &Read{}, fmt.Errorf("empty quality sequence for %q,  got to line %d: %w", seqIdentifier, parser.line, err)
	}
	quality = string(line[:len(line)-1])

	// Parsing ended. Check for inconsistencies.
	if lookingForIdentifier {
		return &Read{}, fmt.Errorf("did not find fastq start '@', got to line %d: %w", parser.line, err)
	}
	fastq := Read{
		Identifier: seqIdentifier,
		Optionals:  optionals,
		Quality:    quality,
		Sequence:   string(newSequence),
	}
	// Gotten to this point err is non-nil only in EOF case.
	// We report this error to note the fastq may be incomplete/corrupt
	// like in the case of using an io.LimitReader wrapping the underlying reader.
	return &fastq, nil
}

// Reset discards all data in buffer and resets state.
func (parser *Parser) Reset(r io.Reader) {
	parser.reader.Reset(r)
	parser.line = 0
}

/******************************************************************************

Start of  Write functions

******************************************************************************/

func (read *Read) WriteTo(w io.Writer) (int64, error) {
	var err error
	var writtenBytes int64
	var newWrittenBytes int
	newWrittenBytes, err = fmt.Fprintf(w, "@%s", read.Identifier)
	writtenBytes += int64(newWrittenBytes)
	if err != nil {
		return writtenBytes, err
	}
	for key, val := range read.Optionals {
		newWrittenBytes, err = fmt.Fprintf(w, " %s=%s", key, val)
		writtenBytes += int64(newWrittenBytes)
		if err != nil {
			return writtenBytes, err
		}
	}

	// fastq doesn't limit at 80 characters, since it is
	// mainly reading big ole' sequencing files without
	// human input.
	newWrittenBytes, err = fmt.Fprintf(w, "\n%s\n+\n%s\n", read.Sequence, read.Quality)
	writtenBytes += int64(newWrittenBytes)
	return writtenBytes, err
}
