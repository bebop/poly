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
	"compress/gzip"
	"errors"
	"fmt"
	"io"
	"math"
	"os"
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

// Fastq is a struct representing a single Fastq file element with an Identifier, its corresponding sequence, its quality score, and any optional pieces of data.
type Fastq struct {
	Identifier string            `json:"identifier"`
	Optionals  map[string]string `json:"optionals"` // Nanopore, for example, carries along data like: read=13956 ch=53 start_time=2020-11-11T01:49:01Z
	Sequence   string            `json:"sequence"`
	Quality    string            `json:"quality"`
}

// Parse parses a given Fastq file into an array of Fastq structs. Internally, it uses ParseFastqConcurrent.
func Parse(r io.Reader) ([]Fastq, error) {
	// 32kB is a magic number often used by the Go stdlib for parsing. We multiply it by two.
	const maxLineSize = 2 * 32 * 1024
	parser := NewParser(r, maxLineSize)
	return parser.ParseAll()
}

// Parser is a flexible parser that provides ample
// control over reading fastq-formatted sequences.
// It is initialized with NewParser.
type Parser struct {
	// reader keeps state of current reader.
	reader bufio.Reader
	line   uint
}

// NewParser returns a Parser that uses r as the source
// from which to parse fastq formatted sequences.
func NewParser(r io.Reader, maxLineSize int) *Parser {
	return &Parser{
		reader: *bufio.NewReaderSize(r, maxLineSize),
	}
}

// ParseAll parses all sequences in underlying reader only returning non-EOF errors.
// It returns all valid fastq sequences up to error if encountered.
func (parser *Parser) ParseAll() ([]Fastq, error) {
	return parser.ParseN(math.MaxInt)
}

// ParseN parses up to maxSequences fastq sequences from the Parser's underlying reader.
// ParseN does not return EOF if encountered.
// If an non-EOF error is encountered it returns it and all correctly parsed sequences up to then.
func (parser *Parser) ParseN(maxSequences int) (fastqs []Fastq, err error) {
	for counter := 0; counter < maxSequences; counter++ {
		fastq, _, err := parser.ParseNext()
		if err != nil {
			if errors.Is(err, io.EOF) {
				err = nil // EOF not treated as parsing error.
			}
			return fastqs, err
		}
		fastqs = append(fastqs, fastq)
	}
	return fastqs, nil
}

// ParseNext reads next fastq genome in underlying reader and returns the result
// and the amount of bytes read during the call.
// ParseNext only returns an error if it:
//   - Attempts to read and fails to find a valid fastq sequence.
//   - Returns reader's EOF if called after reader has been exhausted.
//   - If a EOF is encountered immediately after a sequence with no newline ending.
//     In this case the Fastq up to that point is returned with an EOF error.
//
// It is worth noting the amount of bytes read are always right up to before
// the next fastq starts which means this function can effectively be used
// to index where fastqs start in a file or string.
//
// ParseNext is simplified for fastq files from fasta files. Unlike fasta
// files, fastq always have 4 lines following each other - not variable with
// a line limit of 80 like fasta files have. So instead of a for loop, you
// can just parse 4 lines at once.
func (parser *Parser) ParseNext() (Fastq, int64, error) {
	if _, err := parser.reader.Peek(1); err != nil {
		// Early return on error. Probably will be EOF.
		return Fastq{}, 0, err
	}

	// More general case of error handling.
	handleErr := func(err error) error {
		isEOF := errors.Is(err, io.EOF)
		if errors.Is(err, bufio.ErrBufferFull) {
			// Buffer size too small to read fastq line.
			return fmt.Errorf("line %d too large for buffer, use larger maxLineSize: %w", parser.line+1, err)
		} else if isEOF {
			return fmt.Errorf("line %d failed: unexepcted EOF encountered", parser.line+1)
		}
		return err
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
		totalRead              int64
	)

	// parse identifier
	line, err = parser.reader.ReadSlice('\n')
	totalRead += int64(len(line))
	parser.line++
	if handleErr(err) != nil {
		return Fastq{}, totalRead, handleErr(err)
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
	totalRead += int64(len(line))
	parser.line++
	if handleErr(err) != nil {
		return Fastq{}, totalRead, handleErr(err)
	}
	if len(line) <= 1 { // newline delimiter - actually checking for empty line
		return Fastq{}, totalRead, fmt.Errorf("empty fastq sequence for %q,  got to line %d: %w", seqIdentifier, parser.line, err)
	}
	sequence = line[:len(line)-1] // Exclude newline delimiter.
	var newSequence []byte
	for _, char := range sequence {
		newSequence = append(newSequence, char)
		if !strings.ContainsRune("ATGCN", rune(char)) {
			return Fastq{}, totalRead, errors.New("Only letters ATGCN are allowed for DNA/RNA in fastq file. Got letter: " + string(char))
		}
	}

	// skip +
	_, err = parser.reader.ReadSlice('\n')
	totalRead += int64(len(line))
	parser.line++
	if handleErr(err) != nil {
		return Fastq{}, totalRead, handleErr(err)
	}

	// parse quality
	line, err = parser.reader.ReadSlice('\n')
	totalRead += int64(len(line))
	parser.line++
	if handleErr(err) != nil {
		return Fastq{}, totalRead, handleErr(err)
	}
	if len(line) <= 1 { // newline delimiter - actually checking for empty line
		return Fastq{}, totalRead, fmt.Errorf("empty quality sequence for %q,  got to line %d: %w", seqIdentifier, parser.line, err)
	}
	quality = string(line[:len(line)-1])

	// Parsing ended. Check for inconsistencies.
	if lookingForIdentifier {
		return Fastq{}, totalRead, fmt.Errorf("did not find fastq start '@', got to line %d: %w", parser.line, err)
	}
	fastq := Fastq{
		Identifier: seqIdentifier,
		Optionals:  optionals,
		Quality:    quality,
		Sequence:   string(newSequence),
	}
	// Gotten to this point err is non-nil only in EOF case.
	// We report this error to note the fastq may be incomplete/corrupt
	// like in the case of using an io.LimitReader wrapping the underlying reader.
	return fastq, totalRead, err
}

// Reset discards all data in buffer and resets state.
func (parser *Parser) Reset(r io.Reader) {
	parser.reader.Reset(r)
	parser.line = 0
}

/******************************************************************************

Start of  Read functions

******************************************************************************/

// ReadGz reads a gzipped file into an array of Fastq structs.
func ReadGz(path string) ([]Fastq, error) {
	file, err := os.Open(path)
	if err != nil {
		return nil, err
	}
	defer file.Close()
	reader, err := gzip.NewReader(file)
	if err != nil {
		return nil, err
	}
	defer reader.Close()
	return Parse(reader)
}

// Read reads a  file into an array of Fastq structs
func Read(path string) ([]Fastq, error) {
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

func (fastq *Fastq) Write(w io.Writer) error {
	var err error
	_, err = w.Write([]byte("@"))
	if err != nil {
		return err
	}
	_, err = w.Write([]byte(fastq.Identifier))
	if err != nil {
		return err
	}
	for key, val := range fastq.Optionals {
		_, err = w.Write([]byte(" "))
		if err != nil {
			return err
		}
		_, err = w.Write([]byte(key))
		if err != nil {
			return err
		}
		_, err = w.Write([]byte("="))
		if err != nil {
			return err
		}
		_, err = w.Write([]byte(val))
		if err != nil {
			return err
		}
	}
	_, err = w.Write([]byte("\n"))
	if err != nil {
		return err
	}

	// fastq doesn't limit at 80 characters, since it is
	// mainly reading big ole' sequencing files without
	// human input.
	_, err = w.Write([]byte(fastq.Sequence))
	if err != nil {
		return err
	}
	_, err = w.Write([]byte("\n+\n"))
	if err != nil {
		return err
	}
	_, err = w.Write([]byte(fastq.Quality))
	if err != nil {
		return err
	}
	_, err = w.Write([]byte("\n"))
	if err != nil {
		return err
	}
	return nil
}

// Write writes a fastq array to a file.
func Write(fastqs []Fastq, path string) error {
	file, err := os.OpenFile(path, os.O_WRONLY|os.O_CREATE, 0644)
	if err != nil {
		return err
	}
	defer file.Close()

	for _, fastq := range fastqs {
		err := fastq.Write(file)
		if err != nil {
			return err
		}
	}
	return nil
}
