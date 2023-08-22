package fasta2

import (
	"bufio"
	"bytes"
	"fmt"
	"io"
	"os"
)

// Parser is a fasta parser it is initialized by the NewParser() function.
type Parser struct {
	buff    bytes.Buffer
	header  string
	start   bool
	scanner *bufio.Scanner
	line    int
	more    bool
}

func NewParser(r io.Reader) *Parser {
	return &Parser{
		start:   true,
		more:    true,
		scanner: bufio.NewScanner(r),
	}
}

// Lines returns the number of lines parsed.
func (p *Parser) Lines() int {
	return p.line
}

// HasNext returns true if the parser can continue parsing.
func (p *Parser) HasNext() bool {
	return p.more
}

func (p *Parser) newRecord() Record {
	sequence := p.buff.String()
	record := Record{
		Header:   p.header,
		Sequence: sequence,
	}
	// Reset sequence buffer
	p.buff.Reset()
	return record
}

// Next parsed the next record in the io.Reader and returns it, in case
// something went wrong an error and the partial result is returned.
func (p *Parser) Next() (Record, error) {
	for p.scanner.Scan() {
		line := p.scanner.Bytes()
		p.line++
		switch {
		// if there's nothing on this line skip this iteration of the loop
		case len(line) == 0:
			continue
		// if it's a comment skip this line
		case line[0] == ';':
			continue
		// start of file with no header, error
		case line[0] != '>' && p.start:
			err := fmt.Errorf("invalid input: missing sequence header for sequence starting at line %d", p.line)
			record := p.newRecord()
			return record, err
		// start of a fasta line
		case line[0] != '>':
			p.buff.Write(line)
		// Process normal new lines
		case line[0] == '>' && !p.start:
			record := p.newRecord()
			// New name
			p.header = string(line[1:])
			return record, nil
		// Process first line of file
		case line[0] == '>' && p.start:
			p.header = string(line[1:])
			p.start = false
		}
	}
	p.more = false
	// Add final sequence in file
	record := p.newRecord()
	return record, p.scanner.Err()
}

// ParseAll will parse all the records found in the reader and returns them in
// a slice.
func ParseAll(r io.Reader) ([]Record, error) {
	var (
		ret []Record
		p   = NewParser(r)
	)

	for p.HasNext() {
		rec, err := p.Next()
		if err != nil {
			return ret, err
		}
		ret = append(ret, rec)
	}

	return ret, nil
}

// ReadFile will parse all the records found in the file and returns them in
// a slice.
func ReadFile(path string) ([]Record, error) {
	var ret []Record
	f, err := os.Open(path)
	if err != nil {
		return nil, fmt.Errorf("error while reading file %q: %w", path, err)
	}
	defer f.Close()

	p := NewParser(f)
	for p.HasNext() {
		rec, err := p.Next()
		if err != nil {
			return ret, err
		}
		ret = append(ret, rec)
	}

	return ret, nil
}
