/*
Package fasta contains fasta parsers and writers.

Fasta is a flat text file format developed in 1985 to store nucleotide and
amino acid sequences. It is extremely simple and well-supported across many
languages. However, this simplicity means that annotation of genetic objects
is not supported.

This package provides a parser and writer for working with Fasta formatted
genetic sequences.
*/
package fasta

import (
	"bufio"
	"bytes"
	"fmt"
	"io"
)

/******************************************************************************
Apr 25, 2021

Fasta Parser begins here

Many thanks to Jordan Campbell (https://github.com/0x106) for building the first
parser for Poly and thanks to Tim Stiles (https://github.com/TimothyStiles)
for helping complete that PR. This work expands on the previous work by allowing
for concurrent  parsing and giving Poly a specific parser subpackage,
as well as few bug fixes.

Fasta is a very simple file format for working with DNA, RNA, or protein sequences.
It was first released in 1985 and is still widely used in bioinformatics.

https://en.wikipedia.org/wiki/FASTA_format

One interesting use of the concurrent parser is working with the Uniprot
fasta dump files, which are far too large to fit into RAM. This parser is able
to easily handle those files by doing computation actively while the data dump
is getting parsed.

https://www.uniprot.org/downloads

I have removed the  Parsers from the io.go file and moved them into this
subpackage.

Hack the Planet,

Keoni

******************************************************************************/

// Record is a struct representing a single Record file element with a Identifier and its corresponding Sequence.
type Record struct {
	Identifier string `json:"identifier"`
	Sequence   string `json:"sequence"`
}

// Header is a blank struct, needed for compatibility with bio parsers. It contains nothing.
type Header struct{}

// WriteTo is a blank function, needed for compatibility with bio parsers. It doesn't do anything.
func (header *Header) WriteTo(w io.Writer) (int64, error) {
	return 0, nil
}

// Parser is a flexible parser that provides ample
// control over reading fasta-formatted sequences.
// It is initialized with NewParser.
type Parser struct {
	// scanner keeps state of current reader.
	scanner    bufio.Scanner
	buff       bytes.Buffer
	identifier string
	start      bool
	line       uint
	more       bool
}

// Header returns nil,nil.
func (p *Parser) Header() (*Header, error) {
	return &Header{}, nil
}

// NewParser returns a Parser that uses r as the source
// from which to parse fasta formatted sequences.
func NewParser(r io.Reader, maxLineSize int) *Parser {
	scanner := bufio.NewScanner(r)
	buf := make([]byte, maxLineSize)
	scanner.Buffer(buf, maxLineSize)
	return &Parser{
		scanner: *scanner,
		start:   true,
		more:    true,
	}
}

// Next reads next fasta genome in underlying reader and returns the result
// and the amount of bytes read during the call.
// Next only returns an error if it:
//   - Attempts to read and fails to find a valid fasta sequence.
//   - Returns reader's EOF if called after reader has been exhausted.
//   - If a EOF is encountered immediately after a sequence with no newline ending.
//     In this case the Fasta up to that point is returned with an EOF error.
//
// It is worth noting the amount of bytes read are always right up to before
// the next fasta starts which means this function can effectively be used
// to index where fastas start in a file or string.
func (p *Parser) Next() (*Record, error) {
	if !p.more {
		return &Record{}, io.EOF
	}
	for p.scanner.Scan() {
		line := p.scanner.Bytes()
		if p.scanner.Err() != nil {
			break
		}
		p.line++
		switch {
		// if there's nothing on this line skip this iteration of the loop
		case len(line) == 0:
			continue
		// if it's a comment skip this line
		case line[0] == ';':
			continue
		// start of file with no identifier, error
		case line[0] != '>' && p.start:
			err := fmt.Errorf("invalid input: missing sequence identifier for sequence starting at line %d", p.line)
			record, _ := p.newRecord()
			return &record, err
		// start of a fasta line
		case line[0] != '>':
			p.buff.Write(line)
		// Process normal new lines
		case line[0] == '>' && !p.start:
			record, err := p.newRecord()
			// New name
			p.identifier = string(line[1:])
			return &record, err
		// Process first line of file
		case line[0] == '>' && p.start:
			p.identifier = string(line[1:])
			p.start = false
		}
	}
	p.more = false
	// Add final sequence in file
	record, err := p.newRecord()
	if err != nil {
		return &record, err
	}
	return &record, nil
}

func (p *Parser) newRecord() (Record, error) {
	sequence := p.buff.String()
	if sequence == "" {
		return Record{}, fmt.Errorf("%s has no sequence", p.identifier)
	}
	record := Record{
		Identifier: p.identifier,
		Sequence:   sequence,
	}
	// Reset sequence buffer
	p.buff.Reset()
	return record, nil
}

///******************************************************************************
//
//Start of  Write functions
//
//******************************************************************************/

// WriteTo implements the io.WriterTo interface for fasta records.
func (record *Record) WriteTo(w io.Writer) (int64, error) {
	var writtenBytes int64
	var newWrittenBytes int
	newWrittenBytes, err := w.Write([]byte(">"))
	if err != nil {
		return writtenBytes, err
	}
	writtenBytes += int64(newWrittenBytes)
	newWrittenBytes, err = w.Write([]byte(record.Identifier))
	if err != nil {
		return writtenBytes, err
	}
	writtenBytes += int64(newWrittenBytes)
	newWrittenBytes, err = w.Write([]byte("\n"))
	if err != nil {
		return writtenBytes, err
	}
	writtenBytes += int64(newWrittenBytes)

	lineCount := 0
	// write the fasta sequence 80 characters at a time
	for _, character := range record.Sequence {
		newWrittenBytes, err = w.Write([]byte{byte(character)})
		if err != nil {
			return writtenBytes, err
		}
		writtenBytes += int64(newWrittenBytes)
		lineCount++
		if lineCount == 80 {
			newWrittenBytes, err = w.Write([]byte("\n"))
			if err != nil {
				return writtenBytes, err
			}
			writtenBytes += int64(newWrittenBytes)
			lineCount = 0
		}
	}
	newWrittenBytes, err = w.Write([]byte("\n\n"))
	if err != nil {
		return writtenBytes, err
	}
	writtenBytes += int64(newWrittenBytes)
	return writtenBytes, nil
}
