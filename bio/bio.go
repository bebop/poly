/*
Package bio provides utilities for reading and writing sequence data.
*/
package bio

import (
	"bufio"
	"io"
	"os"

	"github.com/TimothyStiles/poly/bio/fasta"
)

type Format int

const (
	Fasta Format = iota
	Fastq
	Gff
	Genbank
	Slow5
	Pileup
	Uniprot
	Rebase
)

// DefaultMaxLineLength variables are defined for performance reasons. While
// parsing, reading byte-by-byte takes far, far longer than reading many bytes
// into a buffer. In golang, this buffer in bufio is usually 64kb. However,
// many files, especially slow5 files, can contain single lines that are much
// longer. We use the default maxLineLength from bufio unless there is a
// particular reason to use a different number.
const defaultMaxLineLength int = bufio.MaxScanTokenSize // 64kB is a magic number often used by the Go stdlib for parsing.
var (
	Fasta_DefaultMaxLineLength   = defaultMaxLineLength
	Fastq_DefaultMaxLineLength   = 8 * 1024 * 1024 // The longest single nanopore sequencing read so far is 4Mb. A 8mb buffer should be large enough for any sequencing.
	Gff_DefaultMaxLineLength     = defaultMaxLineLength
	Genbank_DefaultMaxLineLength = defaultMaxLineLength
	Slow5_DefaultMaxLineLength   = 128 * 1024 * 1024 // 128mb is used because slow5 lines can be massive, since a single read can be many millions of base pairs.
	Pileup_DefaultMaxLineLength  = defaultMaxLineLength
	Uniprot_DefaultMaxLineLength = defaultMaxLineLength
	Rebase_DefaultMaxLineLength  = defaultMaxLineLength
)

/******************************************************************************
Aug 30, 2023

Lower level interfaces

******************************************************************************/

type LowLevelParser[DataType fasta.Fasta, DataTypeHeader fasta.Header] interface {
	Header() (DataTypeHeader, int64, error)
	Next() (DataType, int64, error)
}

type Writer interface {
	Write(io.Writer) (int, error)
}

/******************************************************************************

Higher level parse

******************************************************************************/

type Parser[DataType fasta.Fasta, DataTypeHeader fasta.Header] struct {
	LowLevelParser LowLevelParser[DataType, DataTypeHeader]
}

var emptyParser Parser[fasta.Fasta, fasta.Header] = Parser[fasta.Fasta, fasta.Header]{}

func NewParser(format Format, r io.Reader) (*Parser[fasta.Fasta, fasta.Header], error) {
	var maxLineLength int
	switch format {
	case Fasta:
		maxLineLength = Fasta_DefaultMaxLineLength
	}
	return NewParserWithMaxLine(format, r, maxLineLength)
}

func NewParserWithMaxLine(format Format, r io.Reader, maxLineLength int) (*Parser[fasta.Fasta, fasta.Header], error) {
	switch format {
	case Fasta:
		return &Parser[fasta.Fasta, fasta.Header]{LowLevelParser: fasta.NewParser(r, maxLineLength)}, nil
	}
	return nil, nil
}

/******************************************************************************

Parser read functions

******************************************************************************/

func ReadWithMaxLine(format Format, path string, maxLineLength int) (*Parser[fasta.Fasta, fasta.Header], error) {
	file, err := os.Open(path)
	if err != nil {
		return &emptyParser, err
	}
	return NewParserWithMaxLine(format, file, maxLineLength)
}

func Read(format Format, path string) (*Parser[fasta.Fasta, fasta.Header], error) {
	file, err := os.Open(path)
	if err != nil {
		return &emptyParser, err
	}
	return NewParser(format, file)
}

/******************************************************************************

Parser higher-level functions

******************************************************************************/

func (p *Parser[DataType, DataTypeHeader]) Next() (DataType, int64, error) {
	return p.LowLevelParser.Next()
}

func (p *Parser[DataType, DataTypeHeader]) Header() (DataTypeHeader, int64, error) {
	return p.LowLevelParser.Header()
}

func (p *Parser[DataType, DataTypeHeader]) ParseN(countN int) ([]DataType, error) {
	return nil, nil
}

func (p *Parser[DataType, DataTypeHeader]) Parse() ([]DataType, error) {
	return nil, nil
}

func (p *Parser[DataType, DataTypeHeader]) ParseAll() ([]DataType, DataTypeHeader, error) {
	header, _, err := p.Header()
	if err != nil {
		return []DataType{}, DataTypeHeader{}, err
	}
	data, err := p.Parse()
	if err != nil {
		return []DataType{}, DataTypeHeader{}, err
	}
	return data, header, nil
}

func (p *Parser[DataType, DataTypeHeader]) ParseConcurrent(channel chan<- DataType) error {
	return nil
}

/******************************************************************************

Writer functions

******************************************************************************/

func WriteAll(data []Writer, header Writer, w io.Writer) error {
	_, err := header.Write(w)
	if err != nil {
		return err
	}
	for _, datum := range data {
		_, err = datum.Write(w)
		if err != nil {
			return err
		}
	}
	return nil
}

func WriteFile(data []Writer, header Writer, path string) error {
	file, err := os.Open(path)
	if err != nil {
		return err
	}
	return WriteAll(data, header, file)
}
