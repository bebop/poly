/*
Package bio provides utilities for reading and writing sequence data.
*/
package bio

import (
	"bufio"
	"compress/gzip"
	"errors"
	"fmt"
	"io"
	"math"
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
var DefaultMaxLengths = map[Format]int{
	Fasta:   defaultMaxLineLength,
	Fastq:   8 * 1024 * 1024, // The longest single nanopore sequencing read so far is 4Mb. A 8mb buffer should be large enough for any sequencing.
	Gff:     defaultMaxLineLength,
	Genbank: defaultMaxLineLength,
	Slow5:   128 * 1024 * 1024, // 128mb is used because slow5 lines can be massive, since a single read can be many millions of base pairs.
	Pileup:  defaultMaxLineLength,
	Uniprot: defaultMaxLineLength,
	Rebase:  defaultMaxLineLength,
}

/******************************************************************************
Aug 30, 2023

Lower level interfaces

******************************************************************************/

type LowLevelParser[DataType fasta.Record, DataTypeHeader fasta.Header] interface {
	Header() (DataTypeHeader, error)
	Next() (DataType, error)
}

/******************************************************************************

Higher level parse

******************************************************************************/

type Parser[DataType fasta.Record, DataTypeHeader fasta.Header] struct {
	LowLevelParser LowLevelParser[DataType, DataTypeHeader]
}

func NewParser(format Format, r io.Reader) (*Parser[fasta.Record, fasta.Header], error) {
	maxLineLength, ok := DefaultMaxLengths[format]
	if !ok {
		return nil, fmt.Errorf("did not find format in default lengths")
	}
	return NewParserWithMaxLine(format, r, maxLineLength)
}

func NewParserWithMaxLine(format Format, r io.Reader, maxLineLength int) (*Parser[fasta.Record, fasta.Header], error) {
	switch format {
	case Fasta:
		return &Parser[fasta.Record, fasta.Header]{LowLevelParser: fasta.NewParser(r, maxLineLength)}, nil
	}
	return nil, nil
}

func NewParserGz(format Format, r io.Reader) (*Parser[fasta.Record, fasta.Header], error) {
	return NewParserCompressed(format, r, func(r io.Reader) (io.Reader, error) { return gzip.NewReader(r) })
}

func NewParserCompressed(format Format, r io.Reader, decoderFunc func(io.Reader) (io.Reader, error)) (*Parser[fasta.Record, fasta.Header], error) {
	decodedReader, err := decoderFunc(r)
	if err != nil {
		return nil, err
	}
	return NewParser(format, decodedReader)
}

/******************************************************************************

Parser read functions

******************************************************************************/

func Read(format Format, path string) (*Parser[fasta.Record, fasta.Header], error) {
	file, err := os.Open(path)
	if err != nil {
		return nil, err
	}
	return NewParser(format, file)
}

func ReadGz(format Format, path string) (*Parser[fasta.Record, fasta.Header], error) {
	return ReadCompressed(format, path, func(r io.Reader) (io.Reader, error) { return gzip.NewReader(r) })
}

func ReadCompressed(format Format, path string, decoderFunc func(io.Reader) (io.Reader, error)) (*Parser[fasta.Record, fasta.Header], error) {
	file, err := os.Open(path)
	if err != nil {
		return nil, err
	}
	return NewParserCompressed(format, file, decoderFunc)
}

/******************************************************************************

Parser higher-level functions

******************************************************************************/

func (p *Parser[DataType, DataTypeHeader]) Next() (DataType, error) {
	return p.LowLevelParser.Next()
}

func (p *Parser[DataType, DataTypeHeader]) Header() (DataTypeHeader, error) {
	return p.LowLevelParser.Header()
}

func (p *Parser[DataType, DataTypeHeader]) ParseN(countN int) ([]DataType, error) {
	var records []DataType
	for counter := 0; counter < countN; counter++ {
		record, err := p.Next()
		if err != nil {
			if errors.Is(err, io.EOF) {
				err = nil // EOF not treated as parsing error.
			}
			return records, err
		}
		records = append(records, record)
	}
	return records, nil
}

func (p *Parser[DataType, DataTypeHeader]) Parse() ([]DataType, error) {
	return p.ParseN(math.MaxInt)
}

func (p *Parser[DataType, DataTypeHeader]) ParseWithHeader() ([]DataType, DataTypeHeader, error) {
	header, err := p.Header()
	if err != nil {
		return []DataType{}, DataTypeHeader{}, err
	}
	data, err := p.Parse()
	if err != nil {
		return []DataType{}, DataTypeHeader{}, err
	}
	return data, header, nil
}

// Parse to channel parses a
func (p *Parser[DataType, DataTypeHeader]) ParseToChannel(channel chan<- DataType) error {
	for {
		record, err := p.Next()
		if err != nil {
			if errors.Is(err, io.EOF) {
				err = nil // EOF not treated as parsing error.
			}
			return err
		}
		channel <- record
	}
}

/******************************************************************************

Writer functions

******************************************************************************/

func WriteAll(data []io.WriterTo, header io.WriterTo, w io.Writer) error {
	_, err := header.WriteTo(w)
	if err != nil {
		return err
	}
	for _, datum := range data {
		_, err = datum.WriteTo(w)
		if err != nil {
			return err
		}
	}
	return nil
}

func WriteFile(data []io.WriterTo, header io.WriterTo, path string) error {
	file, err := os.Open(path)
	if err != nil {
		return err
	}
	return WriteAll(data, header, file)
}
