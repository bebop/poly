/*
Package bio provides utilities for reading and writing sequence data.
*/
package bio

import (
	"bufio"
	"errors"
	"io"
	"math"

	"github.com/TimothyStiles/poly/bio/fasta"
	"github.com/TimothyStiles/poly/bio/fastq"
	"github.com/TimothyStiles/poly/bio/pileup"
	"github.com/TimothyStiles/poly/bio/slow5"
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
}

/******************************************************************************
Aug 30, 2023

Lower level interfaces

******************************************************************************/

type LowLevelParser[DataType fasta.Record | fastq.Read | slow5.Read | pileup.Line, DataTypeHeader fasta.Header | fastq.Header | slow5.Header | pileup.Header] interface {
	Header() (*DataTypeHeader, error)
	Next() (*DataType, error)
}

/******************************************************************************

Higher level parse

******************************************************************************/

type Parser[DataType fasta.Record | fastq.Read | slow5.Read | pileup.Line, DataTypeHeader fasta.Header | fastq.Header | slow5.Header | pileup.Header] struct {
	LowLevelParser LowLevelParser[DataType, DataTypeHeader]
}

func NewFastaParser(r io.Reader) (*Parser[fasta.Record, fasta.Header], error) {
	return NewFastaParserWithMaxLineLength(r, DefaultMaxLengths[Fasta])
}

func NewFastaParserWithMaxLineLength(r io.Reader, maxLineLength int) (*Parser[fasta.Record, fasta.Header], error) {
	return &Parser[fasta.Record, fasta.Header]{LowLevelParser: fasta.NewParser(r, maxLineLength)}, nil
}

func NewFastqParser(r io.Reader) (*Parser[fastq.Read, fastq.Header], error) {
	return NewFastqParserWithMaxLineLength(r, DefaultMaxLengths[Fastq])
}

func NewFastqParserWithMaxLineLength(r io.Reader, maxLineLength int) (*Parser[fastq.Read, fastq.Header], error) {
	return &Parser[fastq.Read, fastq.Header]{LowLevelParser: fastq.NewParser(r, maxLineLength)}, nil
}

func NewSlow5Parser(r io.Reader) (*Parser[slow5.Read, slow5.Header], error) {
	return NewSlow5ParserWithMaxLineLength(r, DefaultMaxLengths[Slow5])
}

func NewSlow5ParserWithMaxLineLength(r io.Reader, maxLineLength int) (*Parser[slow5.Read, slow5.Header], error) {
	parser, err := slow5.NewParser(r, maxLineLength)
	return &Parser[slow5.Read, slow5.Header]{LowLevelParser: parser}, err
}

func NewPileupParser(r io.Reader) (*Parser[pileup.Line, pileup.Header], error) {
	return NewPileupParserWithMaxLineLength(r, DefaultMaxLengths[Pileup])
}

func NewPileupParserWithMaxLineLength(r io.Reader, maxLineLength int) (*Parser[pileup.Line, pileup.Header], error) {
	return &Parser[pileup.Line, pileup.Header]{LowLevelParser: pileup.NewParser(r, maxLineLength)}, nil
}

/******************************************************************************

Parser higher-level functions

******************************************************************************/

func (p *Parser[DataType, DataTypeHeader]) Next() (*DataType, error) {
	return p.LowLevelParser.Next()
}

func (p *Parser[DataType, DataTypeHeader]) Header() (*DataTypeHeader, error) {
	return p.LowLevelParser.Header()
}

func (p *Parser[DataType, DataTypeHeader]) ParseN(countN int) ([]*DataType, error) {
	var records []*DataType
	for counter := 0; counter < countN; counter++ {
		record, err := p.Next()
		if err != nil {
			if errors.Is(err, io.EOF) {
				err = nil // EOF not treated as parsing error.
			}
			if record != nil {
				records = append(records, record)
			}
			return records, err
		}
		records = append(records, record)
	}
	return records, nil
}

func (p *Parser[DataType, DataTypeHeader]) Parse() ([]*DataType, error) {
	return p.ParseN(math.MaxInt)
}

func (p *Parser[DataType, DataTypeHeader]) ParseWithHeader() ([]*DataType, *DataTypeHeader, error) {
	header, headerErr := p.Header()
	data, err := p.Parse()
	if headerErr != nil {
		return data, header, err
	}
	if err != nil {
		return data, header, err
	}
	return data, header, nil
}

func (p *Parser[DataType, DataTypeHeader]) ParseToChannel(channel chan<- DataType) error {
	for {
		record, err := p.Next()
		if err != nil {
			if errors.Is(err, io.EOF) {
				err = nil // EOF not treated as parsing error.
			}
			close(channel)
			return err
		}
		channel <- *record
	}
}
