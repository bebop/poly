/*
Package bio provides utilities for reading and writing sequence data.

There are many different biological file formats for different applications.
The poly bio package provides a consistent way to work with each of these file
formats. The underlying data returned by each parser is as raw as we can return
while still being easy to use for downstream applications, and should be
immediately recognizable as the original format.
*/
package bio

import (
	"bufio"
	"context"
	"errors"
	"io"
	"math"

	"github.com/TimothyStiles/poly/bio/fasta"
	"github.com/TimothyStiles/poly/bio/fastq"
	"github.com/TimothyStiles/poly/bio/genbank"
	"github.com/TimothyStiles/poly/bio/pileup"
	"github.com/TimothyStiles/poly/bio/slow5"
	"golang.org/x/sync/errgroup"
)

// Format is a enum of different parser formats.
type Format int

const (
	Fasta Format = iota
	Fastq
	Genbank
	Slow5
	Pileup
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
	Genbank: defaultMaxLineLength,
	Slow5:   128 * 1024 * 1024, // 128mb is used because slow5 lines can be massive, since a single read can be many millions of base pairs.
	Pileup:  defaultMaxLineLength,
}

/******************************************************************************
Aug 30, 2023

Lower level interfaces

******************************************************************************/

// parserInterface is a generic interface that all parsers must support. It is
// very simple, only requiring two functions, Header() and Next(). Header()
// returns the header of the file if there is one: most files, like fasta,
// fastq, and pileup do not contain headers, while others like sam and slow5 do
// have headers. Next() returns a record/read/line from the file format, and
// terminates on an io.EOF error.
//
// Next() may terminate with an io.EOF error with a nil Data or with a
// full Data, depending on where the EOF is in the actual file. A check
// for this is needed at the last Next(), when it returns an io.EOF error. A
// pointer is used to represent the difference between a null Data and an
// empty Data.
type parserInterface[Data io.WriterTo, Header io.WriterTo] interface {
	Header() (Header, error)
	Next() (Data, error)
}

/******************************************************************************

Higher level parse

******************************************************************************/

// Parser is generic bioinformatics file parser. It contains a LowerLevelParser
// and implements useful functions on top of it: such as Parse(), ParseToChannel(), and
// ParseWithHeader().
type Parser[Data io.WriterTo, Header io.WriterTo] struct {
	parserInterface parserInterface[Data, Header]
}

// NewFastaParser initiates a new FASTA parser from an io.Reader.
func NewFastaParser(r io.Reader) (*Parser[*fasta.Record, *fasta.Header], error) {
	return NewFastaParserWithMaxLineLength(r, DefaultMaxLengths[Fasta])
}

// NewFastaParserWithMaxLineLength initiates a new FASTA parser from an
// io.Reader and a user-given maxLineLength.
func NewFastaParserWithMaxLineLength(r io.Reader, maxLineLength int) (*Parser[*fasta.Record, *fasta.Header], error) {
	return &Parser[*fasta.Record, *fasta.Header]{parserInterface: fasta.NewParser(r, maxLineLength)}, nil
}

// NewFastqParser initiates a new FASTQ parser from an io.Reader.
func NewFastqParser(r io.Reader) (*Parser[*fastq.Read, *fastq.Header], error) {
	return NewFastqParserWithMaxLineLength(r, DefaultMaxLengths[Fastq])
}

// NewFastqParserWithMaxLineLength initiates a new FASTQ parser from an
// io.Reader and a user-given maxLineLength.
func NewFastqParserWithMaxLineLength(r io.Reader, maxLineLength int) (*Parser[*fastq.Read, *fastq.Header], error) {
	return &Parser[*fastq.Read, *fastq.Header]{parserInterface: fastq.NewParser(r, maxLineLength)}, nil
}

// NewGenbankParser initiates a new Genbank parser form an io.Reader.
func NewGenbankParser(r io.Reader) (*Parser[*genbank.Genbank, *genbank.Header], error) {
	return NewGenbankParserWithMaxLineLength(r, DefaultMaxLengths[Genbank])
}

// NewGenbankParserWithMaxLineLength initiates a new Genbank parser from an
// io.Reader and a user-given maxLineLength.
func NewGenbankParserWithMaxLineLength(r io.Reader, maxLineLength int) (*Parser[*genbank.Genbank, *genbank.Header], error) {
	return &Parser[*genbank.Genbank, *genbank.Header]{parserInterface: genbank.NewParser(r, maxLineLength)}, nil
}

// NewSlow5Parser initiates a new SLOW5 parser from an io.Reader.
func NewSlow5Parser(r io.Reader) (*Parser[*slow5.Read, *slow5.Header], error) {
	return NewSlow5ParserWithMaxLineLength(r, DefaultMaxLengths[Slow5])
}

// NewSlow5ParserWithMaxLineLength initiates a new SLOW5 parser from an
// io.Reader and a user-given maxLineLength.
func NewSlow5ParserWithMaxLineLength(r io.Reader, maxLineLength int) (*Parser[*slow5.Read, *slow5.Header], error) {
	parser, err := slow5.NewParser(r, maxLineLength)
	return &Parser[*slow5.Read, *slow5.Header]{parserInterface: parser}, err
}

// NewPileupParser initiates a new Pileup parser from an io.Reader.
func NewPileupParser(r io.Reader) (*Parser[*pileup.Line, *pileup.Header], error) {
	return NewPileupParserWithMaxLineLength(r, DefaultMaxLengths[Pileup])
}

// NewPileupParserWithMaxLineLength initiates a new Pileup parser from an
// io.Reader and a user-given maxLineLength.
func NewPileupParserWithMaxLineLength(r io.Reader, maxLineLength int) (*Parser[*pileup.Line, *pileup.Header], error) {
	return &Parser[*pileup.Line, *pileup.Header]{parserInterface: pileup.NewParser(r, maxLineLength)}, nil
}

/******************************************************************************

Parser higher-level functions

******************************************************************************/

// Next is a parsing primitive that should be used when low-level control is
// needed. It returns the next record/read/line from the parser. On EOF, it
// returns an io.EOF error, though the returned FASTA/FASTQ/SLOW5/Pileup may or
// may not be nil, depending on where the io.EOF is. This should be checked by
// downstream software. Next can only be called as many times as there are
// records in a file, as the parser reads the underlying io.Reader in a
// straight line.
func (p *Parser[Data, Header]) Next() (Data, error) {
	return p.parserInterface.Next()
}

// Header is a parsing primitive that should be used when low-level control is
// needed. It returns the header of the parser, which is usually parsed prior
// to the parser being returned by the "NewXXXParser" functions. Unlike the
// Next() function, Header() can be called as many times as needed. Sometimes
// files have useful headers, while other times they do not.
//
// The following file formats do not have a useful header:
//
//	FASTA
//	FASTQ
//	Pileup
//
// The following file formats do have a useful header:
//
//	SLOW5
func (p *Parser[Data, Header]) Header() (Header, error) {
	return p.parserInterface.Header()
}

// ParseN returns a countN number of records/reads/lines from the parser.
func (p *Parser[Data, Header]) ParseN(countN int) ([]Data, error) {
	var records []Data
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

// Parse returns all records/reads/lines from the parser, but does not include
// the header. It can only be called once on a given parser because it will
// read all the input from the underlying io.Reader before exiting.
func (p *Parser[Data, Header]) Parse() ([]Data, error) {
	return p.ParseN(math.MaxInt)
}

// ParseWithHeader returns all records/reads/lines, plus the header, from the
// parser. It can only be called once on a given parser because it will read
// all the input from the underlying io.Reader before exiting.
func (p *Parser[Data, Header]) ParseWithHeader() ([]Data, Header, error) {
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

/******************************************************************************

Concurrent higher-level functions

******************************************************************************/

// ParseToChannel pipes all records/reads/lines from a parser into a channel,
// then optionally closes that channel. If parsing a single file,
// "keepChannelOpen" should be set to false, which will close the channel once
// parsing is complete. If many files are being parsed to a single channel,
// keepChannelOpen should be set to true, so that an external function will
// close channel once all are done parsing.
//
// Context can be used to close the parser in the middle of parsing - for
// example, if an error is found in another parser elsewhere and all files
// need to close.
func (p *Parser[Data, Header]) ParseToChannel(ctx context.Context, channel chan<- Data, keepChannelOpen bool) error {
	for {
		select {
		case <-ctx.Done():
			return ctx.Err()
		default:
			record, err := p.Next()
			if err != nil {
				if errors.Is(err, io.EOF) {
					err = nil // EOF not treated as parsing error.
				}
				if !keepChannelOpen {
					close(channel)
				}
				return err
			}
			channel <- record
		}
	}
}

// ManyToChannel is a generic function that implements the ManyXXXToChannel
// functions. It properly does concurrent parsing of many parsers to a
// single channel, then closes that channel. If any of the files fail to
// parse, the entire pipeline exits and returns.
func ManyToChannel[Data io.WriterTo, Header io.WriterTo](ctx context.Context, channel chan<- Data, parsers ...*Parser[Data, Header]) error {
	errorGroup, ctx := errgroup.WithContext(ctx)
	// For each parser, start a new goroutine to parse data to the channel
	for _, p := range parsers {
		parser := p // Copy to local variable to avoid loop variable scope issues
		errorGroup.Go(func() error {
			// Use the parser's ParseToChannel function, but keep the channel open
			return parser.ParseToChannel(ctx, channel, true)
		})
	}
	// Wait for all parsers to complete
	err := errorGroup.Wait()
	close(channel)
	return err
}
