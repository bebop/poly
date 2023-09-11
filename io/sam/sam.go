/*
Package sam implements a SAM file parser and writer.

Paper: https://doi.org/10.1093%2Fbioinformatics%2Fbtp352
Update to do date: http://samtools.github.io/hts-specs/SAMv1.pdf
*/
package sam

import (
	"bufio"
	"fmt"
	"io"
	"strconv"
	"strings"
)

const DefaultMaxLineSize int = 1024 * 32 * 2 // // 32kB is a magic number often used by the Go stdlib for parsing. We multiply it by two.

// Each header in a SAM file begins with an @ followed by a two letter record
// code type. Each line is tab delimited, and contains TAG:VALUE pairs. HD, the
// first line, only occurs once, while SQ, RG, and PG can appear multiple
// times. Finally, @CO contains user generated comments.
//
// For more information, check section 1.3 of the reference document.
type Header struct {
	HD map[string]string   // File-level metadata. Optional. If present, there must be only one @HD line and it must be the first line of the file.
	SQ []map[string]string // Reference sequence dictionary. The order of @SQ lines defines the alignment sorting order.
	RG []map[string]string // Read group. Unordered multiple @RG lines are allowed.
	PG []map[string]string // Program.
	CO []string            // One-line text comment. Unordered multiple @CO lines are allowed. UTF-8 encoding may be used.
}

// WriteTo writes a SAM header to an io.Writer.
func (header *Header) WriteTo(w io.Writer) (int64, error) {
	return 0, nil // TODO
}

// Validate validatest that the header has all required information, as
// described in the SAMv1 specification document.
func (header *Header) Validate() error {
	return nil // TODO
}

// Optional fields in SAM alignments are structured as TAG:TYPE:DATA, where
// the type identifiers the typing of the data.
//
// For more information, check section 1.5 of http://samtools.github.io/hts-specs/SAMv1.pdf.
type Optional struct {
	Type rune   // The type may be one of A (character), B (general array), f (real number), H (hexadecimal array), i (integer), or Z (string).
	Data string // Optional data
}

// Each alignment is a single line of a SAM file, representing a linear
// alignment of a segment, consisting of 11 or more tab delimited fields. The
// 11 fields (QNAME -> QUAL) are always available (if the data isn't there, a
// placeholder '0' or '*' is used instead), with additional optional fields
// following.
//
// For more information, check section 1.4 of the reference document.
type Alignment struct {
	QNAME     string              // Query template NAME
	FLAG      uint16              // bitwise FLAG
	RNAME     string              // References sequence NAME
	POS       int32               // 1- based leftmost mapping POSition
	MAPQ      byte                // MAPping Quality
	CIGAR     string              // CIGAR string
	RNEXT     string              // Ref. name of the mate/next read
	PNEXT     int32               // Position of the mate/next read
	TLEN      int32               // observed Template LENgth
	SEQ       string              // segment SEQuence
	QUAL      string              // ASCII of Phred-scaled base QUALity+33
	Optionals map[string]Optional // Map of TAG to {TYPE:DATA}
}

// Parser is a sam file parser that provide sample control over reading sam
// alignments. It should be initialized with NewParser.
type Parser struct {
	reader     bufio.Reader
	line       uint
	FileHeader Header
	firstLine  string
}

// Header returns the parsed sam header.
func (p *Parser) Header() (*Header, error) {
	return &p.FileHeader, nil
}

// NewParser creates a parser from an io.Reader for sam data. For larger
// alignments, you will want to increase the maxLineSize.
func NewParser(r io.Reader, maxLineSize int) (*Parser, Header, error) {
	parser := &Parser{
		reader: *bufio.NewReaderSize(r, maxLineSize),
	}
	var header Header
	var hdParsed bool
	// Initialize header maps
	header.HD = make(map[string]string)
	header.SQ = []map[string]string{}
	header.RG = []map[string]string{}
	header.PG = []map[string]string{}
	header.CO = []string{}

	// We need to first read the header before returning the parser to the
	// user for analyzing alignments.
	for {
		lineBytes, err := parser.reader.ReadSlice('\n')
		if err != nil {
			return parser, Header{}, err
		}
		line := strings.TrimSpace(string(lineBytes))
		parser.line++
		if len(line) == 0 {
			return parser, Header{}, fmt.Errorf("Line %d is empty. Empty lines are not allowed in headers.", parser.line)
		}
		// If this line is the start of the alignments, set the firstLine
		// into memory, and then break this loop.
		if line[0] != '@' {
			parser.firstLine = line
			break
		}
		values := strings.Split(line, "\t")
		if len(values) < 1 {
			return parser, Header{}, fmt.Errorf("Line %d should contain at least 1 value. Got: %d. Line text: %s", parser.line, len(values), line)
		}

		// If we haven't parsed HD, it is always the first line: lets parse it.
		if !hdParsed {
			if values[0] != "@HD" {
				return parser, Header{}, fmt.Errorf("First line (%d) should always contain @HD first. Line text: %s", parser.line, line)
			}
			// Now parse the rest of the HD header
			for _, value := range values[1:] {
				valueSplit := strings.Split(value, ":")
				header.HD[valueSplit[0]] = valueSplit[1]
			}
			hdParsed = true
			continue
		}

		// CO lines are unique in that they are just strings. So we try to parse them
		// first. We include the entire comment line for these.
		if values[0] == "@CO" {
			header.CO = append(header.CO, line)
			continue
		}

		// HD/CO lines have been successfully parsed, now we work on SQ, RG, and PG.
		// Luckily, each one has an identical form ( TAG:DATA ), so we can parse that
		// first and then just apply it to the respect top level tag.
		genericMap := make(map[string]string)
		for _, value := range values[1:] {
			valueSplit := strings.Split(value, ":")
			genericMap[valueSplit[0]] = valueSplit[1]
		}
		switch values[0] {
		case "@SQ":
			header.SQ = append(header.SQ, genericMap)
		case "@RG":
			header.RG = append(header.RG, genericMap)
		case "@PG":
			header.PG = append(header.PG, genericMap)
		default:
			return parser, Header{}, fmt.Errorf("Line %d should contain @SQ, @RG, @PG or @CO as top level tags, but they weren't found. Line text: %s", parser.line, line)
		}
	}
	parser.FileHeader = header
	return parser, header, nil
}

// Next parsers the next read from a parser. Returns an `io.EOF` upon EOF.
func (parser *Parser) Next() (*Alignment, error) {
	var alignment Alignment
	lineBytes, err := parser.reader.ReadSlice('\n')
	if err != nil {
		return nil, err
	}
	parser.line++
	line := strings.TrimSpace(string(lineBytes))
	values := strings.Split(line, "\t")
	if len(values) < 11 {
		return nil, fmt.Errorf("Line %d had error: must have at least 11 tab-delimited values. Had %d", parser.line, len(values))
	}
	alignment.QNAME = values[0]
	flag64, err := strconv.ParseUint(values[1], 10, 16) // convert string to uint16
	if err != nil {
		return nil, fmt.Errorf("Line %d had error: %s", parser.line, err)
	}
	alignment.FLAG = uint16(flag64)
	alignment.RNAME = values[2]
	pos64, err := strconv.ParseInt(values[3], 10, 32) // convert string to int32
	if err != nil {
		return nil, fmt.Errorf("Line %d had error: %s", parser.line, err)
	}
	alignment.POS = int32(pos64)
	mapq64, err := strconv.ParseUint(values[4], 10, 8) // convert string to uint8 (otherwise known as byte)
	if err != nil {
		return nil, fmt.Errorf("Line %d had error: %s", parser.line, err)
	}
	alignment.MAPQ = uint8(mapq64)
	alignment.CIGAR = values[5]
	alignment.RNEXT = values[6]
	pnext64, err := strconv.ParseInt(values[7], 10, 32)
	if err != nil {
		return nil, fmt.Errorf("Line %d had error: %s", parser.line, err)
	}
	alignment.PNEXT = int32(pnext64)
	tlen64, err := strconv.ParseInt(values[8], 10, 32)
	if err != nil {
		return nil, fmt.Errorf("Line %d had error: %s", parser.line, err)
	}
	alignment.TLEN = int32(tlen64)
	alignment.SEQ = values[9]
	alignment.QUAL = values[10]

	optionals := make(map[string]Optional)
	for _, value := range values[11:] {
		valueSplit := strings.Split(value, ":")
		optionals[valueSplit[0]] = Optional{Type: rune(valueSplit[1][0]), Data: valueSplit[2]}
	}
	alignment.Optionals = optionals
	return &alignment, nil
}
