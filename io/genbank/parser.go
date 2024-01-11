package genbank

import (
	"bufio"
	"fmt"
	"io"
	"regexp"
	"strconv"
	"strings"
	"time"
	"unicode"
)

/* PARSER TYPE DEFINITIONS */

// A Parser stores state during parsing of a Genbank
// flatfile.
type Parser struct {
	reader     *bufio.Reader
	line       uint
	currLine   string
	peekedLine string
}

// A context holds the data structures currently being parsed.
type context struct {
	entry *Entry
}

// NewParser instantiates a new Genbank parser to
// parse a flatfile in an io.Reader.
func NewParser(reader io.Reader) *Parser {
	return &Parser{
		reader: bufio.NewReader(reader),
		line:   0,
	}
}

/* PARSER UTILITY FUNCTIONS */

func (p *Parser) makeSyntaxError(msg string, innerError ...error) GenbankSyntaxError {
	res := GenbankSyntaxError{
		Msg:     msg,
		Line:    p.line,
		Context: p.currLine,
	}
	if len(innerError) > 0 {
		res.InnerErr = innerError[0]
	}

	return res
}

// readLine reads a line from the underlying reader.
// Does not return a line's newline.
func (p *Parser) readLine() (string, error) {
	// If we previously peeked a line, make sure we output it
	// before we start reading new lines. Otherwise, read
	// from the underlying reader.
	var res string
	if p.peekedLine != "" {
		res = p.peekedLine
	} else {
		var err error
		res, err = p.reader.ReadString('\n')
		if err != nil && !(err == io.EOF && res != "") { // Don't return an error if we are reading the last line
			return "", err
		}
	}

	p.peekedLine = ""
	p.line++
	p.currLine = res
	trimmed, _ := strings.CutSuffix(res, "\n")
	return trimmed, nil
}

// peekLine returns the next line without consuming it from
// the underlying reader. Repeated calls to peekLine return
// the same result. Does not return a line's newline.
func (p *Parser) peekLine() (string, error) {
	// If we've already peeked this line, simply return it again.
	if p.peekedLine != "" {
		res, _ := strings.CutSuffix(p.peekedLine, "\n")
		return res, nil
	}

	res, err := p.reader.ReadString('\n')
	if err != nil && !(err == io.EOF && res != "") { // Don't return an error if we are reading the last line
		return "", err
	}

	p.peekedLine = res
	res, _ = strings.CutSuffix(res, "\n")
	return res, nil
}

type lineType int

const (
	empty lineType = iota
	subKeyword
	continuation
	featureKey
	sequence
	entryEnd
	keyword
	eof
)

// Determine which type of information the next line encodes
// (see Genbank spec 3.4.2).
func (p *Parser) peekNextLineType() (lineType, error) {
	nextLine, err := p.peekLine()
	if err == io.EOF {
		return eof, nil
	} else if err != nil {
		return empty, err
	}

	switch {
	case len(strings.TrimSpace(nextLine)) == 0:
		return empty, nil
	case unicode.IsLetter(rune(nextLine[0])):
		return keyword, nil
	case strings.HasPrefix(nextLine, "//"):
		return entryEnd, nil
	case unicode.IsLetter(rune(nextLine[2])):
		return subKeyword, nil
	case !unicode.IsSpace(rune(nextLine[5])):
		return featureKey, nil
	}

	if _, err := strconv.Atoi(strings.TrimSpace(nextLine[:9])); err == nil {
		return sequence, nil
	} else if len(nextLine) > 0 {
		return continuation, nil
	}

	return empty, fmt.Errorf("unrecognized line type")
}

func (p *Parser) dispatchByPrefix(ctx context, funcs map[string]func(context) error) (err error) {
	line, err := p.peekLine()
	if err != nil {
		return err
	}

	for prefix, parseFunc := range funcs {
		if strings.HasPrefix(line, prefix) {
			return parseFunc(ctx)
		}
	}

	p.readLine() // We only peeked previously, so we have to read for proper error reporting
	return p.makeSyntaxError("keyword expected, none found")
}

/* MAIN PARSER FUNCTIONS */

// Parse parses a Genbank file from a Parser's
// io.Reader.
func (p *Parser) Parse() (Genbank, error) {
	res := Genbank{}

	header, err := p.parseHeader()
	if err != nil {
		return res, fmt.Errorf("failed to parse header: %w", err)
	}
	res.Header = header

	res.Entries = make([]Entry, 0)
	var entry Entry
	reachedEOF := false
	for !reachedEOF {
		entry, reachedEOF, err = p.parseEntry()
		if err != nil {
			return res, fmt.Errorf("failed to parse entry %v: %w", len(res.Entries), err)
		}
		res.Entries = append(res.Entries, entry)

	}

	return res, nil
}

const headerDateLayout = "January 2 2006"

func (p *Parser) parseHeader() (Header, error) {
	res := Header{}

	fileNameLine, err := p.readLine()
	if err != nil {
		return res, fmt.Errorf("failed to read file name: %w", err)
	}
	fileName := strings.TrimSpace(fileNameLine[:10])
	if fileName == "" {
		return res, p.makeSyntaxError("empty file name")
	}
	res.FileName = fileName

	dateLine, err := p.readLine()
	if err != nil {
		return res, fmt.Errorf("failed to read date: %w", err)
	}
	date, err := time.Parse(headerDateLayout, strings.TrimSpace(dateLine))
	if err != nil {
		return res, p.makeSyntaxError("failed to parse date", err)
	}
	res.Date = date

	p.readLine() // Skip empty line

	releaseLine, err := p.readLine()
	if err != nil {
		return res, fmt.Errorf("failed to read release: %w", err)
	}
	if len(releaseLine) < 48 {
		return res, p.makeSyntaxError("release line too short")
	}
	release := strings.Split(strings.TrimSpace(releaseLine[47:]), ".")
	if len(release) != 2 {
		return res, p.makeSyntaxError("malformed release version")
	}
	res.MajorRelease, err = strconv.Atoi(release[0])
	if err != nil {
		return res, p.makeSyntaxError("failed to parse major release", err)
	}
	res.MinorRelease, err = strconv.Atoi(release[1])
	if err != nil {
		return res, p.makeSyntaxError("failed to parse minor release", err)
	}

	p.readLine() // Skip empty line

	titleLine, err := p.readLine()
	if err != nil {
		return res, fmt.Errorf("failed to read title: %w", err)
	}
	res.Title = strings.TrimSpace(titleLine)

	p.readLine() // Skip empty line

	statsLine, err := p.readLine()
	if err != nil {
		return res, fmt.Errorf("failed to read file statistics: %w", err)
	}
	if len(releaseLine) < 48 {
		return res, p.makeSyntaxError("release line too short")
	}
	res.NumEntries, err = strconv.Atoi(strings.TrimSpace(statsLine[:8]))
	if err != nil {
		return res, p.makeSyntaxError("failed to parse number of entries/loci", err)
	}
	res.NumBases, err = strconv.Atoi(strings.TrimSpace(statsLine[15:26]))
	if err != nil {
		return res, p.makeSyntaxError("failed to parse number of bases", err)
	}
	res.NumSequences, err = strconv.Atoi(strings.TrimSpace(statsLine[39:47]))
	if err != nil {
		return res, p.makeSyntaxError("failed to parse number of sequences", err)
	}

	p.readLine() // Skip empty line

	return res, nil
}

func (p *Parser) parseEntry() (entry Entry, reachedEOF bool, err error) {
	keywordFuncs := map[string]func(context) error{
		"LOCUS": p.parseLocus,
	}

	res := Entry{}
	ctx := context{entry: &res}
	for {
		lineType, err := p.peekNextLineType()
		if err != nil {
			return res, false, fmt.Errorf("could not parse entry: %w", err)
		}

		switch lineType {
		case keyword:
			err = p.dispatchByPrefix(ctx, keywordFuncs)
			if err == io.EOF {
				return res, true, nil
			} else if err != nil {
				return res, false, err
			} else if _, err := p.peekLine(); err == io.EOF {
				return res, true, nil
			}

		case empty:
			p.readLine() // We need to read the empty line to advance the reader
		case eof:
			return res, true, nil
		case subKeyword:
			return res, false, p.makeSyntaxError("subkeyword found outside of keyword context")
		case entryEnd:
			return res, false, nil
		case continuation:
			return res, false, p.makeSyntaxError("continuation found outside of keyword context")
		case featureKey:
			return res, false, p.makeSyntaxError("feature key found outside of FEATURES keyword context")
		case sequence:
			return res, false, p.makeSyntaxError("sequence found outside of ORIGIN keyword context")
		}
	}
}

/* KEYWORD PARSING FUNCTIONS */

const locusDateLayout = "02-Jan-2006"

// See Genbank spec 3.4.4
func (p *Parser) parseLocus(ctx context) error {
	line, err := p.readLine()
	if err != nil {
		return fmt.Errorf("failed to read LOCUS line: %w", err)
	}

	re := regexp.MustCompile(`\s+`)
	tokens := re.Split(line, -1)
	if len(tokens) != 8 {
		return fmt.Errorf("LOCUS line should have 8 tokens, got %v", len(tokens))
	}

	ctx.entry.Name = tokens[1]

	length, err := strconv.Atoi(tokens[2])
	if err != nil {
		return p.makeSyntaxError("failed to parse sequence length", err)
	}
	ctx.entry.Length = length

	molecule := strings.Split(tokens[4], "-")
	if len(molecule) == 2 {
		ctx.entry.Strandedness = Strandedness(molecule[0])
		ctx.entry.MoleculeType = MoleculeType(molecule[1])
	} else {
		ctx.entry.Strandedness = Unknown
		ctx.entry.MoleculeType = MoleculeType(molecule[0])
	}

	ctx.entry.MoleculeToplogy = MoleculeToplogy(tokens[5])

	ctx.entry.DivisionCode = DivisionCode(tokens[6])

	date, err := time.Parse(locusDateLayout, tokens[7])
	if err != nil {
		return p.makeSyntaxError("failed to parse entry date", err)
	}
	ctx.entry.UpdateDate = date

	return nil
}
