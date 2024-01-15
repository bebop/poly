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

// A parseContext holds the data structures currently being parsed.
type parseContext struct {
	entry     *Entry
	reference *Reference
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

func (p *Parser) makeWrongContextError(linetype lineType) GenbankSyntaxError {
	switch linetype {
	case keyword:
		return p.makeSyntaxError("unexpected keyword found")
	case empty:
		return p.makeSyntaxError("unexpected empty line found")
	case eof:
		return p.makeSyntaxError("unexpected end of file found")
	case subKeyword:
		return p.makeSyntaxError("subkeyword found outside of keyword context")
	case entryEnd:
		return p.makeSyntaxError("unexpected end of entry found")
	case continuation:
		return p.makeSyntaxError("continuation found outside of keyword context")
	case featureKey:
		return p.makeSyntaxError("feature key found outside of FEATURES keyword context")
	case sequence:
		return p.makeSyntaxError("sequence found outside of ORIGIN keyword context")
	default:
		return p.makeSyntaxError("unknown linetype found, please report this to repository maintainers")
	}
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

func (p *Parser) dispatchByPrefix(pCtx parseContext, funcs map[string]func(parseContext) error) (err error) {
	line, err := p.peekLine()
	if err != nil {
		return err
	}

	for prefix, parseFunc := range funcs {
		if strings.HasPrefix(line, prefix) {
			return parseFunc(pCtx)
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
	keywordFuncs := map[string]func(parseContext) error{
		"LOCUS":      p.parseLocus,
		"DEFINITION": p.parseDefinition,
		"ACCESSION":  p.parseAccession,
		"VERSION":    p.parseVersion,
		"DBLINK":     p.parseDBLink,
		"KEYWORDS":   p.parseKeywords,
		"SEGMENT":    p.parseSegment,
		"SOURCE":     p.parseSource,
	}

	res := Entry{}
	pCtx := parseContext{entry: &res}
	for {
		linetype, err := p.peekNextLineType()
		if err != nil {
			return res, false, fmt.Errorf("could not parse entry: %w", err)
		}

		switch linetype {
		case keyword:
			err = p.dispatchByPrefix(pCtx, keywordFuncs)
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
		case entryEnd:
			return res, false, nil
		default:
			return res, false, p.makeWrongContextError(linetype)
		}
	}
}

/* KEYWORD PARSING FUNCTIONS */

const locusDateLayout = "02-Jan-2006"

// See Genbank spec 3.4.4
func (p *Parser) parseLocus(pCtx parseContext) error {
	line, err := p.readLine()
	if err != nil {
		return fmt.Errorf("failed to read LOCUS line: %w", err)
	}

	re := regexp.MustCompile(`\s+`)
	tokens := re.Split(line, -1)
	if len(tokens) != 8 {
		return fmt.Errorf("LOCUS line should have 8 tokens, got %v", len(tokens))
	}

	pCtx.entry.Name = tokens[1]

	length, err := strconv.Atoi(tokens[2])
	if err != nil {
		return p.makeSyntaxError("failed to parse sequence length", err)
	}
	pCtx.entry.Length = length

	molecule := strings.Split(tokens[4], "-")
	if len(molecule) == 2 {
		pCtx.entry.Strandedness = Strandedness(molecule[0])
		pCtx.entry.MoleculeType = MoleculeType(molecule[1])
	} else {
		pCtx.entry.Strandedness = Unknown
		pCtx.entry.MoleculeType = MoleculeType(molecule[0])
	}

	pCtx.entry.MoleculeToplogy = MoleculeToplogy(tokens[5])

	pCtx.entry.DivisionCode = DivisionCode(tokens[6])

	date, err := time.Parse(locusDateLayout, tokens[7])
	if err != nil {
		return p.makeSyntaxError("failed to parse entry date", err)
	}
	pCtx.entry.UpdateDate = date

	return nil
}

// See Genbank spec 3.4.5
func (p *Parser) parseDefinition(pCtx parseContext) error {
	line, err := p.readLine()
	if err != nil {
		return fmt.Errorf("failed to read DEFINITION line: %w", err)
	}

	after, _ := strings.CutPrefix(line, "DEFINITION")
	after = strings.TrimSpace(after)
	pCtx.entry.Definition = appendLine(pCtx.entry.Definition, after)

	for {
		linetype, err := p.peekNextLineType()
		if err != nil {
			return err
		}

		switch linetype {
		case empty:
			_, err := p.readLine()
			if err != nil {
				return err
			}
			pCtx.entry.Definition = appendLine(pCtx.entry.Definition, "\n")
		case continuation:
			line, err := p.readLine()
			if err != nil {
				return err
			}
			pCtx.entry.Definition = appendLine(pCtx.entry.Definition, line[10:])
		default:
			return nil
		}
	}
}

// See Genbank spec 3.4.6
func (p *Parser) parseAccession(pCtx parseContext) error {
	line, err := p.readLine()
	if err != nil {
		return fmt.Errorf("failed to read ACCESSION line: %w", err)
	}

	after, _ := strings.CutPrefix(line, "ACCESSION")
	after = strings.TrimSpace(after)
	pCtx.entry.Accession += after

	return nil
}

// See Genbank spec 3.4.7.1
func (p *Parser) parseVersion(pCtx parseContext) error {
	line, err := p.readLine()
	if err != nil {
		return fmt.Errorf("failed to read VERSION line: %w", err)
	}

	after, _ := strings.CutPrefix(line, "VERSION")
	accessionVersion := strings.Split(strings.TrimSpace(after), ".")
	if len(accessionVersion) != 2 {
		return p.makeSyntaxError("could not parse accession version")
	}
	version, err := strconv.Atoi(accessionVersion[1])
	if err != nil {
		return p.makeSyntaxError("accession version should be an integer")
	}
	pCtx.entry.AccessionVersion = version

	return nil
}

// See Genbank spec 3.4.7.2
func (p *Parser) parseDBLink(pCtx parseContext) error {
	line, err := p.readLine()
	if err != nil {
		return fmt.Errorf("failed to read DBLINK line: %w", err)
	}

	if pCtx.entry.DatabaseLinks == nil {
		pCtx.entry.DatabaseLinks = make(map[string][]string)
	}

	after, _ := strings.CutPrefix(line, "DBLINK")
	crossRefType, crossRefs, err := parseRawDBLink(after)
	if err != nil {
		return p.makeSyntaxError("failed to parse DBLINK entry", err)
	}
	pCtx.entry.DatabaseLinks[crossRefType] = append(pCtx.entry.DatabaseLinks[crossRefType], crossRefs...)

	for {
		lineType, err := p.peekNextLineType()
		if err == io.EOF {
			return nil
		} else if err != nil {
			return err
		}

		switch lineType {
		case empty:
			_, err := p.readLine()
			if err != nil {
				return err
			}
		case continuation:
			line, err := p.readLine()
			if err != nil {
				return err
			}
			crossRefType, crossRefs, err := parseRawDBLink(line)
			if err != nil {
				return p.makeSyntaxError("failed to parse DBLINK entry", err)
			}
			pCtx.entry.DatabaseLinks[crossRefType] = append(pCtx.entry.DatabaseLinks[crossRefType], crossRefs...)
		default:
			return nil
		}
	}
}

func parseRawDBLink(data string) (crossRefType string, crossRefs []string, err error) {
	before, after, found := strings.Cut(data, ":")
	if !found {
		return "", nil, fmt.Errorf("should be in the format TYPE:REFERENCE")
	}

	crossRefType = strings.TrimSpace(before)
	for _, crossRef := range strings.Split(after, ",") {
		crossRefs = append(crossRefs, strings.TrimSpace(crossRef))
	}

	return
}

// See Genbank spec section 3.4.8
func (p *Parser) parseKeywords(pCtx parseContext) error {
	line, err := p.readLine()
	if err != nil {
		return fmt.Errorf("failed to read KEYWORDS line: %w", err)
	}

	after, _ := strings.CutPrefix(line, "KEYWORDS")
	after, endsWithPeriod := strings.CutSuffix(after, ".")
	if !endsWithPeriod {
		return p.makeSyntaxError("KEYWORDS line must end with a period")
	}
	keywords := strings.Split(after, ";")
	for _, keyword := range keywords {
		pCtx.entry.Keywords = append(pCtx.entry.Keywords, strings.TrimSpace(keyword))
	}

	return nil
}

// See Genbank spec section 3.4.9
func (p *Parser) parseSegment(pCtx parseContext) error {
	line, err := p.readLine()
	if err != nil {
		return fmt.Errorf("failed to read KEYWORDS line: %w", err)
	}

	after, _ := strings.CutPrefix(line, "SEGMENT")
	segmentStr, totalSegmentsStr, found := strings.Cut(after, " of ")
	if !found {
		return p.makeSyntaxError("malformed SEGMENT line, should be in the form 'n of m'")
	}
	segment, err := strconv.Atoi(strings.TrimSpace(segmentStr))
	if err != nil {
		return p.makeSyntaxError("could not parse segment number", err)
	}
	totalSegments, err := strconv.Atoi(strings.TrimSpace(totalSegmentsStr))
	if err != nil {
		return p.makeSyntaxError("could not parse total number of segments", err)
	}

	pCtx.entry.Segment = segment
	pCtx.entry.TotalSegments = totalSegments

	return nil

}

// See Genbank spec section 3.4.10
func (p *Parser) parseSource(pCtx parseContext) error {
	line, err := p.readLine()
	if err != nil {
		return fmt.Errorf("failed to read SOURCE line: %w", err)
	}

	after, _ := strings.CutPrefix(line, "SOURCE")
	after, endsWithPeriod := strings.CutSuffix(strings.TrimSpace(after), ".")
	if !endsWithPeriod {
		return p.makeSyntaxError("SOURCE line must end with a period")
	}
	pCtx.entry.Source.Name = strings.TrimSpace(after)

	return p.parseOrganism(pCtx)

}

// See Genbank spec section 3.4.10
func (p *Parser) parseOrganism(pCtx parseContext) error {
	organismLine, err := p.readLine()
	if err != nil {
		return fmt.Errorf("failed to read ORGANISM line: %w", err)
	}

	after, organismKeyword := strings.CutPrefix(organismLine, "  ORGANISM")
	if !organismKeyword {
		return p.makeSyntaxError("SOURCE keyword must be followed by ORGANISM subkeyword")
	}
	after, endsWithPeriod := strings.CutSuffix(strings.TrimSpace(after), ".")
	if !endsWithPeriod {
		return p.makeSyntaxError("ORGANISM line must end with a period")
	}

	pCtx.entry.Source.ScientificName = strings.TrimSpace(after)

	return p.parseTaxonomy(pCtx)
}

// See Genbank spec section 3.4.10
func (p *Parser) parseTaxonomy(pCtx parseContext) error {
	for {
		lineType, err := p.peekNextLineType()
		if err != nil {
			return err
		}

		switch lineType {
		case empty:
			p.readLine() // skip empty lines
		case continuation:
			line, err := p.readLine()
			if err != nil {
				return err
			}

			lineType, err := p.peekNextLineType()
			if err != nil {
				return err
			}
			if lineType != continuation && !strings.HasSuffix(line, ".") {
				return p.makeSyntaxError("final line of ORGANISM subkeyword taxonomy must end in a period")
			} else if lineType == continuation && !strings.HasSuffix(line, ";") {
				return p.makeSyntaxError("ORGANISM subkeyword taxonomy lines must end in a semicolon")
			}

			taxa := strings.Split(strings.TrimSpace(line[:len(line)-1]), ";")
			for _, taxon := range taxa {
				pCtx.entry.Source.Taxonomy = append(pCtx.entry.Source.Taxonomy, strings.TrimSpace(taxon))
			}
		default:
			return nil
		}

	}
}

// See Genbank spec section 3.4.11
func (p *Parser) parseReference(pCtx parseContext) error {
	line, err := p.readLine()
	if err != nil {
		return err
	}

	after, _ := strings.CutPrefix(line, "REFERENCE")
	after = strings.TrimSpace(after)

	re := regexp.MustCompile(`(\d+)\s+\(bases (\d)+ to (\d+)`)
	matches := re.FindStringSubmatch(after)
	if matches == nil {
		return p.makeSyntaxError("malfored REFERENCE line")
	}

	refNum, _ := strconv.Atoi(matches[0])
	basesFrom, _ := strconv.ParseUint(matches[1], 10, 0)
	basesTo, _ := strconv.ParseUint(matches[2], 10, 0)

	pCtx.entry.References = append(pCtx.entry.References, Reference{
		Number: refNum,
		BaseRange: Range{
			From: uint(basesFrom),
			To:   uint(basesTo),
		},
	})
	pCtx.reference = &pCtx.entry.References[len(pCtx.entry.References)-1]

	subKeywordFuncs := map[string]func(parseContext) error{}

	for {
		lineType, err := p.peekNextLineType()
		if err != nil {
			return fmt.Errorf("could not parse REFERENCE keyword: %w", err)
		}

		switch lineType {
		case subKeyword:
			err = p.dispatchByPrefix(pCtx, subKeywordFuncs)
			if err != nil {
				return err
			}
		case empty:
			p.readLine() // We need to read the empty line to advance the reader

		default:
			pCtx.reference = nil
			return nil
		}
	}
}

// See Genbank spec section 3.4.11
func (p *Parser) parseAuthors(pCtx parseContext) error {
	line, err := p.readLine()
	if err != nil {
		return err
	}

	authors, _ := strings.CutPrefix(line, "AUTHORS")
	authors = strings.TrimSpace(authors)
	for {
		lineType, err := p.peekNextLineType()
		if err != nil {
			return fmt.Errorf("could not parse AUTHORS subkeyword: %w", err)
		}
		switch lineType {
		case continuation:
			line, err := p.readLine()
			if err != nil {
				return err
			}

			authors += strings.TrimSpace(line[:len(line)-1])
		case empty:
			p.readLine() // We need to read the empty line to advance the reader
		default:
			return nil
		}
	}
}

/* OTHER UTILITY FUNCTIONS */

// appendLine appends a line to s. If s is empty, simply
// returns append. If s is not empty, ensures that s is followed
// by a newline and then what is contained in apend.
func appendLine(s string, append string) string {
	if s == "" {
		return append
	}

	return strings.TrimSuffix(s, "\n") + "\n" + append
}
