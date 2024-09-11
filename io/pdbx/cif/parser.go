package cif

import (
	"bufio"
	"fmt"
	"io"
	"strconv"
	"strings"
)

const whitespaceChars = " \t\n"

// A Parser parses CIF data from an io.Reader.
type Parser struct {
	reader *bufio.Reader

	// State held during parsing.
	line           int
	dataBlockName  string
	saveFrameName  string
	cif            CIF
	lastByteWasEOL bool
}

// NewParser creates a new Parser from an io.Reader.
func NewParser(r io.Reader) *Parser {
	return &Parser{
		reader: bufio.NewReader(r),
	}
}

// Parse parses data in io.Reader the Parser was provided into a CIF.
// Stops parsing when an io.EOF is encountered.
func (p *Parser) Parse() (CIF, error) {
	// Clean up on failure.
	defer func() {
		p.dataBlockName = ""
		p.saveFrameName = ""
	}()

	p.cif = NewCIF()

	// Handle until an error is found.
	for {
		err := p.peekAndHandle(
			map[string]func() error{
				// Skip comments.
				"#": p.skipComment,

				// Handle whitespace and newlines.
				" ":  p.skipWhitespace,
				"\t": p.skipWhitespace,
				"\n": p.skipWhitespace,

				// Handle data headers.
				"data_": p.handleDataBlockHeader,
				"save_": p.handleSaveFrameHeader,

				// Handle data items.
				"loop_": p.handleLoop,
				"_":     p.handleTagValue,
			},
			func() error {
				// EOFs aren't real errors
				if _, err := p.reader.ReadByte(); err == io.EOF {
					return io.EOF
				}

				return p.makeSyntaxError("unrecognized token")
			},
		)

		// Swallow EOF errors.
		if err == io.EOF {
			if p.saveFrameName != "" {
				return p.cif, p.makeSyntaxError("save frame %q was not terminated before EOF", p.saveFrameName)
			}

			return p.cif, nil
		} else if err != nil {
			return p.cif, err
		}
	}
}

/* ----------- Handler functions -------------- */

// handleDataBlockHeader reads in a data block header
// and sets the current data block name.
func (p *Parser) handleDataBlockHeader() error {
	header, err := p.readUntilWhitespace()
	if err != nil {
		return err
	}

	name := header[len("data_"):]

	if len(name) == 0 {
		return p.makeSyntaxError("data block header missing name")
	}

	if p.saveFrameName != "" {
		return p.makeSyntaxError("save frame %q was not terminated before next data block header", p.saveFrameName)
	}

	if _, exists := p.cif.DataBlocks[name]; exists {
		return p.makeSyntaxError(fmt.Sprintf("data block with name %q already encountered", name))
	}

	p.dataBlockName = name
	p.cif.DataBlocks[p.dataBlockName] = NewDataBlock(p.dataBlockName)
	return nil
}

// handleDataBlockHeader reads in a data block header
// and sets the current save frame name.
func (p *Parser) handleSaveFrameHeader() error {
	if p.dataBlockName == "" {
		return p.makeSyntaxError("save frame header found before data block header")
	}

	header, err := p.readUntilWhitespace()
	if err != nil {
		return err
	}

	// Handle closing of a save frame.
	if header == "save_" {
		if p.saveFrameName == "" {
			return p.makeSyntaxError("save frames must be named")
		}
		p.saveFrameName = ""
		return nil
	}

	name := header[len("save_"):]

	if _, exists := p.cif.DataBlocks[p.dataBlockName].SaveFrames[name]; exists {
		return p.makeSyntaxError("save frame with name %q already encountered in data block %q", name, p.dataBlockName)
	}

	p.saveFrameName = name
	p.cif.DataBlocks[p.dataBlockName].SaveFrames[p.saveFrameName] = NewSaveFrame(p.saveFrameName)

	return nil
}

// handleTagValue handles a tag:value data item.
func (p *Parser) handleTagValue() error {
	// Ensure we are in a data block.
	if p.dataBlockName == "" {
		return p.makeSyntaxError("tag:value pairs can only exist within a data block")
	}

	tag, err := p.readTag()
	if err != nil {
		if err, ok := err.(CIFSyntaxError); ok {
			return err.Wrap("could not read tag of tag:value pair")
		}
	}

	p.skipWhitespace()
	if err != nil {
		if err, ok := err.(CIFSyntaxError); ok {
			return err.Wrap("tag of tag:value pair must be followed by whitespace")
		}
	}

	val, err := p.readValue()
	if err != nil {
		if err, ok := err.(CIFSyntaxError); ok {
			return err.Wrap("could not read value of tag:value pair")
		}
	}

	dataItems := p.currDataItems()
	if _, exists := dataItems[tag]; exists {
		return p.makeSyntaxError("tag %q already has a value (%v)", tag, dataItems[tag])
	}

	dataItems[tag] = val
	return nil
}

// handleLoop handles a loop_ entry.
func (p *Parser) handleLoop() error {
	p.reader.Discard(len([]byte("loop_")))

	// Ensure we are in a data block.
	if p.dataBlockName == "" {
		return p.makeSyntaxError("loops can only exist within a data block")
	}

	// Ensure there is nothing after the loop_ token.
	if err := p.skipWhitespace(); err != nil {
		return p.makeSyntaxError("loop_ must be followed by whitespace")
	}

	// Read in tags until we find something that isn't a tag.
	tags := make([]string, 0)
	valueFound := false
	for !valueFound {
		err := p.peekAndHandle(
			map[string]func() error{
				"_": func() error {
					tag, err := p.readTag()
					if err != nil {
						return err
					}

					// Tags must be followed by whitespace
					if err := p.skipWhitespace(); err != nil {
						return err
					}

					tags = append(tags, tag)
					return nil
				},
			},
			func() error {
				valueFound = true
				return nil
			},
		)
		if err != nil {
			if err, ok := err.(CIFSyntaxError); ok {
				return err.Wrap("could not read tags in loop_")
			}
			return err
		}
	}

	// Loops must have at least one tag.
	if len(tags) == 0 {
		return p.makeSyntaxError("loop_ header must have at least one tag")
	}

	// Handle the values.
	values := make([]any, 0)
	valuesRemain := true
	stopLooping := func() error {
		valuesRemain = false
		return nil
	}
	for valuesRemain {
		err := p.peekAndHandle(map[string]func() error{
			// Stop looping if we find a tag or reserved word.
			"save_": stopLooping,
			"data_": stopLooping,
			"loop_": stopLooping,
			"_":     stopLooping,
		}, func() error {
			value, err := p.readValue()
			if err != nil {
				return err
			}

			values = append(values, value)

			err = p.skipWhitespace()
			if err != nil {
				// If we get a syntax error when skipping whitespace, we've reached EOF.
				if _, ok := err.(CIFSyntaxError); ok {
					stopLooping()
				} else {
					return err
				}
			}
			return nil
		})
		if err != nil {
			if err, ok := err.(CIFSyntaxError); ok {
				return err.Wrap("could not read values in loop_")
			}
			return err
		}
	}

	if len(values)%len(tags) != 0 {
		return p.makeSyntaxError("number of values provided in loop_ must be a multiple of number of tags (tags: %v, values: %v)", len(tags), len(values))
	}

	// Make empty slices to hold values.
	dataItems := p.currDataItems()
	for _, tag := range tags {
		if _, exists := dataItems[tag]; exists {
			return p.makeSyntaxError("tag %q already has a value (%v)", tag, dataItems[tag])
		}
		dataItems[tag] = make([]any, len(values)/len(tags))
	}

	// Store the values.
	row := 0
	for i, value := range values {
		tag := tags[i%len(tags)]
		dataItems[tag].([]any)[row] = value

		// Increment the row counter at the end of a row.
		if (i+1)%len(tags) == 0 {
			row++
		}
	}

	return nil
}

// skipComment skips a comment.
//
// Sets p.lastByteWasEOL to true.
func (p *Parser) skipComment() error {
	_, err := p.readUntil([]byte{'\n'})
	return err
}

// skipWhitespace skips whitespace.
//
// Sets p.lastByteWasEOL to true if the last byte it skips
// is an EOL.
func (p *Parser) skipWhitespace() error {
	foundWhitespace := false
	for {
		b, err := p.reader.ReadByte()
		if err == io.EOF {
			break
		} else if err != nil {
			return err
		}

		if !isWhitespace(b) {
			err := p.reader.UnreadByte()
			if err != nil {
				return err
			}
			break
		}

		if b == '\n' {
			p.line++
		}
		foundWhitespace = true
	}

	if !foundWhitespace {
		return p.makeSyntaxError("no whitespace found when whitespace was expected")
	}
	return nil
}

/* ----------- Reader functions -------------- */

// readTextField reads in a text field.
//
// Can return an empty string. Returns a syntax error if
// an io.EOF is encountered.
//
// Assumes whitespace before the text field has been skipped,
// including the EOL required immediately before the semicolon.
func (p *Parser) readTextField() (string, error) {
	firstChar, err := p.reader.ReadByte()
	if err != nil {
		return "", err
	} else if firstChar != ';' {
		return "", p.makeSyntaxError("text field must begin with ';', not '%v'", firstChar)
	}

	res, err := p.readUntil([]byte{';'})
	if err != nil {
		if err, ok := err.(CIFSyntaxError); ok {
			return "", err.Wrap("could not read text field")
		}
		return "", err
	}

	return res, nil
}

// readValue reads a value and returns it.
//
// Appropriately handles inapplicable ('.'), unknown ('?'),
// unquoted string, quoted string, numeric, and text field values.
// Does not include string/text field delimeters in the
// return value.
func (p *Parser) readValue() (any, error) {
	var res any

	err := p.peekAndHandle(
		map[string]func() error{
			";": func() error {
				// Check if the last read byte was an EOL to disambiguate
				// between a text field and an unquoted string.
				var err error
				if p.lastByteWasEOL {
					res, err = p.readTextField()
					if err != nil {
						return err
					}
				}

				res, err = p.readUnquotedValue()
				return err
			},
			"'": func() error {
				var err error
				res, err = p.readQuotedString('\'')
				return err
			},
			"\"": func() error {
				var err error
				res, err = p.readQuotedString('"')
				return err
			},
		},
		func() error {
			var err error
			res, err = p.readUnquotedValue()
			return err
		},
	)
	if err != nil {
		if err, ok := err.(CIFSyntaxError); ok {
			return "", err.Wrap("could not read value")
		}
		return nil, err
	}

	return res, nil
}

// readUnquotedValue reads an unquoted string, numeric value, or
// the special '.' and '?' characters and returns it.
func (p *Parser) readUnquotedValue() (any, error) {
	res, err := p.readUntilWhitespace()
	if err != nil {
		if err, ok := err.(CIFSyntaxError); ok {
			return "", err.Wrap("could not read unquoted value")
		}
	}

	// Handle special values.
	if res == "." {
		return Inapplicable, nil
	} else if res == "?" {
		return Unknown, nil
	}

	if numeric, ok := parseNumeric(res); ok {
		return numeric, nil
	}

	return res, nil
}

// readQuotedString reads a quoted string and returns its value.
func (p *Parser) readQuotedString(quote byte) (string, error) {
	// Drop the opening quote.
	_, err := p.reader.Discard(1)
	if err != nil {
		return "", err
	}

	res, err := p.readUntil([]byte{quote})
	if err != nil {
		if err, ok := err.(CIFSyntaxError); ok {
			return "", err.Wrap("could not read %v quoted string", quote)
		}
		return "", err
	}

	// Skip the ending quote.
	_, err = p.reader.ReadByte()
	if err != nil {
		return "", err
	}

	return res, nil
}

// readTag reads a tag and returns it.
func (p *Parser) readTag() (string, error) {
	res, err := p.readUntilWhitespace()
	if err != nil {
		if err, ok := err.(CIFSyntaxError); ok {
			return "", err.Wrap("could not read tag")
		}
		return "", err
	}

	if !strings.HasPrefix(res, "_") || res == "_" {
		return "", p.makeSyntaxError("invalid tag name %q", res)
	}

	return res, nil
}

// readUntilWhitespace reads until whitespace or an io.EOF.
//
// Returns a syntax error if no characters could be read.
//
// Appropriately sets p.lastByteWasEOL.
func (p *Parser) readUntilWhitespace() (string, error) {
	return p.readUntil([]byte(whitespaceChars))
}

// readUntil reads until any character in delimSet or an io.EOF
// is encountered.
//
// Returns a syntax error if no characters could be read.
//
// Appropriately sets p.lastByteWasEOL.
func (p *Parser) readUntil(delimSet []byte) (string, error) {
	res := make([]byte, 0)
	delimFound := ""
	for {
		b, err := p.reader.ReadByte()
		if err == io.EOF {
			delimFound = "EOF"
			goto done
		} else if err != nil {
			return "", err
		}

		for _, delim := range delimSet {
			if b == delim {
				p.reader.UnreadByte()
				if err != nil {
					return "", err
				}
				delimFound = string([]byte{delim})
				goto done
			}
		}

		if b == '\n' {
			p.line++
			p.lastByteWasEOL = true
		} else {
			p.lastByteWasEOL = false
		}
		res = append(res, b)
	}

done:
	if len(res) == 0 {
		return "", p.makeSyntaxError("%q encountered before any characters could be read", delimFound)
	}
	return string(res), nil
}

/* ----------- Utility functions -------------- */

// currDataItems returns the current DataItems map.
//
// Must be called from within a data block or save frame,
// panics if not.
func (p *Parser) currDataItems() map[string]any {
	if p.saveFrameName != "" {
		return p.cif.DataBlocks[p.dataBlockName].SaveFrames[p.saveFrameName].DataItems
	}

	return p.cif.DataBlocks[p.dataBlockName].DataItems
}

// peekAndHandle peeks until a matching handler function is found in handlers.
// If no matching handler is found or an io.EOF is reached, calls defaultHandler.
// Does not consume any bytes from the reader.
//
// Matching is performed on the keys of the handlers map. Matches the shortest
// possible key, i.e. having "a" and "asdf" as keys in the map will result in
// the "a" handler shadowing the "asdf" handler.
//
// Because peekAndHandle peeks from the reader, handlers may not unread a byte
// before first performing a read operation.
func (p *Parser) peekAndHandle(handlers map[string]func() error, defaultHandler func() error) error {
	// Find what the longest handler key is.
	limit := 0
	for k := range handlers {
		if len(k) > limit {
			limit = len(k)
		}
	}

	// Find the matching handler function.
	n := 0
	var handlerFound bool
	for n <= limit || handlerFound {
		peek, err := p.reader.Peek(n)
		if err == io.EOF {
			defaultHandler()
		} else if err != nil {
			return err
		}

		handler, exists := handlers[strings.ToLower(string(peek))]
		if exists {
			return handler()
		}

		n++
	}

	return defaultHandler()
}

// makeSyntaxError returns an error message tagged with the
// line at which the syntax error was encountered.
func (p *Parser) makeSyntaxError(format string, a ...any) error {
	return CIFSyntaxError{
		Msg:  fmt.Sprintf(format, a...),
		Line: p.line,
	}
}

// isWhitespace returns whether or not a character is whitespace.
func isWhitespace(b byte) bool {
	for _, char := range whitespaceChars {
		if char == rune(b) {
			return true
		}
	}
	return false
}

// parseNumeric attempts to convert a string to a numeric value.
func parseNumeric(s string) (any, bool) {
	var err error
	var val any

	// First, try to parse an int.
	val, err = strconv.ParseInt(s, 10, 64)
	if err == nil {
		return val, true
	}

	// If that fails, try to parse a uint.
	val, err = strconv.ParseUint(s, 10, 64)
	if err == nil {
		return val, true
	}

	// If that fails, try to parse a float.
	val, err = strconv.ParseFloat(s, 64)
	if err == nil {
		return val, true
	}

	return nil, false
}
