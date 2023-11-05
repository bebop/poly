/*
Package genbank provides genbank parsers and writers.

GenBank is a flat text file format developed in the 1980s to annotate genetic
sequences, and has since become the standard for sharing annotated genetic
sequences.

This package provides a parser and writer to convert between the GenBank file
format and the more general Genbank struct.
*/
package genbank

import (
	"bufio"
	"bytes"
	"fmt"
	"io"
	"os"
	"regexp"
	"strconv"
	"strings"

	"github.com/TimothyStiles/poly/transform"
	"github.com/lunny/log"
	"github.com/mitchellh/go-wordwrap"
)

/******************************************************************************

GBK specific IO related things begin here.

******************************************************************************/

var (
	readFileFn        = os.ReadFile
	parseMultiNthFn   = parseMultiNth
	parseReferencesFn = parseReferences
)

// Genbank is the main struct for the Genbank file format.
type Genbank struct {
	Meta     Meta
	Features []Feature
	Sequence string // will be changed and include reader, writer, and byte slice.
}

// Meta holds the meta data for Genbank and other annotated sequence files.
type Meta struct {
	Date                 string            `json:"date"`
	Definition           string            `json:"definition"`
	Accession            string            `json:"accession"`
	Version              string            `json:"version"`
	Keywords             string            `json:"keywords"`
	Organism             string            `json:"organism"`
	Source               string            `json:"source"`
	Taxonomy             []string          `json:"taxonomy"`
	Origin               string            `json:"origin"`
	Locus                Locus             `json:"locus"`
	References           []Reference       `json:"references"`
	BaseCount            []BaseCount       `json:"base_count"`
	Other                map[string]string `json:"other"`
	Name                 string            `json:"name"`
	SequenceHash         string            `json:"sequence_hash"`
	SequenceHashFunction string            `json:"hash_function"`
}

// Feature holds the information for a feature in a Genbank file and other annotated sequence files.
type Feature struct {
	Type                 string              `json:"type"`
	Description          string              `json:"description"`
	Attributes           map[string][]string `json:"attributes"`
	SequenceHash         string              `json:"sequence_hash"`
	SequenceHashFunction string              `json:"hash_function"`
	Sequence             string              `json:"sequence"`
	Location             Location            `json:"location"`
	ParentSequence       *Genbank            `json:"-"`
}

// Reference holds information for one reference in a Meta struct.
type Reference struct {
	Authors    string `json:"authors"`
	Title      string `json:"title"`
	Journal    string `json:"journal"`
	PubMed     string `json:"pub_med"`
	Remark     string `json:"remark"`
	Range      string `json:"range"`
	Consortium string `json:"consortium"`
}

// Locus holds Locus information in a Meta struct.
type Locus struct {
	Name             string `json:"name"`
	SequenceLength   string `json:"sequence_length"`
	MoleculeType     string `json:"molecule_type"`
	GenbankDivision  string `json:"genbank_division"`
	ModificationDate string `json:"modification_date"`
	SequenceCoding   string `json:"sequence_coding"`
	Circular         bool   `json:"circular"`
}

// Location is a struct that holds the location of a feature.
type Location struct {
	Start             int        `json:"start"`
	End               int        `json:"end"`
	Complement        bool       `json:"complement"`
	Join              bool       `json:"join"`
	FivePrimePartial  bool       `json:"five_prime_partial"`
	ThreePrimePartial bool       `json:"three_prime_partial"`
	GbkLocationString string     `json:"gbk_location_string"`
	SubLocations      []Location `json:"sub_locations"`
}

// BaseCount is a struct that holds the base counts for a sequence.
type BaseCount struct {
	Base  string
	Count int
}

// Precompiled regular expressions:
var (
	basePairRegex         = regexp.MustCompile(` \d* \w{2} `)
	circularRegex         = regexp.MustCompile(` circular `)
	modificationDateRegex = regexp.MustCompile(`\d{2}-[A-Z]{3}-\d{4}`)
	partialRegex          = regexp.MustCompile("<|>")
	sequenceRegex         = regexp.MustCompile("[^a-zA-Z]+")
)

// StoreFeatureSequences calls StoreSequence on all features.
// The resulting JSON is guaranteed to have useful Feature.Sequence values.
// Useful when exporting for downstream analysis, such as with json.Marshal.
func (sequence *Genbank) StoreFeatureSequences() error {
	for i := range sequence.Features {
		_, err := sequence.Features[i].StoreSequence()
		if err != nil {
			return err
		}
	}
	return nil
}

// AddFeature adds a feature to a Genbank struct.
// NOTE: This method assumes feature is not referenced in another location
// as this only creates a shallow copy.
// If you intend to duplicate a feature from another Genbank and plan
// to modify in either location, it is recommended you first call feature.Copy()
// before passing as input to save yourself trouble.
func (sequence *Genbank) AddFeature(feature *Feature) error {
	feature.ParentSequence = sequence
	sequence.Features = append(sequence.Features, *feature)
	return nil
}

// GetSequence returns the sequence of a feature.
func (feature Feature) GetSequence() (string, error) {
	return getFeatureSequence(feature, feature.Location)
}

// StoreSequence infers and assigns the value of feature.Sequence
// if currently an empty string.
func (feature *Feature) StoreSequence() (string, error) {
	if feature.Sequence != "" {
		return feature.Sequence, nil
	}
	seq, err := getFeatureSequence(*feature, feature.Location)
	if err == nil {
		feature.Sequence = seq
	}
	return seq, err
}

// Copy creates deep copy of Feature, which supports safe duplication.
func (feature *Feature) Copy() Feature {
	copy := *feature
	copy.Location = CopyLocation(feature.Location)
	copy.Attributes = NewMultiMap[string, string]()
	ForEachKey(feature.Attributes, func(k string, v []string) {
		copy.Attributes[k] = MapSlice(v, identity[string])
	})
	return copy
}

// CopyLocation creates deep copy of Location, which supports safe duplication
func CopyLocation(location Location) Location {
	location.SubLocations = MapSlice(location.SubLocations, CopyLocation)
	return location
}

// getFeatureSequence takes a feature and location object and returns a sequence string.
func getFeatureSequence(feature Feature, location Location) (string, error) {
	var sequenceBuffer bytes.Buffer
	parentSequence := feature.ParentSequence.Sequence

	if len(location.SubLocations) == 0 {
		sequenceBuffer.WriteString(parentSequence[location.Start:location.End])
	} else {
		for _, subLocation := range location.SubLocations {
			sequence, err := getFeatureSequence(feature, subLocation)
			if err != nil {
				return "", err
			}

			sequenceBuffer.WriteString(sequence)
		}
	}

	// reverse complements resulting string if needed.
	sequenceString := sequenceBuffer.String()
	if location.Complement {
		sequenceString = transform.ReverseComplement(sequenceString)
	}

	return sequenceString, nil
}

// WriteTo implements the io.WriterTo interface on genbank records.
func (sequence *Genbank) WriteTo(w io.Writer) (int64, error) {
	var writtenBytes int64
	var newWrittenBytes int
	var err error

	locus := sequence.Meta.Locus
	var shape string

	if locus.Circular {
		shape = "circular"
	} else {
		shape = "linear"
	}

	fivespace := generateWhiteSpace(subMetaIndex)

	// building locus
	locusData := locus.Name + fivespace + locus.SequenceLength + " bp" + fivespace + locus.MoleculeType + fivespace + shape + fivespace + locus.GenbankDivision + fivespace + locus.ModificationDate
	locusString := "LOCUS       " + locusData + "\n"
	newWrittenBytes, err = w.Write([]byte(locusString))
	writtenBytes += int64(newWrittenBytes)
	if err != nil {
		return writtenBytes, err
	}

	// building other standard meta features
	definitionString := buildMetaString("DEFINITION", sequence.Meta.Definition)
	newWrittenBytes, err = w.Write([]byte(definitionString))
	writtenBytes += int64(newWrittenBytes)
	if err != nil {
		return writtenBytes, err
	}

	accessionString := buildMetaString("ACCESSION", sequence.Meta.Accession)
	newWrittenBytes, err = w.Write([]byte(accessionString))
	writtenBytes += int64(newWrittenBytes)
	if err != nil {
		return writtenBytes, err
	}

	versionString := buildMetaString("VERSION", sequence.Meta.Version)
	newWrittenBytes, err = w.Write([]byte(versionString))
	writtenBytes += int64(newWrittenBytes)
	if err != nil {
		return writtenBytes, err
	}

	keywordsString := buildMetaString("KEYWORDS", sequence.Meta.Keywords)
	newWrittenBytes, err = w.Write([]byte(keywordsString))
	writtenBytes += int64(newWrittenBytes)
	if err != nil {
		return writtenBytes, err
	}

	sourceString := buildMetaString("SOURCE", sequence.Meta.Source)
	newWrittenBytes, err = w.Write([]byte(sourceString))
	writtenBytes += int64(newWrittenBytes)
	if err != nil {
		return writtenBytes, err
	}

	organismString := buildMetaString("  ORGANISM", sequence.Meta.Organism)
	newWrittenBytes, err = w.Write([]byte(organismString))
	writtenBytes += int64(newWrittenBytes)
	if err != nil {
		return writtenBytes, err
	}

	if len(sequence.Meta.Taxonomy) > 0 {
		var taxonomyString strings.Builder
		for i, taxonomyData := range sequence.Meta.Taxonomy {
			taxonomyString.WriteString(taxonomyData)
			if len(sequence.Meta.Taxonomy) == i+1 {
				taxonomyString.WriteString(".")
			} else {
				taxonomyString.WriteString("; ")
			}
		}
		newWrittenBytes, err = w.Write([]byte(buildMetaString("", taxonomyString.String())))
		writtenBytes += int64(newWrittenBytes)
		if err != nil {
			return writtenBytes, err
		}
	}

	// building references
	// TODO: could use reflection to get keys and make more general.
	for referenceIndex, reference := range sequence.Meta.References {
		referenceString := buildMetaString("REFERENCE", fmt.Sprintf("%d  %s", referenceIndex+1, reference.Range))
		newWrittenBytes, err = w.Write([]byte(referenceString))
		writtenBytes += int64(newWrittenBytes)
		if err != nil {
			return writtenBytes, err
		}

		if reference.Authors != "" {
			authorsString := buildMetaString("  AUTHORS", reference.Authors)
			newWrittenBytes, err = w.Write([]byte(authorsString))
			writtenBytes += int64(newWrittenBytes)
			if err != nil {
				return writtenBytes, err
			}
		}

		if reference.Title != "" {
			titleString := buildMetaString("  TITLE", reference.Title)
			newWrittenBytes, err = w.Write([]byte(titleString))
			writtenBytes += int64(newWrittenBytes)
			if err != nil {
				return writtenBytes, err
			}
		}

		if reference.Journal != "" {
			journalString := buildMetaString("  JOURNAL", reference.Journal)
			newWrittenBytes, err = w.Write([]byte(journalString))
			writtenBytes += int64(newWrittenBytes)
			if err != nil {
				return writtenBytes, err
			}
		}

		if reference.PubMed != "" {
			pubMedString := buildMetaString("  PUBMED", reference.PubMed)
			newWrittenBytes, err = w.Write([]byte(pubMedString))
			writtenBytes += int64(newWrittenBytes)
			if err != nil {
				return writtenBytes, err
			}
		}
		if reference.Consortium != "" {
			consrtmString := buildMetaString("  CONSRTM", reference.Consortium)
			newWrittenBytes, err = w.Write([]byte(consrtmString))
			writtenBytes += int64(newWrittenBytes)
			if err != nil {
				return writtenBytes, err
			}
		}
	}

	// building other meta fields that are catch all
	otherKeys := make([]string, 0, len(sequence.Meta.Other))
	for key := range sequence.Meta.Other {
		otherKeys = append(otherKeys, key)
	}

	for _, otherKey := range otherKeys {
		otherString := buildMetaString(otherKey, sequence.Meta.Other[otherKey])
		newWrittenBytes, err = w.Write([]byte(otherString))
		writtenBytes += int64(newWrittenBytes)
		if err != nil {
			return writtenBytes, err
		}
	}

	// start writing features section.
	newWrittenBytes, err = w.Write([]byte("FEATURES             Location/Qualifiers\n"))
	writtenBytes += int64(newWrittenBytes)
	if err != nil {
		return writtenBytes, err
	}
	for _, feature := range sequence.Features {
		newWrittenBytes, err = w.Write([]byte(BuildFeatureString(feature)))
		writtenBytes += int64(newWrittenBytes)
		if err != nil {
			return writtenBytes, err
		}
	}

	// start writing base count
	if len(sequence.Meta.BaseCount) > 0 {
		newWrittenBytes, err = w.Write([]byte("BASE COUNT    "))
		writtenBytes += int64(newWrittenBytes)
		if err != nil {
			return writtenBytes, err
		}
		for _, baseCount := range sequence.Meta.BaseCount {
			newWrittenBytes, err = w.Write([]byte(strconv.Itoa(baseCount.Count) + " " + baseCount.Base + "   "))
			writtenBytes += int64(newWrittenBytes)
			if err != nil {
				return writtenBytes, err
			}
		}
		newWrittenBytes, err = w.Write([]byte("\n"))
		writtenBytes += int64(newWrittenBytes)
		if err != nil {
			return writtenBytes, err
		}
	}

	// start writing sequence section.
	newWrittenBytes, err = w.Write([]byte("ORIGIN\n"))
	writtenBytes += int64(newWrittenBytes)
	if err != nil {
		return writtenBytes, err
	}

	// iterate over every character in sequence range.
	for index, base := range sequence.Sequence {
		// if 60th character add newline then whitespace and index number and space before adding next base.
		if index%60 == 0 {
			if index != 0 {
				newWrittenBytes, err = w.Write([]byte("\n"))
				writtenBytes += int64(newWrittenBytes)
				if err != nil {
					return writtenBytes, err
				}
			}
			lineNumberString := strconv.Itoa(index + 1)          // genbank indexes at 1 for some reason
			leadingWhiteSpaceLength := 9 - len(lineNumberString) // <- I wish I was kidding
			for i := 0; i < leadingWhiteSpaceLength; i++ {
				newWrittenBytes, err = w.Write([]byte(" "))
				writtenBytes += int64(newWrittenBytes)
				if err != nil {
					return writtenBytes, err
				}
			}
			newWrittenBytes, err = w.Write([]byte(lineNumberString + " "))
			writtenBytes += int64(newWrittenBytes)
			if err != nil {
				return writtenBytes, err
			}
			newWrittenBytes, err = w.Write([]byte{byte(base)})
			writtenBytes += int64(newWrittenBytes)
			if err != nil {
				return writtenBytes, err
			}
			// if base index is divisible by ten add a space (genbank convention)
		} else if index%10 == 0 {
			newWrittenBytes, err = w.Write([]byte(" "))
			writtenBytes += int64(newWrittenBytes)
			if err != nil {
				return writtenBytes, err
			}
			newWrittenBytes, err = w.Write([]byte{byte(base)})
			writtenBytes += int64(newWrittenBytes)
			if err != nil {
				return writtenBytes, err
			}
			// else just add the base.
		} else {
			newWrittenBytes, err = w.Write([]byte{byte(base)})
			writtenBytes += int64(newWrittenBytes)
			if err != nil {
				return writtenBytes, err
			}
		}
	}
	// finish genbank file with "//" on newline (again a genbank convention)
	newWrittenBytes, err = w.Write([]byte("\n//\n"))
	writtenBytes += int64(newWrittenBytes)
	if err != nil {
		return writtenBytes, err
	}

	return writtenBytes, nil
}

// ParseError represents failures encountered while parsing,
// and pointers to it are fully compatable with the error interface.
type ParseError struct {
	file   string // the file origin
	line   string // the offending line
	before bool   // whether the error occurred before or on this line
	lineNo int    // the line number, 0 indexed
	info   string `default:"syntax error"` // description of the error type
	wraps  error  // stores the error that led to this, if any
}

func (e ParseError) Error() string {
	var out, loc string
	if e.wraps == io.EOF {
		out = "unexpected EOF"
		if e.file != "" {
			return fmt.Sprintf("%s in %s", out, e.file)
		} else {
			return out
		}
	}
	if e.file == "" {
		loc = fmt.Sprintf("line %d", e.lineNo)
	} else {
		loc = fmt.Sprintf("%s:%d", e.file, e.lineNo)
	}
	if e.before {
		out = fmt.Sprintf("%s encountered before %s", e.info, loc)
	} else {
		out = fmt.Sprintf("%s encountered on %s: %s", e.info, loc, e.line)
	}
	if e.wraps != nil {
		out = fmt.Sprintf("%s\nfrom %v", out, e.wraps)
	}
	return out
}

// defines state for the parser, and utility methods to modify
type parseLoopParameters struct {
	newLocation      bool
	attribute        string
	attributeValue   string
	emptyAttribute   bool
	sequenceBuilder  strings.Builder
	parseStep        string
	genbank          Genbank // since we are scanning lines we need a Genbank struct to store the data outside the loop.
	feature          Feature
	features         []Feature
	metadataTag      string
	metadataData     []string //this stutters but will remain to make it easier to batch rename variables when compared to parameters.metadataTag.
	genbankStarted   bool
	currentLine      string
	prevline         string
	multiLineFeature bool
}

// method to init loop parameters
func (params *parseLoopParameters) init() {
	params.newLocation = true
	params.feature.Attributes = NewMultiMap[string, string]()
	params.parseStep = "metadata"
	params.genbankStarted = false
	params.genbank.Meta.Other = make(map[string]string)
}

// save our completed attribute / qualifier string to the current feature
// useful as a wrap-up step from multiple states
func (params *parseLoopParameters) saveLastAttribute() {
	newValue := params.attributeValue != ""
	emptyType := params.feature.Type != ""
	if newValue || emptyType {
		if newValue {
			Put(params.feature.Attributes, params.attribute, params.attributeValue)
		}
		params.features = append(params.features, params.feature)

		// reset attribute state
		params.attributeValue = ""
		params.attribute = ""
		params.feature = Feature{}
		params.feature.Attributes = NewMultiMap[string, string]()
	}
}

// Header is a blank struct, needed for compatibility with bio parsers. It contains nothing.
type Header struct{}

// WriteTo is a blank function, needed for compatibility with bio parsers. It doesn't do anything.
func (header *Header) WriteTo(w io.Writer) (int64, error) {
	return 0, nil
}

// Parser is a genbank parser created on an io.Reader.
type Parser struct {
	scanner    bufio.Scanner
	parameters parseLoopParameters
}

// Header returns nil,nil.
func (parser *Parser) Header() (*Header, error) {
	return &Header{}, nil
}

// NewParser returns a Parser that uses r as the source
// from which to parse genbank formatted sequences.
func NewParser(r io.Reader, maxLineSize int) *Parser {
	scanner := bufio.NewScanner(r)
	buf := make([]byte, maxLineSize)
	scanner.Buffer(buf, maxLineSize)
	return &Parser{
		scanner: *scanner,
	}
}

// Next takes in a reader representing a multi gbk/gb/genbank file and outputs the next record
func (parser *Parser) Next() (*Genbank, error) {
	parser.parameters.init()
	// Loop through each line of the file
	for lineNum := 0; parser.scanner.Scan(); lineNum++ {
		// get line from scanner and split it
		line := parser.scanner.Text()
		splitLine := strings.Split(strings.TrimSpace(line), " ")

		prevline := parser.parameters.currentLine
		parser.parameters.currentLine = line
		parser.parameters.prevline = prevline

		// keep scanning until we find the start of the first record
		if !parser.parameters.genbankStarted {
			// We detect the beginning of a new genbank file with "LOCUS"
			locusFlag := strings.Contains(line, "LOCUS")

			if locusFlag {
				parser.parameters = parseLoopParameters{}
				parser.parameters.init()
				parser.parameters.genbank.Meta.Locus = parseLocus(line)
				parser.parameters.genbankStarted = true
			}
			continue
		}

		// define parser state machine
		switch parser.parameters.parseStep {
		case "metadata":
			// Handle empty lines
			if len(line) == 0 {
				return &Genbank{}, &ParseError{line: line, lineNo: lineNum, info: "unexpected empty metadata"}
			}

			// If we are currently reading a line, we need to figure out if it is a new meta line.
			if string(line[0]) != " " || parser.parameters.metadataTag == "FEATURES" {
				// If this is true, it means we are beginning a new meta tag. In that case, let's save
				// the older data, and then continue along.
				switch parser.parameters.metadataTag {
				case "DEFINITION":
					parser.parameters.genbank.Meta.Definition = parseMetadata(parser.parameters.metadataData)
				case "ACCESSION":
					parser.parameters.genbank.Meta.Accession = parseMetadata(parser.parameters.metadataData)
				case "VERSION":
					parser.parameters.genbank.Meta.Version = parseMetadata(parser.parameters.metadataData)
				case "KEYWORDS":
					parser.parameters.genbank.Meta.Keywords = parseMetadata(parser.parameters.metadataData)
				case "SOURCE":
					parser.parameters.genbank.Meta.Source, parser.parameters.genbank.Meta.Organism, parser.parameters.genbank.Meta.Taxonomy = getSourceOrganism(parser.parameters.metadataData)
				case "REFERENCE":
					// parseReferencesFn = parseReferences in genbank_test. We use Fn for testing purposes.
					reference, err := parseReferencesFn(parser.parameters.metadataData)
					if err != nil {
						return &Genbank{}, &ParseError{line: line, lineNo: lineNum, before: true, wraps: err, info: "failed in parsing reference"}
					}
					parser.parameters.genbank.Meta.References = append(parser.parameters.genbank.Meta.References, reference)

				case "FEATURES":
					parser.parameters.parseStep = "features"

					// We know that we are now parsing features, so lets initialize our first feature
					parser.parameters.feature.Type = strings.TrimSpace(splitLine[0])
					parser.parameters.feature.Location.GbkLocationString = strings.TrimSpace(splitLine[len(splitLine)-1])
					parser.parameters.newLocation = true

					continue

				default:
					if parser.parameters.metadataTag != "" {
						parser.parameters.genbank.Meta.Other[parser.parameters.metadataTag] = parseMetadata(parser.parameters.metadataData)
					}
				}

				parser.parameters.metadataTag = strings.TrimSpace(splitLine[0])
				parser.parameters.metadataData = []string{strings.TrimSpace(line[len(parser.parameters.metadataTag):])}
			} else {
				parser.parameters.metadataData = append(parser.parameters.metadataData, line)
			}
		case "features":
			baseCountFlag := strings.Contains(line, "BASE COUNT") // example string for BASE COUNT: "BASE COUNT    67070277 a   48055043 c   48111528 g   67244164 t   18475410 n"
			if baseCountFlag {
				fields := strings.Fields(line)
				for countIndex := 2; countIndex < len(fields)-1; countIndex += 2 { // starts at two because we don't want to include "BASE COUNT" in our fields
					count, err := strconv.Atoi(fields[countIndex])
					if err != nil {
						return &Genbank{}, &ParseError{line: line, lineNo: lineNum, wraps: err, info: "invalid base count"}
					}

					baseCount := BaseCount{
						Base:  fields[countIndex+1],
						Count: count,
					}
					parser.parameters.genbank.Meta.BaseCount = append(parser.parameters.genbank.Meta.BaseCount, baseCount)
				}
				break
			}
			// Switch to sequence parsing
			originFlag := strings.Contains(line, "ORIGIN") // we detect the beginning of the sequence with "ORIGIN"
			contigFlag := strings.Contains(line, "CONTIG")
			if originFlag || contigFlag {
				parser.parameters.parseStep = "sequence"

				parser.parameters.saveLastAttribute()

				// add our features to the genbank
				for _, feature := range parser.parameters.features {
					// TODO: parse location when line is read, or track line number so error is localized
					location, err := parseLocation(feature.Location.GbkLocationString)
					if err != nil {
						return &Genbank{}, &ParseError{before: true, line: line, lineNo: lineNum, wraps: err, info: "invalid feature location"}
					}
					feature.Location = location
					err = parser.parameters.genbank.AddFeature(&feature)
					if err != nil {
						return &Genbank{}, &ParseError{before: true, line: line, lineNo: lineNum, wraps: err, info: "problem adding feature"}
					}
				}

				if contigFlag {
					parser.parameters.genbank.Meta.Other["CONTIG"] = parseMetadata(splitLine[1:])
				}
				continue
			}

			// check if current line contains anything but whitespace
			trimmedLine := strings.TrimSpace(line)
			if len(trimmedLine) == 0 {
				continue
			}

			indent := countLeadingSpaces(parser.parameters.currentLine)
			// determine if current line is a new top level feature
			if indent == 0 {
				return &Genbank{}, &ParseError{line: line, lineNo: lineNum, info: "unexpected metadata when parsing feature"}
			} else if indent < countLeadingSpaces(parser.parameters.prevline) || parser.parameters.prevline == "FEATURES" {
				parser.parameters.saveLastAttribute()

				parser.parameters.feature = Feature{}
				parser.parameters.feature.Attributes = NewMultiMap[string, string]()

				// An initial feature line looks like this: `source          1..2686` with a type separated by its location
				if len(splitLine) < 2 {
					return &Genbank{}, &ParseError{line: line, lineNo: lineNum, info: "malformed feature"}
				}
				parser.parameters.feature.Type = strings.TrimSpace(splitLine[0])
				parser.parameters.feature.Location.GbkLocationString = strings.TrimSpace(splitLine[len(splitLine)-1])
				parser.parameters.multiLineFeature = false // without this we can't tell if something is a multiline feature or multiline qualifier
			} else if !strings.Contains(parser.parameters.currentLine, "/") { // current line is continuation of a feature or qualifier (sub-constituent of a feature)
				// if it's a continuation of the current feature, add it to the location
				if !strings.Contains(parser.parameters.currentLine, "\"") && (countLeadingSpaces(parser.parameters.currentLine) > countLeadingSpaces(parser.parameters.prevline) || parser.parameters.multiLineFeature) {
					parser.parameters.feature.Location.GbkLocationString += strings.TrimSpace(line)
					parser.parameters.multiLineFeature = true // without this we can't tell if something is a multiline feature or multiline qualifier
				} else { // it's a continued line of a qualifier
					removeAttributeValueQuotes := strings.Replace(trimmedLine, "\"", "", -1)

					parser.parameters.attributeValue = parser.parameters.attributeValue + removeAttributeValueQuotes
				}
			} else if strings.Contains(parser.parameters.currentLine, "/") { // current line is a new qualifier
				trimmedCurrentLine := strings.TrimSpace(parser.parameters.currentLine)
				if trimmedCurrentLine[0] != '/' { // if we have an exception case, like (adenine(1518)-N(6)/adenine(1519)-N(6))-
					parser.parameters.attributeValue = parser.parameters.attributeValue + trimmedCurrentLine
					continue
				}
				// save our completed attribute / qualifier string to the current feature
				if parser.parameters.attributeValue != "" || parser.parameters.emptyAttribute {
					Put(parser.parameters.feature.Attributes, parser.parameters.attribute, parser.parameters.attributeValue)
					parser.parameters.emptyAttribute = false
				}
				parser.parameters.attributeValue = ""
				splitAttribute := strings.Split(line, "=")
				trimmedSpaceAttribute := strings.TrimSpace(splitAttribute[0])
				removedForwardSlashAttribute := strings.Replace(trimmedSpaceAttribute, "/", "", 1)

				parser.parameters.attribute = removedForwardSlashAttribute

				var removeAttributeValueQuotes string
				if len(splitAttribute) == 1 { // handle case of ` /pseudo `, which has no text
					removeAttributeValueQuotes = ""
					parser.parameters.emptyAttribute = true
				} else { // this is normally triggered
					removeAttributeValueQuotes = strings.Replace(splitAttribute[1], "\"", "", -1)
				}
				parser.parameters.attributeValue = removeAttributeValueQuotes
				parser.parameters.multiLineFeature = false // without this we can't tell if something is a multiline feature or multiline qualifier
			} else {
				return &Genbank{}, &ParseError{line: line, lineNo: lineNum, info: "invalid feature"}
			}

		case "sequence":
			if len(line) < 2 { // throw error if line is malformed
				return &Genbank{}, &ParseError{line: line, lineNo: lineNum, info: "too short line found while parsing genbank sequence"}
			} else if line[0:2] == "//" { // end of sequence
				parser.parameters.genbank.Sequence = parser.parameters.sequenceBuilder.String()

				parser.parameters.genbankStarted = false
				parser.parameters.sequenceBuilder.Reset()
				return &parser.parameters.genbank, nil
			} else { // add line to total sequence
				parser.parameters.sequenceBuilder.WriteString(sequenceRegex.ReplaceAllString(line, ""))
			}
		default:
			log.Warnf("Unknown parse step: %s", parser.parameters.parseStep)
			parser.parameters.genbankStarted = false
		}
	}
	return &Genbank{}, io.EOF
}

func countLeadingSpaces(line string) int {
	return len(line) - len(strings.TrimLeft(line, " "))
}

func parseMetadata(metadataData []string) string {
	var outputMetadata string
	if len(metadataData) == 0 {
		return "."
	}
	for _, data := range metadataData {
		outputMetadata = outputMetadata + strings.TrimSpace(data) + " "
	}
	outputMetadata = outputMetadata[:len(outputMetadata)-1] // Remove trailing whitespace
	return outputMetadata
}

func parseReferences(metadataData []string) (Reference, error) {
	var reference Reference
	var err error
	rangeIndex := strings.Index(metadataData[0], "(")
	if rangeIndex != -1 {
		reference.Range = metadataData[0][rangeIndex:]
	}
	var referenceKey string
	var referenceValue string

	if len(metadataData) == 1 {
		return Reference{}, fmt.Errorf("Got reference with no additional information")
	}

	referenceKey = strings.Split(strings.TrimSpace(metadataData[1]), " ")[0]
	referenceValue = strings.TrimSpace(metadataData[1][len(referenceKey)+2:])
	for index := 2; index < len(metadataData); index++ {
		if len(metadataData[index]) > 3 {
			if metadataData[index][3] != ' ' {
				err = reference.addKey(referenceKey, referenceValue)
				if err != nil {
					return reference, err
				}
				referenceKey = strings.Split(strings.TrimSpace(metadataData[index]), " ")[0]
				referenceValue = strings.TrimSpace(metadataData[index][len(referenceKey)+2:])
			} else {
				// Otherwise, simply append the next metadata.
				referenceValue = referenceValue + " " + strings.TrimSpace(metadataData[index])
			}
		}
	}
	err = reference.addKey(referenceKey, referenceValue)
	if err != nil {
		return reference, err
	}

	return reference, nil
}

func (reference *Reference) addKey(referenceKey string, referenceValue string) error {
	switch referenceKey {
	case "AUTHORS":
		reference.Authors = referenceValue
	case "TITLE":
		reference.Title = referenceValue
	case "JOURNAL":
		reference.Journal = referenceValue
	case "PUBMED":
		reference.PubMed = referenceValue
	case "REMARK":
		reference.Remark = referenceValue
	case "CONSRTM":
		reference.Consortium = referenceValue
	default:
		return fmt.Errorf("ReferenceKey not in [AUTHORS, TITLE, JOURNAL, PUBMED, REMARK, CONSRTM]. Got: %s", referenceKey)
	}
	return nil
}

var genBankMoleculeTypes = []string{
	"DNA",
	"genomic DNA",
	"genomic RNA",
	"mRNA",
	"tRNA",
	"rRNA",
	"other RNA",
	"other DNA",
	"transcribed RNA",
	"viral cRNA",
	"unassigned DNA",
	"unassigned RNA",
}

// used in parseLocus function though it could be useful elsewhere.
var genbankDivisions = []string{
	"PRI", //primate sequences
	"ROD", //rodent sequences
	"MAM", //other mamallian sequences
	"VRT", //other vertebrate sequences
	"INV", //invertebrate sequences
	"PLN", //plant, fungal, and algal sequences
	"BCT", //bacterial sequences
	"VRL", //viral sequences
	"PHG", //bacteriophage sequences
	"SYN", //synthetic sequences
	"UNA", //unannotated sequences
	"EST", //EST sequences (expressed sequence tags)
	"PAT", //patent sequences
	"STS", //STS sequences (sequence tagged sites)
	"GSS", //GSS sequences (genome survey sequences)
	"HTG", //HTG sequences (high-throughput genomic sequences)
	"HTC", //unfinished high-throughput cDNA sequencing
	"ENV", //environmental sampling sequences
}

// TODO rewrite with proper error handling.
// parses locus from provided string.
func parseLocus(locusString string) Locus {
	locus := Locus{}

	locusSplit := strings.Split(strings.TrimSpace(locusString), " ")

	var filteredLocusSplit []string
	for i := range locusSplit {
		if locusSplit[i] != "" {
			filteredLocusSplit = append(filteredLocusSplit, locusSplit[i])
		}
	}

	locus.Name = filteredLocusSplit[1]

	// sequence length and coding
	baseSequenceLength := basePairRegex.FindString(locusString)
	if baseSequenceLength != "" {
		splitBaseSequenceLength := strings.Split(strings.TrimSpace(baseSequenceLength), " ")
		if len(splitBaseSequenceLength) == 2 {
			locus.SequenceLength = splitBaseSequenceLength[0]
			locus.SequenceCoding = splitBaseSequenceLength[1]
		}
	}

	// molecule type
	for _, moleculeType := range genBankMoleculeTypes {
		moleculeRegex, _ := regexp.Compile(moleculeType)
		match := string(moleculeRegex.Find([]byte(locusString)))
		if match != "" {
			locus.MoleculeType = match
			break
		}
	}

	// circularity flag
	if circularRegex.Match([]byte(locusString)) {
		locus.Circular = true
	}

	// genbank division
	for _, genbankDivision := range genbankDivisions {
		genbankDivisionRegex, _ := regexp.Compile(genbankDivision)
		match := string(genbankDivisionRegex.Find([]byte(locusString)))
		if match != "" {
			locus.GenbankDivision = match
			break
		}
	}

	// ModificationDate
	locus.ModificationDate = modificationDateRegex.FindString(locusString)

	return locus
}

// indices for random points of interests on a gbk line.
const subMetaIndex = 5
const qualifierIndex = 21

func getSourceOrganism(metadataData []string) (string, string, []string) {
	source := strings.TrimSpace(metadataData[0])
	var organism string
	var taxonomy []string
	for iterator := 1; iterator < len(metadataData); iterator++ {
		dataLine := metadataData[iterator]
		headString := strings.Split(strings.TrimSpace(dataLine), " ")[0]
		if headString == "ORGANISM" {
			index := strings.Index(dataLine, `ORGANISM`)
			organism = strings.TrimSpace(dataLine[index+len("ORGANISM"):])
			continue
		}
		for _, taxonomyData := range strings.Split(strings.TrimSpace(dataLine), ";") {
			taxonomyDataTrimmed := strings.TrimSpace(taxonomyData)
			// Taxonomy ends with a ".", which we check for here
			if len(taxonomyDataTrimmed) > 1 {
				if taxonomyDataTrimmed[len(taxonomyDataTrimmed)-1] == '.' {
					taxonomyDataTrimmed = taxonomyDataTrimmed[:len(taxonomyDataTrimmed)-1]
				}
				taxonomy = append(taxonomy, taxonomyDataTrimmed)
			}
		}
	}
	return source, organism, taxonomy
}

func parseLocation(locationString string) (Location, error) {
	var location Location
	location.GbkLocationString = locationString
	if !strings.ContainsAny(locationString, "(") { // Case checks for simple expression of x..x
		if !strings.ContainsAny(locationString, ".") { //Case checks for simple expression x
			position, err := strconv.Atoi(partialRegex.ReplaceAllString(locationString, ""))
			if err != nil {
				return Location{}, err
			}
			location = Location{Start: position, End: position}
		} else {
			// to remove FivePrimePartial and ThreePrimePartial indicators from start and end before converting to int.
			startEndSplit := strings.Split(locationString, "..")
			start, err := strconv.Atoi(partialRegex.ReplaceAllString(startEndSplit[0], ""))
			if err != nil {
				return Location{}, err
			}
			end, err := strconv.Atoi(partialRegex.ReplaceAllString(startEndSplit[1], ""))
			if err != nil {
				return Location{}, err
			}
			location = Location{Start: start - 1, End: end}
		}
	} else {
		firstOuterParentheses := strings.Index(locationString, "(")
		expression := locationString[firstOuterParentheses+1 : strings.LastIndex(locationString, ")")]
		switch command := locationString[0:firstOuterParentheses]; command {
		case "join":
			location.Join = true
			// This case checks for join(complement(x..x),complement(x..x)), or any more complicated derivatives
			if strings.ContainsAny(expression, "(") {
				firstInnerParentheses := strings.Index(expression, "(")
				ParenthesesCount := 1
				prevSubLocationStart := 0
				for i := firstInnerParentheses + 1; i < len(expression); i++ { // "(" is at 0, so we start at 1
					switch expression[i] {
					case '(':
						ParenthesesCount++
					case ')':
						ParenthesesCount--
					case ',':
						if ParenthesesCount == 0 {
							parsedSubLocation, err := parseLocation(expression[prevSubLocationStart:i])
							if err != nil {
								return Location{}, err
							}
							parsedSubLocation.GbkLocationString = locationString
							location.SubLocations = append(location.SubLocations, parsedSubLocation)
							prevSubLocationStart = i + 1
						}
					}
				}
				if ParenthesesCount != 0 {
					return Location{}, fmt.Errorf("Unbalanced parentheses")
				}
				parsedSubLocation, err := parseLocation(expression[prevSubLocationStart:])
				if err != nil {
					return Location{}, err
				}
				parsedSubLocation.GbkLocationString = locationString
				location.SubLocations = append(location.SubLocations, parsedSubLocation)
			} else { // This is the default join(x..x,x..x)
				for _, numberRange := range strings.Split(expression, ",") {
					joinLocation, err := parseLocation(numberRange)
					if err != nil {
						return Location{}, err
					}
					location.SubLocations = append(location.SubLocations, joinLocation)
				}
			}

		case "complement":
			// location.Complement = true
			subLocation, err := parseLocation(expression)
			if err != nil {
				return Location{}, err
			}
			subLocation.Complement = true
			subLocation.GbkLocationString = locationString
			location.SubLocations = append(location.SubLocations, subLocation)
		}
	}

	if strings.Contains(locationString, "<") {
		location.FivePrimePartial = true
	}

	if strings.Contains(locationString, ">") {
		location.ThreePrimePartial = true
	}

	// if excess root node then trim node. Maybe should just be handled with second arg?
	if location.Start == 0 && location.End == 0 && !location.Join && !location.Complement {
		location = location.SubLocations[0]
	}

	return location, nil
}

// buildMetaString is a helper function to build the meta section of genbank files.
func buildMetaString(name string, data string) string {
	keyWhitespaceTrailLength := 12 - len(name) // I wish I was kidding.
	var keyWhitespaceTrail string
	for i := 0; i < keyWhitespaceTrailLength; i++ {
		keyWhitespaceTrail += " "
	}
	name += keyWhitespaceTrail
	wrappedData := wordwrap.WrapString(data, 68)
	splitData := strings.Split(wrappedData, "\n")
	var returnData string
	for index, datum := range splitData {
		if index == 0 {
			returnData = name + datum + "\n"
		} else {
			returnData += generateWhiteSpace(12) + datum + "\n"
		}
	}

	return returnData
}

// BuildLocationString is a recursive function that takes a location object and creates a gbk location string for Build()
func BuildLocationString(location Location) string {
	var locationString string

	if location.Complement {
		location.Complement = false
		locationString = "complement(" + BuildLocationString(location) + ")"
	} else if location.Join {
		locationString = "join("
		for _, sublocation := range location.SubLocations {
			locationString += BuildLocationString(sublocation) + ","
		}
		locationString = strings.TrimSuffix(locationString, ",") + ")"
	} else {
		locationString = strconv.Itoa(location.Start+1) + ".." + strconv.Itoa(location.End)
		if location.FivePrimePartial {
			locationString = "<" + locationString
		}

		if location.ThreePrimePartial {
			locationString += ">"
		}
	}
	return locationString
}

// BuildFeatureString is a helper function to build gbk feature strings for Build()
func BuildFeatureString(feature Feature) string {
	whiteSpaceTrailLength := 16 - len(feature.Type) // I wish I was kidding.
	whiteSpaceTrail := generateWhiteSpace(whiteSpaceTrailLength)
	var location string

	if feature.Location.GbkLocationString != "" {
		location = feature.Location.GbkLocationString
	} else {
		location = BuildLocationString(feature.Location)
	}
	featureHeader := generateWhiteSpace(subMetaIndex) + feature.Type + whiteSpaceTrail + location + "\n"
	returnString := featureHeader

	if feature.Attributes != nil {
		ForEachValue(feature.Attributes, func(key string, value string) {
			returnString += generateWhiteSpace(qualifierIndex) + "/" + key + "=\"" + value + "\"\n"
		})
	}
	return returnString
}

func generateWhiteSpace(length int) string {
	var spaceBuilder strings.Builder

	for i := 0; i < length; i++ {
		spaceBuilder.WriteString(" ")
	}

	return spaceBuilder.String()
}

/******************************************************************************

GBK specific IO related things end here.

******************************************************************************/

/******************************************************************************
Old functions for testing here.

We used to have integrated read files, before we had the generic parser
interface. These are still here because of the need to switch over our tests.
******************************************************************************/

// read reads a GBK file from path and returns a Genbank struct.
func read(path string) (Genbank, error) {
	genbankSlice, err := readMultiNth(path, 1)
	if err != nil {
		return Genbank{}, err
	}
	genbank := genbankSlice[0]
	return genbank, err
}

// readMulti reads a multi Gbk from path and parses it into a slice of Genbank structs.
func readMulti(path string) ([]Genbank, error) {
	return readMultiNth(path, -1)
}

// readMultiNth reads a multi Gbk from path and parses N entries into a slice of Genbank structs.
func readMultiNth(path string, count int) ([]Genbank, error) {
	file, err := os.Open(path)
	if err != nil {
		return []Genbank{}, err
	}

	sequence, perr := parseMultiNthFn(file, count)
	if perr != nil {
		perr.file = path
		return []Genbank{}, perr
	}

	return sequence, nil
}

func parseMultiNth(r io.Reader, count int) ([]Genbank, *ParseError) {
	parser := NewParser(r, bufio.MaxScanTokenSize)
	var genbanks []Genbank
	for i := 0; i < count; i++ {
		gb, err := parser.Next()
		if err != nil {
			var perr *ParseError
			if err == io.EOF {
				perr = &ParseError{wraps: io.EOF}
			} else if err != nil {
				perr = err.(*ParseError)
			}
			return genbanks, perr
		}
		genbanks = append(genbanks, *gb)
	}
	return genbanks, nil
}

func parse(r io.Reader) (Genbank, error) {
	parser := NewParser(r, bufio.MaxScanTokenSize)
	gb, err := parser.Next()
	return *gb, err
}

func write(gb Genbank, path string) error {
	file, err := os.OpenFile(path, os.O_CREATE|os.O_WRONLY, 0644)
	if err != nil {
		return err
	}
	defer file.Close()
	_, err = gb.WriteTo(file)
	return err
}

func writeMulti(gbs []Genbank, path string) error {
	file, err := os.OpenFile(path, os.O_CREATE|os.O_WRONLY, 0644)
	if err != nil {
		return err
	}
	defer file.Close()
	for _, gb := range gbs {
		_, err = gb.WriteTo(file)
		if err != nil {
			return err
		}
	}
	return nil
}
