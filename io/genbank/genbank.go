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
	parseMultiNthFn   = ParseMultiNth
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
	Type                 string            `json:"type"`
	Description          string            `json:"description"`
	Attributes           map[string]string `json:"attributes"`
	SequenceHash         string            `json:"sequence_hash"`
	SequenceHashFunction string            `json:"hash_function"`
	Sequence             string            `json:"sequence"`
	Location             Location          `json:"location"`
	ParentSequence       *Genbank          `json:"-"`
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

// AddFeature adds a feature to a Genbank struct.
func (sequence *Genbank) AddFeature(feature *Feature) error {
	feature.ParentSequence = sequence
	sequence.Features = append(sequence.Features, *feature)
	return nil
}

// GetSequence returns the sequence of a feature.
func (feature Feature) GetSequence() (string, error) {
	return getFeatureSequence(feature, feature.Location)
}

// getFeatureSequence takes a feature and location object and returns a sequence string.
func getFeatureSequence(feature Feature, location Location) (string, error) {
	var sequenceBuffer bytes.Buffer
	var sequenceString string
	parentSequence := feature.ParentSequence.Sequence

	if len(location.SubLocations) == 0 {
		sequenceBuffer.WriteString(parentSequence[location.Start:location.End])
	} else {
		for _, subLocation := range location.SubLocations {
			sequence, _ := getFeatureSequence(feature, subLocation)

			sequenceBuffer.WriteString(sequence)
		}
	}

	// reverse complements resulting string if needed.
	if location.Complement {
		sequenceString = transform.ReverseComplement(sequenceBuffer.String())
	} else {
		sequenceString = sequenceBuffer.String()
	}

	return sequenceString, nil
}

// Read reads a GBK file from path and returns a Genbank struct.
func Read(path string) (Genbank, error) {
	genbankSlice, err := ReadMultiNth(path, 1)
	if err != nil {
		return Genbank{}, err
	}
	genbank := genbankSlice[0]
	return genbank, err
}

// ReadMulti reads a multi Gbk from path and parses it into a slice of Genbank structs.
func ReadMulti(path string) ([]Genbank, error) {
	return ReadMultiNth(path, -1)
}

// ReadMultiNth reads a multi Gbk from path and parses N entries into a slice of Genbank structs.
func ReadMultiNth(path string, count int) ([]Genbank, error) {
	file, err := os.Open(path)
	if err != nil {
		return []Genbank{}, err
	}

	sequence, err := parseMultiNthFn(file, count)
	if err != nil {
		return []Genbank{}, err
	}

	return sequence, nil
}

// Write takes an Genbank list and a path string and writes out a genbank record to that path.
func Write(sequences Genbank, path string) error {
	// build function always returns nil error.
	// This is for API consistency in case we need to
	// add error handling in the future.
	gbk, _ := Build(sequences)

	err := os.WriteFile(path, gbk, 0644)
	return err
}

// WriteMulti takes a slice of Genbank structs and a path string and writes out a multi genbank record to that path.
func WriteMulti(sequences []Genbank, path string) error {
	// buildmulti function always returns nil error.
	// This is for API consistency in case we need to
	// add error handling in the future.
	gbk, _ := BuildMulti(sequences)

	err := os.WriteFile(path, gbk, 0644)
	return err
}

// Build builds a GBK byte slice to be written out to db or file.
func Build(gbk Genbank) ([]byte, error) {
	gbkSlice := []Genbank{gbk}
	multiGBK, err := BuildMulti(gbkSlice)
	return multiGBK, err
}

// BuildMulti builds a MultiGBK byte slice to be written out to db or file.
func BuildMulti(sequences []Genbank) ([]byte, error) {
	var gbkString bytes.Buffer
	for _, sequence := range sequences {
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
		gbkString.WriteString(locusString)

		// building other standard meta features
		definitionString := buildMetaString("DEFINITION", sequence.Meta.Definition)
		gbkString.WriteString(definitionString)

		accessionString := buildMetaString("ACCESSION", sequence.Meta.Accession)
		gbkString.WriteString(accessionString)

		versionString := buildMetaString("VERSION", sequence.Meta.Version)
		gbkString.WriteString(versionString)

		keywordsString := buildMetaString("KEYWORDS", sequence.Meta.Keywords)
		gbkString.WriteString(keywordsString)

		sourceString := buildMetaString("SOURCE", sequence.Meta.Source)
		gbkString.WriteString(sourceString)

		organismString := buildMetaString("  ORGANISM", sequence.Meta.Organism)
		gbkString.WriteString(organismString)

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
			gbkString.WriteString(buildMetaString("", taxonomyString.String()))
		}

		// building references
		// TODO: could use reflection to get keys and make more general.
		for referenceIndex, reference := range sequence.Meta.References {
			referenceString := buildMetaString("REFERENCE", fmt.Sprintf("%d  %s", referenceIndex+1, reference.Range))
			gbkString.WriteString(referenceString)

			if reference.Authors != "" {
				authorsString := buildMetaString("  AUTHORS", reference.Authors)
				gbkString.WriteString(authorsString)
			}

			if reference.Title != "" {
				titleString := buildMetaString("  TITLE", reference.Title)
				gbkString.WriteString(titleString)
			}

			if reference.Journal != "" {
				journalString := buildMetaString("  JOURNAL", reference.Journal)
				gbkString.WriteString(journalString)
			}

			if reference.PubMed != "" {
				pubMedString := buildMetaString("  PUBMED", reference.PubMed)
				gbkString.WriteString(pubMedString)
			}
			if reference.Consortium != "" {
				consrtmString := buildMetaString("  CONSRTM", reference.Consortium)
				gbkString.WriteString(consrtmString)
			}
		}

		// building other meta fields that are catch all
		otherKeys := make([]string, 0, len(sequence.Meta.Other))
		for key := range sequence.Meta.Other {
			otherKeys = append(otherKeys, key)
		}

		for _, otherKey := range otherKeys {
			otherString := buildMetaString(otherKey, sequence.Meta.Other[otherKey])
			gbkString.WriteString(otherString)
		}

		// start writing features section.
		gbkString.WriteString("FEATURES             Location/Qualifiers\n")
		for _, feature := range sequence.Features {
			gbkString.WriteString(BuildFeatureString(feature))
		}

		if len(sequence.Meta.BaseCount) > 0 {
			gbkString.WriteString("BASE COUNT    ")
			for _, baseCount := range sequence.Meta.BaseCount {
				gbkString.WriteString(strconv.Itoa(baseCount.Count) + " " + baseCount.Base + "   ")
			}
			gbkString.WriteString("\n")
		}
		// start writing sequence section.
		gbkString.WriteString("ORIGIN\n")

		// iterate over every character in sequence range.
		for index, base := range sequence.Sequence {
			// if 60th character add newline then whitespace and index number and space before adding next base.
			if index%60 == 0 {
				if index != 0 {
					gbkString.WriteString("\n")
				}
				lineNumberString := strconv.Itoa(index + 1)          // genbank indexes at 1 for some reason
				leadingWhiteSpaceLength := 9 - len(lineNumberString) // <- I wish I was kidding
				for i := 0; i < leadingWhiteSpaceLength; i++ {
					gbkString.WriteString(" ")
				}
				gbkString.WriteString(lineNumberString + " ")
				gbkString.WriteRune(base)
				// if base index is divisible by ten add a space (genbank convention)
			} else if index%10 == 0 {
				gbkString.WriteString(" ")
				gbkString.WriteRune(base)
				// else just add the base.
			} else {
				gbkString.WriteRune(base)
			}
		}
		// finish genbank file with "//" on newline (again a genbank convention)
		gbkString.WriteString("\n//\n")
	}

	return gbkString.Bytes(), nil
}

// Parse takes in a reader representing a single gbk/gb/genbank file and parses it into a Genbank struct.
func Parse(r io.Reader) (Genbank, error) {
	genbankSlice, err := parseMultiNthFn(r, 1)

	if err != nil {
		return Genbank{}, err
	}

	return genbankSlice[0], err
}

// ParseMulti takes in a reader representing a multi gbk/gb/genbank file and parses it into a slice of Genbank structs.
func ParseMulti(r io.Reader) ([]Genbank, error) {
	genbankSlice, err := parseMultiNthFn(r, -1)

	if err != nil {
		return []Genbank{}, err
	}

	return genbankSlice, err
}

type parseLoopParameters struct {
	newLocation      bool
	quoteActive      bool
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
	params.feature.Attributes = make(map[string]string)
	params.parseStep = "metadata"
	params.genbankStarted = false
	params.genbank.Meta.Other = make(map[string]string)
}

// ParseMultiNth takes in a reader representing a multi gbk/gb/genbank file and parses the first n records into a slice of Genbank structs.
func ParseMultiNth(r io.Reader, count int) ([]Genbank, error) {
	scanner := bufio.NewScanner(r)
	var genbanks []Genbank

	// Sequence setup

	var parameters parseLoopParameters
	parameters.init()

	// Loop through each line of the file
	for lineNum := 0; scanner.Scan(); lineNum++ {
		// get line from scanner and split it
		line := scanner.Text()
		splitLine := strings.Split(strings.TrimSpace(line), " ")

		prevline := parameters.currentLine
		parameters.currentLine = line
		parameters.prevline = prevline

		// keep scanning until we find the start of the first record
		if !parameters.genbankStarted {
			// We detect the beginning of a new genbank file with "LOCUS"
			locusFlag := strings.Contains(line, "LOCUS")

			if locusFlag {
				parameters = parseLoopParameters{}
				parameters.init()
				parameters.genbank.Meta.Locus = parseLocus(line)
				parameters.genbankStarted = true
			}
			continue
		}

		switch parameters.parseStep {
		case "metadata":
			// Handle empty lines
			if len(line) == 0 {
				return genbanks, fmt.Errorf("Empty metadata line on line %d", lineNum)
			}

			// If we are currently reading a line, we need to figure out if it is a new meta line.
			if string(line[0]) != " " || parameters.metadataTag == "FEATURES" {
				// If this is true, it means we are beginning a new meta tag. In that case, let's save
				// the older data, and then continue along.
				switch parameters.metadataTag {
				case "DEFINITION":
					parameters.genbank.Meta.Definition = parseMetadata(parameters.metadataData)
				case "ACCESSION":
					parameters.genbank.Meta.Accession = parseMetadata(parameters.metadataData)
				case "VERSION":
					parameters.genbank.Meta.Version = parseMetadata(parameters.metadataData)
				case "KEYWORDS":
					parameters.genbank.Meta.Keywords = parseMetadata(parameters.metadataData)
				case "SOURCE":
					parameters.genbank.Meta.Source, parameters.genbank.Meta.Organism, parameters.genbank.Meta.Taxonomy = getSourceOrganism(parameters.metadataData)
				case "REFERENCE":
					reference, err := parseReferencesFn(parameters.metadataData)
					if err != nil {
						return []Genbank{}, fmt.Errorf("Failed in parsing reference above line %d. Got error: %s", lineNum, err)
					}
					parameters.genbank.Meta.References = append(parameters.genbank.Meta.References, reference)

				case "FEATURES":
					parameters.parseStep = "features"

					// We know that we are now parsing features, so lets initialize our first feature
					parameters.feature.Type = strings.TrimSpace(splitLine[0])
					parameters.feature.Location.GbkLocationString = strings.TrimSpace(splitLine[len(splitLine)-1])
					parameters.newLocation = true

					continue

				default:
					if parameters.metadataTag != "" {
						parameters.genbank.Meta.Other[parameters.metadataTag] = parseMetadata(parameters.metadataData)
					}
				}

				parameters.metadataTag = strings.TrimSpace(splitLine[0])
				parameters.metadataData = []string{strings.TrimSpace(line[len(parameters.metadataTag):])}
			} else {
				parameters.metadataData = append(parameters.metadataData, line)
			}
		case "features":

			baseCountFlag := strings.Contains(line, "BASE COUNT") // example string for BASE COUNT: "BASE COUNT    67070277 a   48055043 c   48111528 g   67244164 t   18475410 n"
			if baseCountFlag {
				fields := strings.Fields(line)
				for countIndex := 2; countIndex < len(fields)-1; countIndex += 2 { // starts at two because we don't want to include "BASE COUNT" in our fields
					count, err := strconv.Atoi(fields[countIndex])
					if err != nil {
						return []Genbank{}, err
					}

					baseCount := BaseCount{
						Base:  fields[countIndex+1],
						Count: count,
					}
					parameters.genbank.Meta.BaseCount = append(parameters.genbank.Meta.BaseCount, baseCount)
				}
				break
			}
			// Switch to sequence parsing
			originFlag := strings.Contains(line, "ORIGIN") // we detect the beginning of the sequence with "ORIGIN"
			if originFlag {
				parameters.parseStep = "sequence"

				// save our completed attribute / qualifier string to the current feature
				if parameters.attributeValue != "" {
					parameters.feature.Attributes[parameters.attribute] = parameters.attributeValue
					parameters.features = append(parameters.features, parameters.feature)
					parameters.attributeValue = ""
					parameters.attribute = ""
					parameters.feature = Feature{}
					parameters.feature.Attributes = make(map[string]string)
				} else {
					parameters.features = append(parameters.features, parameters.feature)
				}

				// add our features to the genbank
				for _, feature := range parameters.features {
					location, err := parseLocation(feature.Location.GbkLocationString)
					if err != nil {
						return []Genbank{}, err
					}
					feature.Location = location
					err = parameters.genbank.AddFeature(&feature)
					if err != nil {
						return []Genbank{}, err
					}
				}
				continue
			} // end sequence parsing flag logic

			// check if current line contains anything but whitespace
			trimmedLine := strings.TrimSpace(line)
			if len(trimmedLine) < 1 {
				continue
			}

			// determine if current line is a new top level feature
			if countLeadingSpaces(parameters.currentLine) < countLeadingSpaces(parameters.prevline) || parameters.prevline == "FEATURES" {
				// save our completed attribute / qualifier string to the current feature
				if parameters.attributeValue != "" {
					parameters.feature.Attributes[parameters.attribute] = parameters.attributeValue
					parameters.features = append(parameters.features, parameters.feature)
					parameters.attributeValue = ""
					parameters.attribute = ""
					parameters.feature = Feature{}
					parameters.feature.Attributes = make(map[string]string)
				}

				// }
				// checks for empty types
				if parameters.feature.Type != "" {
					parameters.features = append(parameters.features, parameters.feature)
				}

				parameters.feature = Feature{}
				parameters.feature.Attributes = make(map[string]string)

				// An initial feature line looks like this: `source          1..2686` with a type separated by its location
				if len(splitLine) < 2 {
					return genbanks, fmt.Errorf("Feature line malformed on line %d. Got line: %s", lineNum, line)
				}
				parameters.feature.Type = strings.TrimSpace(splitLine[0])
				parameters.feature.Location.GbkLocationString = strings.TrimSpace(splitLine[len(splitLine)-1])
				parameters.multiLineFeature = false // without this we can't tell if something is a multiline feature or multiline qualifier
			} else if !strings.Contains(parameters.currentLine, "/") { // current line is continuation of a feature or qualifier (sub-constituent of a feature)
				// if it's a continuation of the current feature, add it to the location
				if !strings.Contains(parameters.currentLine, "\"") && (countLeadingSpaces(parameters.currentLine) > countLeadingSpaces(parameters.prevline) || parameters.multiLineFeature) {
					parameters.feature.Location.GbkLocationString += strings.TrimSpace(line)
					parameters.multiLineFeature = true // without this we can't tell if something is a multiline feature or multiline qualifier
				} else { // it's a continued line of a qualifier
					removeAttributeValueQuotes := strings.Replace(trimmedLine, "\"", "", -1)

					parameters.attributeValue = parameters.attributeValue + removeAttributeValueQuotes
				}
			} else if strings.Contains(parameters.currentLine, "/") { // current line is a new qualifier
				trimmedCurrentLine := strings.TrimSpace(parameters.currentLine)
				if trimmedCurrentLine[0] != '/' { // if we have an exception case, like (adenine(1518)-N(6)/adenine(1519)-N(6))-
					parameters.attributeValue = parameters.attributeValue + trimmedCurrentLine
					continue
				}
				// save our completed attribute / qualifier string to the current feature
				if parameters.attributeValue != "" || parameters.emptyAttribute {
					parameters.feature.Attributes[parameters.attribute] = parameters.attributeValue
					parameters.emptyAttribute = false
				}
				parameters.attributeValue = ""
				splitAttribute := strings.Split(line, "=")
				trimmedSpaceAttribute := strings.TrimSpace(splitAttribute[0])
				removedForwardSlashAttribute := strings.Replace(trimmedSpaceAttribute, "/", "", 1)

				parameters.attribute = removedForwardSlashAttribute

				var removeAttributeValueQuotes string
				if len(splitAttribute) == 1 { // handle case of ` /pseudo `, which has no text
					removeAttributeValueQuotes = ""
					parameters.emptyAttribute = true
				} else { // this is normally triggered
					removeAttributeValueQuotes = strings.Replace(splitAttribute[1], "\"", "", -1)
				}
				parameters.attributeValue = removeAttributeValueQuotes
				parameters.multiLineFeature = false // without this we can't tell if something is a multiline feature or multiline qualifier
			}

		case "sequence":
			if len(line) < 2 { // throw error if line is malformed
				return genbanks, fmt.Errorf("Too short line found while parsing genbank sequence on line %d. Got line: %s", lineNum, line)
			} else if line[0:2] == "//" { // end of sequence
				parameters.genbank.Sequence = parameters.sequenceBuilder.String()

				genbanks = append(genbanks, parameters.genbank)
				parameters.genbankStarted = false
				parameters.sequenceBuilder.Reset()
			} else { // add line to total sequence
				parameters.sequenceBuilder.WriteString(sequenceRegex.ReplaceAllString(line, ""))
			}
		default:
			log.Warnf("Unknown parse step: %s", parameters.parseStep)
			parameters.genbankStarted = false
		}
	}
	return genbanks, nil
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
			position, err := strconv.Atoi(locationString)
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

	qualifierKeys := make([]string, 0, len(feature.Attributes))
	for key := range feature.Attributes {
		qualifierKeys = append(qualifierKeys, key)
	}

	for _, qualifier := range qualifierKeys {
		returnString += generateWhiteSpace(qualifierIndex) + "/" + qualifier + "=\"" + feature.Attributes[qualifier] + "\"\n"
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
