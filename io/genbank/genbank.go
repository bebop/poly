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
	"io/ioutil"
	"os"
	"regexp"
	"strconv"
	"strings"

	"github.com/TimothyStiles/poly/transform"
	"github.com/mitchellh/go-wordwrap"
)

/******************************************************************************

GBK specific IO related things begin here.

******************************************************************************/

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
	Authors string `json:"authors"`
	Title   string `json:"title"`
	Journal string `json:"journal"`
	PubMed  string `json:"pub_med"`
	Remark  string `json:"remark"`
	Range   string `json:"range"`
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
	Linear           bool   `json:"linear"`
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
			sequence, err := getFeatureSequence(feature, subLocation)
			if err != nil {
				return sequenceBuffer.String(), err
			}
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

// Parse takes in a string representing a gbk/gb/genbank file and parses it into an Sequence object.
func Parse(r io.Reader) ([]Genbank, error) {
	scanner := bufio.NewScanner(r)
	var sequence = Genbank{}

	// Create meta struct
	var meta = Meta{}

	meta.Other = make(map[string]string)

	// Create features struct
	var features = []Feature{}
	var feature = Feature{}
	feature.Attributes = make(map[string]string)

	// Metadata setup
	var metadataTag string
	var metadataData []string

	// Feature setup
	var newLocation bool
	var quoteActive bool

	var attribute string
	var attributeValue string

	// Sequence setup
	reg, _ := regexp.Compile("[^a-zA-Z]+")
	var sequenceBuilder strings.Builder

	// Sequence breaking
	parseStep := "metadata" // we can be parsing metadata, features, or sequence

	genbankStarted := false
	for lineNum := 0; scanner.Scan(); lineNum++ {
		line := scanner.Text()
		splitLine := strings.Split(strings.TrimSpace(line), " ")
		if !genbankStarted {
			// We detect the beginning of a new genbank file with "LOCUS"
			if line[:5] == "LOCUS" {
				meta.Locus = parseLocus(line)
				genbankStarted = true
				continue
			}
			// Otherwise, we continue scanning
			continue
		}
		switch parseStep {
		case "metadata":
			// Handle empty lines
			if len(line) == 0 {
				return []Genbank{}, fmt.Errorf("Empty metadata line on line %d", lineNum)
			}

			// If we are currently reading a line, we need to figure out if it is a new meta line.
			if string(line[0]) != " " || metadataTag == "FEATURES" {
				// If this is true, it means we are beginning a new meta tag. In that case, let's save
				// the older data, and then continue along.
				switch metadataTag {
				case "DEFINITION":
					meta.Definition = parseMetadata(metadataData)
				case "ACCESSION":
					meta.Accession = parseMetadata(metadataData)
				case "VERSION":
					meta.Version = parseMetadata(metadataData)
				case "KEYWORDS":
					meta.Keywords = parseMetadata(metadataData)
				case "SOURCE":
					meta.Source, meta.Organism, meta.Taxonomy = getSourceOrganism(metadataData)
				case "REFERENCE":
					reference, err := parseReferences(metadataData)
					if err != nil {
						return []Genbank{}, fmt.Errorf("Failed in parsing reference above line %d. Got error: %s", lineNum, err)
					}
					meta.References = append(meta.References, reference)
				case "FEATURES":
					parseStep = "features"
					sequence.Meta = meta
					// We know that we are now parsing features, so lets initialize our first
					// feature
					feature.Type = strings.TrimSpace(splitLine[0])
					feature.Location.GbkLocationString = strings.TrimSpace(splitLine[len(splitLine)-1])
					newLocation = true

					continue
				default:
					if metadataTag != "" {
						meta.Other[metadataTag] = parseMetadata(metadataData)
					}
				}

				metadataTag = strings.TrimSpace(splitLine[0])
				metadataData = []string{strings.TrimSpace(line[len(metadataTag):])}
			} else {
				metadataData = append(metadataData, line)
			}
		case "features":
			// Switch to sequence parsing
			if splitLine[0] == "ORIGIN" {
				parseStep = "sequence"
				// This checks for our initial feature
				if feature.Type != "" {
					feature.Location = parseLocation(feature.Location.GbkLocationString)
					features = append(features, feature)
				}
				for _, feature := range features {
					err := sequence.AddFeature(&feature)
					if err != nil {
						return []Genbank{}, err
					}
				}
				continue
			}

			// If there is no active quote and the line does not have a / and newLocation == false, we know this is a top level feature.
			if !quoteActive && !strings.Contains(line, "/") && !newLocation {
				// Check for empty types
				if feature.Type != "" {
					feature.Location = parseLocation(feature.Location.GbkLocationString)
					features = append(features, feature)
				}
				feature = Feature{}
				feature.Attributes = make(map[string]string)
				// An initial feature line looks like this: `source          1..2686` with a type separated by its location
				if len(splitLine) < 2 {
					return []Genbank{}, fmt.Errorf("Feature line malformed on line %d. Got line: %s", lineNum, line)
				}
				feature.Type = strings.TrimSpace(splitLine[0])
				feature.Location.GbkLocationString = strings.TrimSpace(splitLine[len(splitLine)-1])

				// Set newLocation = true to so that we can check the next line for a multi-line location string
				newLocation = true
				continue
			}

			// If not a newFeature, check if the next line does not contain "/". If it does not, then it is a multi-line location string.
			if !strings.Contains(line, "/") && newLocation {
				feature.Location.GbkLocationString = feature.Location.GbkLocationString + strings.TrimSpace(line)
				continue
			}
			newLocation = false

			// First, let's check if we have a quote active (basically, if a attribute is using multiple lines)
			if quoteActive {
				trimmedLine := strings.TrimSpace(line)
				if len(trimmedLine) < 1 {
					continue
				}
				// If the trimmed line ends, append to the attribute and set quoteActive set to false
				if trimmedLine[len(trimmedLine)-1] == '"' {
					quoteActive = false
					attributeValue = attributeValue + trimmedLine[:len(trimmedLine)-1]
					feature.Attributes[attribute] = attributeValue
					continue
				}

				// If there is still lines to go, just append and continue
				attributeValue = attributeValue + trimmedLine
				continue
			}

			// We know a quote isn't active, that we are not parsing location data, and that we are not at a new top level feature.
			attribute = strings.TrimSpace(strings.Split(line, "=")[0])
			if attribute[0] != '/' {
				return []Genbank{}, fmt.Errorf("Feature attribute does not start with a / on line %d. Got line: %s", lineNum, line)
			}
			attribute = attribute[1:]

			// We have the attribute, now we need attributeValue. We do the following in case '=' is in the attribute value
			index := strings.Index(line, `"`)
			if index == -1 {
				// If `"` is not -1, check for = sign.
				index = strings.Index(line, `=`)
				if index == -1 {
					return []Genbank{}, fmt.Errorf("Double-quote and = not found on feature attribute on line %d. Got line: %s", lineNum, line)
				}
			}
			if len(line) < index+1 {
				return []Genbank{}, fmt.Errorf("Attribute has no data after initial double quote or = on line %d. Got line: %s", lineNum, line)
			}
			attributeValue = strings.TrimSpace(line[index+1:])
			// If true, we have completed the attribute string
			if attributeValue[len(attributeValue)-1] == '"' {
				attributeValue = attributeValue[:len(attributeValue)-1]
				feature.Attributes[attribute] = attributeValue
				continue
			}
			// If false, we'll have to continue with a quoteActive
			quoteActive = true
		case "sequence":
			if len(line) < 2 {
				return []Genbank{}, fmt.Errorf("Too short line found while parsing genbank sequence on line %d. Got line: %s", lineNum, line)
			}
			if line[0:2] == "//" {
				sequence.Sequence = sequenceBuilder.String()
				return []Genbank{sequence}, nil
			}
			sequenceBuilder.WriteString(reg.ReplaceAllString(line, ""))
		}
	}
	return []Genbank{}, fmt.Errorf("No termination of file")
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
	default:
		return fmt.Errorf("ReferenceKey not in [AUTHORS, TITLE, JOURNAL, PUBMED, REMARK]. Got: %s", referenceKey)
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

	basePairRegex, _ := regexp.Compile(` \d* \w{2} `)
	circularRegex, _ := regexp.Compile(` circular `)
	linearRegex, _ := regexp.Compile(` linear `)

	ModificationDateRegex, _ := regexp.Compile(`\d{2}-[A-Z]{3}-\d{4}`)

	locusSplit := strings.Split(strings.TrimSpace(locusString), " ")

	var filteredLocusSplit []string
	for i := range locusSplit {
		if locusSplit[i] != "" {
			filteredLocusSplit = append(filteredLocusSplit, locusSplit[i])
		}
	}

	locus.Name = filteredLocusSplit[1]

	// sequence length and coding
	baseSequenceLength := string(basePairRegex.FindString(locusString))
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

	if linearRegex.Match([]byte(locusString)) {
		locus.Linear = true
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
	locus.ModificationDate = ModificationDateRegex.FindString(locusString)

	return locus
}

// Build builds a GBK string to be written out to db or file.
func Build(sequences []Genbank) ([]byte, error) {
	var gbkString bytes.Buffer
	for _, sequence := range sequences {
		locus := sequence.Meta.Locus
		var shape string

		if locus.Circular {
			shape = "circular"
		} else if locus.Linear {
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

		var taxonomyString strings.Builder
		for i, taxonomyData := range sequence.Meta.Taxonomy {
			taxonomyString.WriteString(taxonomyData)
			if len(taxonomyData) == i+1 {
				taxonomyString.WriteString(".")
			} else {
				taxonomyString.WriteString("; ")
			}
		}
		gbkString.WriteString(buildMetaString("", taxonomyString.String()))

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

// Read reads a Gbk from path and parses into an Annotated sequence struct.
func Read(path string) ([]Genbank, error) {
	file, err := os.Open(path)
	if err != nil {
		return []Genbank{}, err
	}

	sequence, err := Parse(file)
	if err != nil {
		return []Genbank{}, err
	}

	return sequence, nil
}

// Write takes an Genbank list and a path string and writes out a gff to that path.
func Write(sequences []Genbank, path string) error {
	gbk, err := Build(sequences)
	if err != nil {
		return err
	}
	err = ioutil.WriteFile(path, gbk, 0644)
	return err
}

// indeces for random points of interests on a gbk line.
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

func parseLocation(locationString string) Location {
	var location Location
	location.GbkLocationString = locationString
	if !(strings.ContainsAny(locationString, "(")) { // Case checks for simple expression of x..x
		if !(strings.ContainsAny(locationString, ".")) { //Case checks for simple expression x
			position, _ := strconv.Atoi(locationString)
			location = Location{Start: position, End: position}
		} else {
			// to remove FivePrimePartial and ThreePrimePartial indicators from start and end before converting to int.
			partialRegex, _ := regexp.Compile("<|>")
			startEndSplit := strings.Split(locationString, "..")
			start, _ := strconv.Atoi(partialRegex.ReplaceAllString(startEndSplit[0], ""))
			end, _ := strconv.Atoi(partialRegex.ReplaceAllString(startEndSplit[1], ""))
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
				comma := 0
				for i := 1; ParenthesesCount > 0; i++ { // "(" is at 0, so we start at 1
					comma = i
					switch expression[firstInnerParentheses+i] {
					case []byte("(")[0]:
						ParenthesesCount++
					case []byte(")")[0]:
						ParenthesesCount--
					}
				}
				location.SubLocations = append(location.SubLocations, parseLocation(expression[:firstInnerParentheses+comma+1]), parseLocation(expression[2+firstInnerParentheses+comma:]))
			} else { // This is the default join(x..x,x..x)
				for _, numberRange := range strings.Split(expression, ",") {
					location.SubLocations = append(location.SubLocations, parseLocation(numberRange))
				}
			}

		case "complement":
			subLocation := parseLocation(expression)
			subLocation.Complement = true
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

	return location
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
