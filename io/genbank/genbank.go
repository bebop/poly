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
	"compress/gzip"
	"io"
	"io/ioutil"
	"log"
	"regexp"
	"strconv"
	"strings"

	"github.com/TimothyStiles/poly/transform"
	"github.com/mitchellh/go-wordwrap"
	"lukechampine.com/blake3"
)

/******************************************************************************

GBK specific IO related things begin here.

******************************************************************************/

type Genbank struct {
	Meta     Meta
	Features []Feature
	Sequence string // will be changed and include reader, writer, and byte slice.
}

type Meta struct {
	Date                 string            `json:"date"`
	Definition           string            `json:"definition"`
	Accession            string            `json:"accession"`
	Version              string            `json:"version"`
	Keywords             string            `json:"keywords"`
	Organism             string            `json:"organism"`
	Source               string            `json:"source"`
	Origin               string            `json:"origin"`
	Locus                Locus             `json:"locus"`
	References           []Reference       `json:"references"`
	Other                map[string]string `json:"other"`
	Name                 string            `json:"name"`
	SequenceHash         string            `json:"sequence_hash"`
	SequenceHashFunction string            `json:"hash_function"`
	CheckSum             [32]byte          `json:"checkSum"` // blake3 checksum of the parsed file itself. Useful for if you want to check if incoming genbank/gff files are different.
}

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

// Reference holds information one reference in a Meta struct.
type Reference struct {
	Index   string `json:"index"`
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

func (sequence *Genbank) AddFeature(feature *Feature) error {
	feature.ParentSequence = sequence
	sequence.Features = append(sequence.Features, *feature)
	return nil
}

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
func Parse(file []byte) (Genbank, error) {

	gbk := string(file)
	lines := strings.Split(gbk, "\n")

	// Create meta struct
	meta := Meta{}
	meta.Other = make(map[string]string)

	// Create features struct
	features := []Feature{}

	// Create sequence struct
	sequence := Genbank{}

	// Add the CheckSum to sequence (blake3)
	sequence.Meta.CheckSum = blake3.Sum256(file)

	for numLine := 0; numLine < len(lines); numLine++ {
		line := lines[numLine]
		splitLine := strings.Split(line, " ")
		subLines := lines[numLine+1:]

		// This is to keep the cursor from scrolling to the bottom another time after GetSequence() is called.
		// Break has to be in scope and can't be called within switch statement.
		// Otherwise it will just break the switch which is redundant.
		sequenceBreakFlag := false
		if sequenceBreakFlag {
			break
		}

		switch strings.TrimSpace(splitLine[0]) {

		case "":
			continue
		case "LOCUS":
			meta.Locus = parseLocus(line)
		case "DEFINITION":
			meta.Definition = joinSubLines(splitLine, subLines)
		case "ACCESSION":
			meta.Accession = joinSubLines(splitLine, subLines)
		case "VERSION":
			meta.Version = joinSubLines(splitLine, subLines)
		case "KEYWORDS":
			meta.Keywords = joinSubLines(splitLine, subLines)
		case "SOURCE":
			meta.Source, meta.Organism = getSourceOrganism(splitLine, subLines)
		case "REFERENCE":
			meta.References = append(meta.References, getReference(splitLine, subLines))
			continue
		case "FEATURES":
			features = getFeatures(subLines)
		case "ORIGIN":
			sequence.Sequence = getSequence(subLines)
			sequenceBreakFlag = true
		default:
			if quickMetaCheck(line) {
				key := strings.TrimSpace(splitLine[0])
				meta.Other[key] = joinSubLines(splitLine, subLines)
			}
		}

	}

	// add meta to genbank struct
	sequence.Meta = meta

	// add features to annotated sequence with pointer to annotated sequence in each feature
	for _, feature := range features {
		sequence.AddFeature(&feature)
	}

	return sequence, nil
}

// Build builds a GBK string to be written out to db or file.
func Build(sequence Genbank) ([]byte, error) {
	var gbkString bytes.Buffer
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

	// building references
	// TODO: could use reflection to get keys and make more general.
	for referenceIndex, reference := range sequence.Meta.References {
		referenceData := strconv.Itoa(referenceIndex+1) + "  " + reference.Range
		referenceString := buildMetaString("REFERENCE", referenceData)
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
	gbkString.WriteString("\n//")

	return gbkString.Bytes(), nil
}

// Read reads a Gbk from path and parses into an Annotated sequence struct.
func Read(path string) (Genbank, error) {
	file, err := ioutil.ReadFile(path)
	if err != nil {
		return Genbank{}, err
	}

	sequence, err := Parse(file)
	if err != nil {
		return Genbank{}, err
	}

	return sequence, nil
}

// Write takes an Sequence struct and a path string and writes out a gff to that path.
func Write(sequence Genbank, path string) error {
	gbk, err := Build(sequence)
	if err != nil {
		return err
	}
	err = ioutil.WriteFile(path, gbk, 0644)
	return err
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

// used in feature check functions.
var genbankTopLevelFeatures = []string{
	"LOCUS",
	"DEFINITION",
	"ACCESSION",
	"VERSION",
	"KEYWORDS",
	"SOURCE",
	"REFERENCE",
	"FEATURES",
	"ORIGIN",
}

// indeces for random points of interests on a gbk line.
const metaIndex = 0
const subMetaIndex = 5
const qualifierIndex = 21

// benchling actually uses this space between start of the line and /
const benchlingQualifierIndex = 29

func quickMetaCheck(line string) bool {
	flag := false
	// Without line length check, this function
	// panics on Genbank flat files - KG 19 Dec 2020
	if len(line) == 0 {
		return flag
	}
	if string(line[metaIndex]) != " " && string(line[0:2]) != "//" {
		flag = true
	}
	return flag
}

func quickSubMetaCheck(line string) bool {
	flag := false
	// Without line length check, this function
	// panics on Genbank flat files - KG 19 Dec 2020
	if len(line) == 0 {
		return flag
	}
	if string(line[metaIndex]) == " " && string(line[subMetaIndex]) != " " {
		flag = true
	}
	return flag
}

func quickFeatureCheck(line string) bool {
	flag := false

	if string(line[metaIndex]) == " " && string(line[subMetaIndex]) != " " {
		flag = true
	}
	return flag
}

func quickQualifierCheck(line string) bool {
	flag := false

	if string(line[metaIndex]) == " " && string(line[subMetaIndex]) == " " && (string(line[qualifierIndex]) == "/" || string(line[benchlingQualifierIndex]) == "/") {
		flag = true
	}
	return flag

}

func quickQualifierSubLineCheck(line string) bool {
	flag := false

	if string(line[metaIndex]) == " " && string(line[subMetaIndex]) == " " && string(line[qualifierIndex]) != "/" && string(line[qualifierIndex-1]) == " " {
		flag = true
	}
	return flag
}

// checks for only top level features in genbankTopLevelFeatures array
func topLevelFeatureCheck(featureString string) bool {
	flag := false
	cleanedFeatureString := strings.TrimSpace(featureString)
	for _, feature := range genbankTopLevelFeatures {
		if feature == cleanedFeatureString {
			flag = true
			break
		}
	}
	return flag
}

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

// really important helper function. It finds sublines of a feature and joins them.
func joinSubLines(splitLine, subLines []string) string {
	base := strings.TrimSpace(strings.Join(splitLine[1:], " "))

	for _, subLine := range subLines {
		if !quickMetaCheck(subLine) && !quickSubMetaCheck(subLine) {
			base = strings.TrimSpace(strings.TrimSpace(base) + " " + strings.TrimSpace(subLine))
		} else {
			break
		}
	}
	return base
}

// get organism name and source. Doesn't use joinSubLines.
func getSourceOrganism(splitLine, subLines []string) (string, string) {
	source := strings.TrimSpace(strings.Join(splitLine[1:], " "))
	var organism string
	for numSubLine, subLine := range subLines {
		headString := strings.Split(strings.TrimSpace(subLine), " ")[0]
		if string(subLine[0]) == " " && headString != "ORGANISM" {
			source = strings.TrimSpace(strings.TrimSpace(source) + " " + strings.TrimSpace(subLine))
		} else {
			organismSubLines := subLines[numSubLine+1:]
			organismSplitLine := strings.Split(strings.TrimSpace(subLine), " ")
			organism = joinSubLines(organismSplitLine, organismSubLines)
			break
		}
	}
	return source, organism
}

// gets a single reference. Parses headstring and the joins sub lines based on feature.
func getReference(splitLine, subLines []string) Reference {
	base := strings.TrimSpace(strings.Join(splitLine[1:], " "))
	reference := Reference{}
	reference.Index = strings.Split(base, " ")[0]
	if len(base) > 1 {
		reference.Range = strings.TrimSpace(strings.Join(strings.Split(base, " ")[1:], " "))
	}

	for numSubLine, subLine := range subLines {
		featureSubLines := subLines[numSubLine+1:]
		featureSplitLine := strings.Split(strings.TrimSpace(subLine), " ")
		headString := featureSplitLine[0]
		if topLevelFeatureCheck(headString) {
			break
		}
		switch headString {
		case "AUTHORS":
			reference.Authors = joinSubLines(featureSplitLine, featureSubLines)
		case "TITLE":
			reference.Title = joinSubLines(featureSplitLine, featureSubLines)
		case "JOURNAL":
			reference.Journal = joinSubLines(featureSplitLine, featureSubLines)
		case "PUBMED":
			reference.PubMed = joinSubLines(featureSplitLine, featureSubLines)
		case "REMARK":
			reference.Remark = joinSubLines(featureSplitLine, featureSubLines)
		default:
			break
		}

	}
	return reference
}

func getFeatures(lines []string) []Feature {
	lineIndex := 0
	features := []Feature{}

	// regex to remove quotes and slashes from qualifiers
	reg, _ := regexp.Compile("[\"/\n]+")

	// go through every line.
	for lineIndex < len(lines) {
		line := lines[lineIndex]
		// This is a break to ensure that cursor doesn't go beyond ORIGIN which is the last top level feature.
		// This could pick up random sequence strings that aren't helpful and will mess with parser.
		// DO NOT MOVE/REMOVE WITHOUT CAUSE AND CONSIDERATION
		if quickMetaCheck(line) || !quickFeatureCheck(line) {
			break
		}

		feature := Feature{}

		// split the current line for feature type and location fields.
		splitLine := strings.Split(strings.TrimSpace(line), " ")

		// assign type and location to feature.
		feature.Type = strings.TrimSpace(splitLine[0])
		// feature.GbkLocationString is the string used by GBK to denote location
		feature.Location.GbkLocationString = strings.TrimSpace(splitLine[len(splitLine)-1])

		// Check if the location string is multiple lines.
		nextLineNum := 0
		for {
			nextLineNum++
			nextLine := lines[lineIndex+nextLineNum]
			// Check if the next line is not a qualifier, it is part of the
			// GbkLocation String
			if !strings.Contains(nextLine, "/") {
				feature.Location.GbkLocationString = feature.Location.GbkLocationString + strings.TrimSpace(nextLine)
			} else {
				break
			}
		}
		feature.Location = parseLocation(feature.Location.GbkLocationString)

		// initialize attributes.
		feature.Attributes = make(map[string]string)

		// end of feature declaration line. Bump to next line and begin looking for qualifiers.
		lineIndex++
		line = lines[lineIndex+nextLineNum-1]

		// loop through potential qualifiers. Break if not a qualifier or sub line.
		// Definition of qualifiers here: http://www.insdc.org/files/feature_table.html#3.3
		for {
			// make sure what we're parsing is a qualifier. Break if not.
			// keeping out of normal if else pattern because of phantom brackets that are hard to trace.
			if !quickQualifierCheck(line) {
				break
			}

			qualifier := line
			qualifierKey := strings.TrimSpace(strings.Split(line, "=")[0])

			// end of qualifier declaration line. Bump to next line and begin looking for qualifier sublines.
			lineIndex++
			line = lines[lineIndex]

			// loop through any potential continuing lines of qualifiers. Break if not.
			for {
				// keeping out of normal if else pattern because of phantom brackets that are hard to trace.
				if !quickQualifierSubLineCheck(line) {
					break
				}
				//append to current qualifier
				// qualifier += strings.TrimSpace(line)
				if qualifierKey != "/translation" {
					qualifier += " " + strings.TrimSpace(line)
				} else {
					qualifier += strings.TrimSpace(line)
				}

				// nextline
				lineIndex++
				line = lines[lineIndex]
			}
			//add qualifier to feature.
			attributeSplit := strings.Split(reg.ReplaceAllString(qualifier, ""), "=")
			attributeLabel := strings.TrimSpace(attributeSplit[0])
			var attributeValue string
			// This if-statement is not tested, and panics when run. Not sure why
			// it is still here - KG 19 Dec 2020
			if len(attributeSplit) < 2 {
				attributeValue = ""
			} else {
				attributeValue = strings.TrimSpace(attributeSplit[1])
			}
			feature.Attributes[attributeLabel] = attributeValue
		}

		//append the parsed feature to the features list to be returned.
		features = append(features, feature)

	}
	return features
}

// takes every line after origin feature and removes anything that isn't in the alphabet. Returns sequence string.
func getSequence(subLines []string) string {
	var sequenceBuffer bytes.Buffer
	reg, err := regexp.Compile("[^a-zA-Z]+")
	if err != nil {
		log.Fatal(err)
	}
	for _, subLine := range subLines {
		sequenceBuffer.WriteString(subLine)
	}
	sequence := reg.ReplaceAllString(sequenceBuffer.String(), "")
	return sequence
}

func parseLocation(locationString string) Location {
	var location Location
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
			returnData += generateWhiteSpace(11) + datum + "\n"
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

/******************************************************************************

Genbank Flat specific IO related things begin here.

******************************************************************************/

// ParseMulti parses multiple Genbank files in a byte array to multiple sequences
func ParseMulti(file []byte) []Genbank {
	r := bytes.NewReader(file)
	sequences := make(chan Genbank)
	go ParseConcurrent(r, sequences)

	var outputGenbanks []Genbank
	for sequence := range sequences {
		outputGenbanks = append(outputGenbanks, sequence)
	}
	return outputGenbanks
}

// ParseFlat specifically takes the output of a Genbank Flat file that from
// the genbank ftp dumps. These files have 10 line headers, which are entirely
// removed
func ParseFlat(file []byte) []Genbank {
	r := bytes.NewReader(file)
	sequences := make(chan Genbank)
	go ParseFlatConcurrent(r, sequences)
	var outputGenbanks []Genbank
	for sequence := range sequences {
		outputGenbanks = append(outputGenbanks, sequence)
	}
	return outputGenbanks
}

// ReadMulti reads multiple genbank files from a single file
func ReadMulti(path string) []Genbank {
	file, _ := ioutil.ReadFile(path)
	sequences := ParseMulti(file)
	return sequences
}

// ReadFlat reads flat genbank files, like the ones provided by the NCBI FTP server (after decompression)
func ReadFlat(path string) []Genbank {
	file, _ := ioutil.ReadFile(path)
	sequences := ParseFlat(file)
	return sequences
}

// ReadFlatGz reads flat gzip'd genbank files, like the ones provided by the NCBI FTP server
func ReadFlatGz(path string) []Genbank {
	file, _ := ioutil.ReadFile(path)
	rdata := bytes.NewReader(file)
	r, _ := gzip.NewReader(rdata)
	s, _ := ioutil.ReadAll(r)
	sequences := ParseFlat(s)
	return sequences
}

/******************************************************************************

Genbank Flat specific IO related things end here.

******************************************************************************/

/******************************************************************************

Genbank Concurrent specific IO related things begin here.

******************************************************************************/

// ParseConcurrent concurrently parses a given multi-Genbank file in an io.Reader into a channel of Genbank.
func ParseConcurrent(r io.Reader, sequences chan<- Genbank) {
	var gbkStr string
	var gbk Genbank

	// Start a new scanner
	scanner := bufio.NewScanner(r)
	for scanner.Scan() {
		line := scanner.Text()
		if line == "//" {
			gbkStr = gbkStr + "//"
			// Parse the genbank string and send it to the channel
			gbk, _ = Parse([]byte(gbkStr)) // TODO: Ask Keoni how to handle this error
			sequences <- gbk
			// Reset the genbank string
			gbkStr = ""
		} else {
			// Append new lines of the Genbank file to a growing string
			gbkStr = gbkStr + line + "\n"
		}
	}
	close(sequences)
}

// ParseFlatConcurrent concurrently parses a given flat-Genbank file in an io.Reader into a channel of poly.Sequnce.
func ParseFlatConcurrent(r io.Reader, sequences chan<- Genbank) {
	// Start a new reader
	reader := bufio.NewReader(r)
	// Read 10 lines, or the header of a flat file
	// Header data is not needed to parse the Genbank files, though it may contain useful information.
	for i := 0; i < 10; i++ {
		_, _, _ = reader.ReadLine()
	}
	go ParseConcurrent(reader, sequences)
}

/******************************************************************************

Genbank Concurrent specific IO related things end here.

******************************************************************************/
