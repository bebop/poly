package genbank

import (
	"bytes"
	"compress/gzip"
	"io/ioutil"
	"log"
	"regexp"
	"strconv"
	"strings"

	"github.com/TimothyStiles/poly"
	"github.com/mitchellh/go-wordwrap"
)

/******************************************************************************

GBK specific IO related things begin here.

******************************************************************************/

// ParseGbk takes in a string representing a gbk/gb/genbank file and parses it into an Sequence object.
func ParseGbk(file []byte) poly.Sequence {

	gbk := string(file)
	lines := strings.Split(gbk, "\n")

	// Create meta struct
	meta := poly.Meta{}
	meta.Other = make(map[string]string)

	// Create features struct
	features := []poly.Feature{}

	// Create sequence struct
	sequence := poly.Sequence{}

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

	// add meta to annotated sequence
	sequence.Meta = meta

	// add features to annotated sequence with pointer to annotated sequence in each feature
	for _, feature := range features {
		sequence.AddFeature(feature)
	}

	return sequence
}

// BuildGbk builds a GBK string to be written out to db or file.
func BuildGbk(sequence poly.Sequence) []byte {
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
		gbkString.WriteString(buildGbkFeatureString(feature))
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

	return gbkString.Bytes()
}

// ReadGbk reads a Gbk from path and parses into an Annotated sequence struct.
func ReadGbk(path string) poly.Sequence {
	file, _ := ioutil.ReadFile(path)
	sequence := ParseGbk(file)
	return sequence
}

// WriteGbk takes an Sequence struct and a path string and writes out a gff to that path.
func WriteGbk(sequence poly.Sequence, path string) {
	gbk := BuildGbk(sequence)
	_ = ioutil.WriteFile(path, gbk, 0644)
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

	if string(line[metaIndex]) == " " && string(line[subMetaIndex]) == " " && string(line[qualifierIndex]) == "/" {
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
func parseLocus(locusString string) poly.Locus {
	locus := poly.Locus{}

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
func getReference(splitLine, subLines []string) poly.Reference {
	base := strings.TrimSpace(strings.Join(splitLine[1:], " "))
	reference := poly.Reference{}
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

func getFeatures(lines []string) []poly.Feature {
	lineIndex := 0
	features := []poly.Feature{}

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

		feature := poly.Feature{}

		// split the current line for feature type and location fields.
		splitLine := strings.Split(strings.TrimSpace(line), " ")

		// assign type and location to feature.
		feature.Type = strings.TrimSpace(splitLine[0])
		// feature.GbkLocationString is the string used by GBK to denote location
		feature.GbkLocationString = strings.TrimSpace(splitLine[len(splitLine)-1])

		// Check if the location string is multiple lines.
		nextLineNum := 0
		for {
			nextLineNum++
			nextLine := lines[lineIndex+nextLineNum]
			// Check if the next line is not a qualifier, it is part of the
			// GbkLocation String
			if !strings.Contains(nextLine, "/") {
				feature.GbkLocationString = feature.GbkLocationString + strings.TrimSpace(nextLine)
			} else {
				break
			}
		}
		feature.SequenceLocation = parseGbkLocation(feature.GbkLocationString)

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

func parseGbkLocation(locationString string) poly.Location {
	var location poly.Location
	if !(strings.ContainsAny(locationString, "(")) { // Case checks for simple expression of x..x
		if !(strings.ContainsAny(locationString, ".")) { //Case checks for simple expression x
			position, _ := strconv.Atoi(locationString)
			location = poly.Location{Start: position, End: position}
		} else {
			// to remove FivePrimePartial and ThreePrimePartial indicators from start and end before converting to int.
			partialRegex, _ := regexp.Compile("<|>")
			startEndSplit := strings.Split(locationString, "..")
			start, _ := strconv.Atoi(partialRegex.ReplaceAllString(startEndSplit[0], ""))
			end, _ := strconv.Atoi(partialRegex.ReplaceAllString(startEndSplit[1], ""))
			location = poly.Location{Start: start - 1, End: end}
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
				location.SubLocations = append(location.SubLocations, parseGbkLocation(expression[:firstInnerParentheses+comma+1]), parseGbkLocation(expression[2+firstInnerParentheses+comma:]))
			} else { // This is the default join(x..x,x..x)
				for _, numberRange := range strings.Split(expression, ",") {
					location.SubLocations = append(location.SubLocations, parseGbkLocation(numberRange))
				}
			}

		case "complement":
			subLocation := parseGbkLocation(expression)
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

// buildGbkLocationString is a recursive function that takes a location object and creates a gbk location string for BuildGbk()
func buildGbkLocationString(location poly.Location) string {

	var locationString string

	if location.Complement {
		location.Complement = false
		locationString = "complement(" + buildGbkLocationString(location) + ")"

	} else if location.Join {
		locationString = "join("
		for _, sublocation := range location.SubLocations {
			locationString += buildGbkLocationString(sublocation) + ","
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

// buildGbkFeatureString is a helper function to build gbk feature strings for BuildGbk()
func buildGbkFeatureString(feature poly.Feature) string {
	whiteSpaceTrailLength := 16 - len(feature.Type) // I wish I was kidding.
	whiteSpaceTrail := generateWhiteSpace(whiteSpaceTrailLength)
	var location string

	if feature.GbkLocationString != "" {
		location = feature.GbkLocationString
	} else {
		location = buildGbkLocationString(feature.SequenceLocation)
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

// ParseGbkMulti parses multiple Genbank files in a byte array to multiple sequences
func ParseGbkMulti(file []byte) []poly.Sequence {

	gbk := string(file)
	genbankFiles := strings.SplitAfter(gbk, "//\n")

	//Remove last genbankFile in list. The real file terminates with //, which
	//will be interpreted as an empty genbankFile.
	genbankFiles = genbankFiles[:len(genbankFiles)-1]

	//Iterate through each genbankFile in genbankFiles list and Parse it
	//using the ParseGbk function. Return output.
	var sequences []poly.Sequence
	for _, f := range genbankFiles {
		sequences = append(sequences, ParseGbk([]byte(f)))
	}

	return sequences

}

// ParseGbkFlat specifically takes the output of a Genbank Flat file that from
// the genbank ftp dumps. These files have 10 line headers, which are entirely
// removed
func ParseGbkFlat(file []byte) []poly.Sequence {

	gbk := string(file)

	// This code removes the header, which is 10 lines long. This is inefficient
	// and gets rid of the data in the header, which may be useful for some
	// application. Header data is not needed to parse the Genbank files, though
	gbkWithoutHeader := []byte(strings.Join(strings.Split(gbk, "\n")[10:], "\n"))

	// Pass gbkWithoutHeader to ParseGbkMulti, which should handle
	// the rest of the parsing just fine
	sequences := ParseGbkMulti(gbkWithoutHeader)
	return sequences
}

// ReadGbkMulti reads multiple genbank files from a single file
func ReadGbkMulti(path string) []poly.Sequence {
	file, _ := ioutil.ReadFile(path)
	sequences := ParseGbkMulti(file)
	return sequences
}

// ReadGbkFlat reads flat genbank files, like the ones provided by the NCBI FTP server (after decompression)
func ReadGbkFlat(path string) []poly.Sequence {
	file, _ := ioutil.ReadFile(path)
	sequences := ParseGbkFlat(file)
	return sequences
}

// ReadGbkFlatGz reads flat gzip'd genbank files, like the ones provided by the NCBI FTP server
func ReadGbkFlatGz(path string) []poly.Sequence {
	file, _ := ioutil.ReadFile(path)
	rdata := bytes.NewReader(file)
	r, _ := gzip.NewReader(rdata)
	s, _ := ioutil.ReadAll(r)
	sequences := ParseGbkFlat(s)
	return sequences
}

/******************************************************************************

Genbank Flat specific IO related things end here.

******************************************************************************/
