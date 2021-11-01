package genbank

import (
	"bufio"
	"fmt"
	"github.com/TimothyStiles/poly/io/poly"
	"io"
	"regexp"
	"strings"
)

func Parse(r io.Reader) (poly.Sequence, error) {
	scanner := bufio.NewScanner(r)

	// Create meta struct
	var meta = poly.Meta{}
	meta.Other = make(map[string]string)

	// Create features struct
	var features = []poly.Feature{}
	var feature = poly.Feature{}
	feature.Attributes = make(map[string]string)

	// Create sequence struct
	var sequence = poly.Sequence{}

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
				metadataTag = strings.TrimSpace(splitLine[0])
				metadataData = []string{strings.TrimSpace(line[len(metadataTag):])}
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
				return poly.Sequence{}, fmt.Errorf("Empty metadata line on line %d", lineNum)
			}

			// If we are currently reading a line, we need to figure out if it is a new meta line.
			if string(line[0]) != " " || metadataTag == "FEATURES" {
				// Handle empty tag lines
				if len(splitLine) < 2 {
					return poly.Sequence{}, fmt.Errorf("Metadata value empty on line %d. Got line: %s", lineNum, line)
				}
				// If this is true, it means we are beginning a new meta tag. In that case, let's save
				// the older data, and then continue along.
				switch metadataTag {
				case "LOCUS":
					meta.Locus = parseLocus(line)
				case "DEFINITION":
					meta.Definition = parseMetadata(metadataData)
				case "ACCESSION":
					meta.Accession = parseMetadata(metadataData)
				case "VERSION":
					meta.Version = parseMetadata(metadataData)
				case "KEYWORDS":
					meta.Keywords = parseMetadata(metadataData)
				case "SOURCE":
					meta.Source, meta.Organism = getSourceOrganism(metadataData)
				case "REFERENCE":
					reference, err := parseReferences(metadataData)
					if err != nil {
						return poly.Sequence{}, fmt.Errorf("Failed in parsing reference above line %d. Got error: %s", lineNum, err)
					}
					meta.References = append(meta.References, reference)
				case "FEATURES":
					parseStep = "features"
					sequence.Meta = meta
				default:
					meta.Other[metadataTag] = parseMetadata(metadataData)
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
					features = append(features, feature)
				}
				for _, feature := range features {
					sequence.AddFeature(&feature)
				}
				continue
			}

			// If there is no active quote and the line does not have a / and newLocation == false, we know this is a top level feature.
			if !quoteActive && !strings.Contains(line, "/") && !newLocation {
				// This checks for our initial feature
				if feature.Type != "" {
					features = append(features, feature)
				}
				feature = poly.Feature{}
				feature.Attributes = make(map[string]string)
				// An initial feature line looks like this: `source          1..2686` with a type separated by its location
				if len(splitLine) < 2 {
					return poly.Sequence{}, fmt.Errorf("Feature line malformed on line %d. Got line: %s", lineNum, line)
				}
				feature.Type = strings.TrimSpace(splitLine[0])
				feature.GbkLocationString = strings.TrimSpace(splitLine[len(splitLine)-1])

				// Set newLocation = true to so that we can check the next line for a multi-line location string
				newLocation = true
				continue
			}

			// If not a newFeature, check if the next line does not contain "/". If it does not, then it is a multi-line location string.
			if !strings.Contains(line, "/") && newLocation {
				feature.GbkLocationString = feature.GbkLocationString + strings.TrimSpace(line)
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
				return poly.Sequence{}, fmt.Errorf("Feature attribute does not start with a / on line %d. Got line: %s", lineNum, line)
			}
			attribute = attribute[1:]

			// We have the attribute, now we need attributeValue. We do the following in case '=' is in the attribute value
			index := strings.Index(line, `"`)
			if index == -1 {
				return poly.Sequence{}, fmt.Errorf("Double-quote not found on feature attribute on line %d. Got line: %s", lineNum, line)
			}
			if len(line) < index+1 {
				return poly.Sequence{}, fmt.Errorf("Attribute has no data after initial double quote on line %d. Got line: %s", lineNum, line)
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
				return poly.Sequence{}, fmt.Errorf("Too short line found while parsing genbank sequence on line %d. Got line: %s", lineNum, line)
			}
			if line[0:2] == "//" {
				sequence.Sequence = sequenceBuilder.String()
				return sequence, nil
			}
			sequenceBuilder.WriteString(reg.ReplaceAllString(line, ""))
		}
	}
	return poly.Sequence{}, fmt.Errorf("No termination of file")
}

func parseMetadata(metadataData []string) string {
	var outputMetadata string
	for _, data := range metadataData {
		outputMetadata = outputMetadata + strings.TrimSpace(strings.Join(strings.Split(data, " ")[1:], " "))
	}
	return outputMetadata
}

type genbankReference poly.Reference

func parseReferences(metadataData []string) (poly.Reference, error) {
	var reference genbankReference
	reference.Index = metadataData[0]
	var err error
	var referenceKey string
	var referenceValue string

	if len(metadataData) == 1 {
		return poly.Reference{}, fmt.Errorf("Got reference with no additional information")
	}

	// Preprocess the first line
	if len(metadataData) > 1 {
		if metadataData[1][3] == ' ' {
			return poly.Reference{}, fmt.Errorf("3rd character in initial GenBank reference has too many spaces. Got: %s", metadataData[1])
		}
		referenceKey = strings.Split(strings.TrimSpace(metadataData[1]), " ")[0]
		referenceValue = metadataData[1][len(referenceKey)+2:]
		// If this is the only data, save it
		if len(metadataData) == 2 {
			err = reference.referenceSwitch(referenceKey, referenceValue)
			if err != nil {
				return poly.Reference{}, err
			}
			return poly.Reference(reference), nil
		}
	}

	for index := 2; index < len(metadataData); index++ {
		// The only way to tell if a reference is top level is if there are only 2 spaces. We use this to know if we should enter our switch.
		if metadataData[index][3] != ' ' {
			err = reference.referenceSwitch(referenceKey, referenceValue)
			if err != nil {
				return poly.Reference{}, err
			}
		} else {
			// Otherwise, simply append the next metadata.
			referenceValue = referenceValue + " " + strings.TrimSpace(metadataData[index])
		}
	}
	err = reference.referenceSwitch(referenceKey, referenceValue)
	if err != nil {
		return poly.Reference{}, err
	}
	return poly.Reference(reference), nil
}

func (ref *genbankReference) referenceSwitch(referenceKey string, referenceValue string) error {
	switch referenceKey {
	case "AUTHORS":
		ref.Authors = referenceValue
	case "TITLE":
		ref.Title = referenceValue
	case "JOURNAL":
		ref.Journal = referenceValue
	case "PUBMED":
		ref.PubMed = referenceValue
	case "REMARK":
		ref.Remark = referenceValue
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

func getSourceOrganism(metadataData []string) (string, string) {
	source := strings.TrimSpace(metadataData[0])
	var organism string
	for iterator := 1; iterator < len(metadataData); iterator++ {
		headString := strings.Split(strings.TrimSpace(metadataData[iterator]), " ")[0]
		if string(metadataData[iterator][0]) == " " && headString != "ORGANISM" {
			source = strings.TrimSpace(strings.TrimSpace(source) + " " + strings.TrimSpace(metadataData[iterator]))
		} else {
			organism = parseMetadata(append(strings.Split(strings.TrimSpace(metadataData[iterator]), " ")[1:], metadataData[iterator+1:]...))
			break
		}
	}
	return source, organism
}
