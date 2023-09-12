/*
Package gff provides gff parsers and writers.

GFF stands for "general feature format". It is an alternative to GenBank for
storing data about genomic sequences. While not often used in synthetic biology
research, it is more commonly used in bioinformatics for digesting features of
genomic sequences.

This package provides a parser and writer to convert between the gff file
format and the more general poly.Sequence struct.
*/
package gff

import (
	"bytes"
	"errors"
	"io"
	"os"
	"sort"
	"strconv"
	"strings"

	"lukechampine.com/blake3"

	"github.com/TimothyStiles/poly/transform"
)

var (
	readAllFn = io.ReadAll
	atoiFn    = strconv.Atoi
	openFn    = os.Open
)

// Gff is a struct that represents a gff file.
type Gff struct {
	Meta     Meta
	Features []Feature // will need a GetFeatures interface to standardize
	Sequence string
}

// Meta holds meta information about a gff file.
type Meta struct {
	Name                 string   `json:"name"`
	Description          string   `json:"description"`
	Version              string   `json:"gff_version"`
	RegionStart          int      `json:"region_start"`
	RegionEnd            int      `json:"region_end"`
	Size                 int      `json:"size"`
	SequenceHash         string   `json:"sequence_hash"`
	SequenceHashFunction string   `json:"hash_function"`
	CheckSum             [32]byte `json:"checkSum"` // blake3 checksum of the parsed file itself. Useful for if you want to check if incoming genbank/gff files are different.
}

// Feature is a struct that represents a feature in a gff file.
type Feature struct {
	Name           string            `json:"name"`
	Source         string            `json:"source"`
	Type           string            `json:"type"`
	Score          string            `json:"score"`
	Strand         string            `json:"strand"`
	Phase          string            `json:"phase"`
	Attributes     map[string]string `json:"attributes"`
	Location       Location          `json:"location"`
	ParentSequence *Gff              `json:"-"`
}

// Location is a struct that represents a location in a gff file.
type Location struct {
	Start             int        `json:"start"`
	End               int        `json:"end"`
	Complement        bool       `json:"complement"`
	Join              bool       `json:"join"`
	FivePrimePartial  bool       `json:"five_prime_partial"`
	ThreePrimePartial bool       `json:"three_prime_partial"`
	SubLocations      []Location `json:"sub_locations"`
}

// AddFeature takes a feature and adds it to the Gff struct.
func (sequence *Gff) AddFeature(feature *Feature) error {
	feature.ParentSequence = sequence
	featureCopy := *feature
	sequence.Features = append(sequence.Features, featureCopy)
	return nil
}

// GetSequence takes a feature and returns a sequence string for that feature.
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

// Parse Takes in a string representing a gffv3 file and parses it into an Sequence object.
func Parse(file io.Reader) (Gff, error) {
	fileBytes, err := readAllFn(file)
	if err != nil {
		return Gff{}, err
	}

	gffString := string(fileBytes)
	gff := Gff{}

	// Add the CheckSum to sequence (blake3)
	gff.Meta.CheckSum = blake3.Sum256(fileBytes)

	lines := strings.Split(gffString, "\n")
	regionStringArray, endOfMetaInfo, err := extractInfoFromField(lines, "##sequence-region")
	metaString := lines[0:endOfMetaInfo]
	versionString := metaString[0]
	if err != nil {
		return Gff{}, err
	}
	// get name for general meta
	meta := Meta{}
	meta.Name = regionStringArray[1] // Formally region name, but changed to name here for generality/interoperability.

	// get meta info only specific to GFF files
	meta.Version = strings.Split(versionString, " ")[1]
	meta.RegionStart, err = atoiFn(regionStringArray[2])
	if err != nil {
		return Gff{}, err
	}
	meta.RegionEnd, err = atoiFn(regionStringArray[3])
	if err != nil {
		return Gff{}, err
	}
	meta.Size = meta.RegionEnd - meta.RegionStart

	var sequenceBuffer bytes.Buffer
	fastaFlag := false
	for _, line := range lines {
		if line == "##FASTA" {
			fastaFlag = true
		} else if len(line) == 0 {
			continue
		} else if line[0:2] == "##" || line[0:2] == "#!" {
			continue
		} else if fastaFlag && line[0:1] != ">" {
			// sequence.Sequence = sequence.Sequence + line
			sequenceBuffer.WriteString(line)
		} else if fastaFlag && line[0:1] == ">" {
			gff.Meta.Description = line
		} else {
			record := Feature{}
			fields := strings.Split(line, "\t")
			record.Name = fields[0]
			record.Source = fields[1]
			record.Type = fields[2]

			// Indexing starts at 1 for gff so we need to shift down for Sequence 0 index.
			record.Location.Start, err = atoiFn(fields[3])
			if err != nil {
				return Gff{}, err
			}

			record.Location.Start--
			record.Location.End, err = atoiFn(fields[4])
			if err != nil {
				return Gff{}, err
			}

			record.Score = fields[5]
			record.Strand = fields[6]
			record.Phase = fields[7]
			record.Attributes = make(map[string]string)
			attributes := fields[8]
			// var eqIndex int
			attributeSlice := strings.Split(attributes, ";")

			for _, attribute := range attributeSlice {
				attributeSplit := strings.Split(attribute, "=")
				key := attributeSplit[0]
				value := attributeSplit[1]
				record.Attributes[key] = value
			}
			err = gff.AddFeature(&record)
			if err != nil {
				return Gff{}, err
			}
		}
	}
	gff.Sequence = sequenceBuffer.String()
	gff.Meta = meta

	return gff, err
}

// regionString takes in the lines array,fieldName that is needed in gff file, and
// returns the region containing fieldName if found
// throws error if not found
func extractInfoFromField(lines []string, fieldName string) ([]string, int, error) {
	index := 0
	endOfMetaInfo := 0
	for lineIndex, line := range lines {
		if strings.Contains(line, "#") {
			if strings.Contains(line, fieldName) {
				index = lineIndex
			}
			continue
		}
		endOfMetaInfo = lineIndex
		break
	}
	if index == 0 && fieldName != "gff-version" {
		return nil, 0, errors.New("the given file does not have any meta information")
	}
	return strings.Split(lines[index], " "), endOfMetaInfo, nil
}

// Build takes an Annotated sequence and returns a byte array representing a gff to be written out.
func Build(sequence Gff) ([]byte, error) {
	var gffBuffer bytes.Buffer

	versionString := "##gff-version 3 \n"
	if sequence.Meta.Version != "" {
		versionString = "##gff-version " + sequence.Meta.Version + "\n"
	}

	gffBuffer.WriteString(versionString)

	name := "Sequence"
	start := "1"
	end := strconv.Itoa(sequence.Meta.RegionEnd)

	if sequence.Meta.Name != "" {
		name = sequence.Meta.Name
	}

	if sequence.Meta.RegionStart != 0 {
		start = strconv.Itoa(sequence.Meta.RegionStart)
	}

	regionString := "##sequence-region " + name + " " + start + " " + end + "\n"
	gffBuffer.WriteString(regionString)

	for _, feature := range sequence.Features {
		var featureString string

		featureSource := "feature"
		if feature.Source != "" {
			featureSource = feature.Source
		}

		featureType := "unknown"
		if feature.Type != "" {
			featureType = feature.Type
		}

		// Indexing starts at 1 for gff so we need to shift up from Sequence 0 index.
		featureStart := strconv.Itoa(feature.Location.Start + 1)
		featureEnd := strconv.Itoa(feature.Location.End)

		featureScore := feature.Score
		featureStrand := feature.Strand
		featurePhase := feature.Phase
		var featureAttributes string

		keys := make([]string, 0, len(feature.Attributes))
		for key := range feature.Attributes {
			keys = append(keys, key)
		}
		sort.Strings(keys)

		for _, key := range keys {
			attributeString := key + "=" + feature.Attributes[key] + ";"
			featureAttributes += attributeString
		}

		if len(featureAttributes) > 0 {
			featureAttributes = featureAttributes[0 : len(featureAttributes)-1]
		}
		TAB := "\t"
		featureString = feature.Name + TAB + featureSource + TAB + featureType + TAB + featureStart + TAB + featureEnd + TAB + featureScore + TAB + featureStrand + TAB + featurePhase + TAB + featureAttributes + "\n"
		gffBuffer.WriteString(featureString)
	}

	gffBuffer.WriteString("###\n")
	gffBuffer.WriteString("##FASTA\n")
	gffBuffer.WriteString(">" + sequence.Meta.Name + "\n")

	for letterIndex, letter := range sequence.Sequence {
		letterIndex++
		if letterIndex%70 == 0 && letterIndex != 0 && letterIndex != sequence.Meta.RegionEnd {
			gffBuffer.WriteRune(letter)
			gffBuffer.WriteString("\n")
		} else {
			gffBuffer.WriteRune(letter)
		}
	}
	gffBuffer.WriteString("\n")
	return gffBuffer.Bytes(), nil
}

// Read takes in a filepath for a .gffv3 file and parses it into an Annotated poly.Sequence struct.
func Read(path string) (Gff, error) {
	file, err := openFn(path)
	if err != nil {
		return Gff{}, err
	}

	sequence, err := Parse(file)
	return sequence, err
}

// Write takes an poly.Sequence struct and a path string and writes out a gff to that path.
func Write(sequence Gff, path string) error {
	gff, _ := Build(sequence)
	return os.WriteFile(path, gff, 0644)
}
