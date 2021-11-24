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
	"io/ioutil"
	"sort"
	"strconv"
	"strings"

	"lukechampine.com/blake3"

	"github.com/TimothyStiles/poly/io/poly"
)

type Gff struct {
	Meta     Meta
	Features []Feature // will need a GetFeatures interface to standardize
	Sequence string
}

type Feature struct {
	Name           string            `json:"name"`
	Source         string            `json:"source"`
	Type           string            `json:"type"`
	Score          string            `json:"score"`
	Strand         string            `json:"strand"`
	Phase          string            `json:"phase"`
	Attributes     map[string]string `json:"attributes"`
	ParentSequence *Gff              `json:"-"`
	Location       poly.Location
}
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

func (sequence *Gff) AddFeature(feature *Feature) error {
	feature.ParentSequence = sequence
	var featureCopy Feature = *feature
	sequence.Features = append(sequence.Features, featureCopy)
	return nil
}

// Parse Takes in a string representing a gffv3 file and parses it into an Sequence object.
func Parse(file []byte) Gff {

	gffString := string(file)
	gff := Gff{}

	// Add the CheckSum to sequence (blake3)
	gff.Meta.CheckSum = blake3.Sum256(file)

	lines := strings.Split(gffString, "\n")
	metaString := lines[0:2]
	versionString := metaString[0]
	regionStringArray := strings.Split(metaString[1], " ")

	// get name for general meta
	meta := Meta{}
	meta.Name = regionStringArray[1] // Formally region name, but changed to name here for generality/interoperability.

	// get meta info only specific to GFF files
	meta.Version = strings.Split(versionString, " ")[1]
	meta.RegionStart, _ = strconv.Atoi(regionStringArray[2])
	meta.RegionEnd, _ = strconv.Atoi(regionStringArray[3])
	meta.Size = meta.RegionEnd - meta.RegionStart

	var sequenceBuffer bytes.Buffer
	fastaFlag := false
	for _, line := range lines {
		if line == "##FASTA" {
			fastaFlag = true
		} else if len(line) == 0 {
			continue
		} else if line[0:2] == "##" {
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
			record.Location.Start, _ = strconv.Atoi(fields[3])
			record.Location.Start--
			record.Location.End, _ = strconv.Atoi(fields[4])

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
			gff.AddFeature(&record)
		}
	}
	gff.Sequence = sequenceBuffer.String()
	gff.Meta = meta

	return gff
}

// Build takes an Annotated sequence and returns a byte array representing a gff to be written out.
func Build(sequence Gff) []byte {
	var gffBuffer bytes.Buffer

	var versionString string
	if sequence.Meta.Version != "" {
		versionString = "##gff-version " + sequence.Meta.Version + "\n"
	} else {
		versionString = "##gff-version 3 \n"
	}
	gffBuffer.WriteString(versionString)

	var regionString string
	var name string
	var start string
	var end string

	if sequence.Meta.RegionStart != 0 {
		start = strconv.Itoa(sequence.Meta.RegionStart)
	} else {
		start = "1"
	}

	end = strconv.Itoa(sequence.Meta.RegionEnd)

	regionString = "##sequence-region " + name + " " + start + " " + end + "\n"
	gffBuffer.WriteString(regionString)

	for _, feature := range sequence.Features {
		var featureString string
		var featureSource string
		if feature.Source != "" {
			featureSource = feature.Source
		} else {
			featureSource = "feature"
		}

		var featureType string
		if feature.Type != "" {
			featureType = feature.Type
		} else {
			featureType = "unknown"
		}

		// Indexing starts at 1 for gff so we need to shift up from Sequence 0 index.
		featureStart := strconv.Itoa(feature.Location.Start + 1)
		featureEnd := strconv.Itoa(feature.Location.End)

		featureScore := feature.Score
		featureStrand := string(feature.Strand)
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
	return gffBuffer.Bytes()
}

// Read takes in a filepath for a .gffv3 file and parses it into an Annotated poly.Sequence struct.
func Read(path string) Gff {
	file, _ := ioutil.ReadFile(path)
	sequence := Parse(file)
	return sequence
}

// Write takes an poly.Sequence struct and a path string and writes out a gff to that path.
func Write(sequence Gff, path string) {
	gff := Build(sequence)
	_ = ioutil.WriteFile(path, gff, 0644)
}
