package main

import (
	"bytes"
	"io/ioutil"
	"log"
	"regexp"
	"sort"
	"strconv"
	"strings"
)

// ParseGff Takes in a string representing a gffv3 file and parses it into an AnnotatedSequence object.
func ParseGff(gff string) AnnotatedSequence {
	lines := strings.Split(gff, "\n")
	metaString := lines[0:2]
	versionString := metaString[0]
	regionStringArray := strings.Split(metaString[1], " ")

	meta := Meta{}
	meta.GffVersion = strings.Split(versionString, " ")[1]
	meta.Name = regionStringArray[1] // Formally region name, but changed to name here for generality/interoperability.
	meta.RegionStart, _ = strconv.Atoi(regionStringArray[2])
	meta.RegionEnd, _ = strconv.Atoi(regionStringArray[3])
	meta.Size = meta.RegionEnd - meta.RegionStart

	records := []Feature{}
	sequence := Sequence{}
	var sequenceBuffer bytes.Buffer
	fastaFlag := false
	for _, line := range lines {
		if line == "##FASTA" {
			fastaFlag = true
		} else if len(line) == 0 {
			continue
		} else if line[0:2] == "##" {
			continue
		} else if fastaFlag == true && line[0:1] != ">" {
			// sequence.Sequence = sequence.Sequence + line
			sequenceBuffer.WriteString(line)
		} else if fastaFlag == true && line[0:1] == ">" {
			sequence.Description = line
		} else {
			record := Feature{}
			fields := strings.Split(line, "\t")
			record.Name = fields[0]
			record.Source = fields[1]
			record.Type = fields[2]
			record.Start, _ = strconv.Atoi(fields[3])
			record.End, _ = strconv.Atoi(fields[4])
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
			records = append(records, record)
		}
	}
	sequence.Sequence = sequenceBuffer.String()
	annotatedSequence := AnnotatedSequence{}
	annotatedSequence.Meta = meta
	annotatedSequence.Features = records
	annotatedSequence.Sequence = sequence

	return annotatedSequence
}

// BuildGff takes an Annotated sequence and returns a byte array representing a gff to be written out.
func BuildGff(annotatedSequence AnnotatedSequence) []byte {
	var gffBuffer bytes.Buffer

	var versionString string
	if annotatedSequence.Meta.GffVersion != "" {
		versionString = "##gff-version " + annotatedSequence.Meta.GffVersion + "\n"
	} else {
		versionString = "##gff-version 3 \n"
	}
	gffBuffer.WriteString(versionString)

	var regionString string
	var name string
	var start string
	var end string

	if annotatedSequence.Meta.Name != "" {
		name = annotatedSequence.Meta.Name
	} else if annotatedSequence.Meta.Locus.Name != "" {
		name = annotatedSequence.Meta.Locus.Name
	} else if annotatedSequence.Meta.Accession != "" {
		name = annotatedSequence.Meta.Accession
	} else {
		name = "unknown"
	}

	if annotatedSequence.Meta.RegionStart != 0 {
		start = strconv.Itoa(annotatedSequence.Meta.RegionStart)
	} else {
		start = "1"
	}

	if annotatedSequence.Meta.RegionEnd != 0 {
		end = strconv.Itoa(annotatedSequence.Meta.RegionEnd)
	} else if annotatedSequence.Meta.Locus.SequenceLength != "" {
		reg, err := regexp.Compile("[^0-9]+")
		if err != nil {
			log.Fatal(err)
		}
		end = reg.ReplaceAllString(annotatedSequence.Meta.Locus.SequenceLength, "")
	} else {
		end = "1"
	}

	regionString = "##sequence-region " + name + " " + start + " " + end + "\n"
	gffBuffer.WriteString(regionString)

	for _, feature := range annotatedSequence.Features {
		var featureString string

		var featureName string
		if feature.Name != "" {
			featureName = feature.Name
		} else {
			featureName = annotatedSequence.Meta.Name
		}

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

		// really really really need to make a genbank parser util for getting start and stop of region.
		featureStart := strconv.Itoa(feature.Start)

		featureEnd := strconv.Itoa(feature.End)
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
		featureString = featureName + TAB + featureSource + TAB + featureType + TAB + featureStart + TAB + featureEnd + TAB + featureScore + TAB + featureStrand + TAB + featurePhase + TAB + featureAttributes + "\n"
		gffBuffer.WriteString(featureString)
	}

	gffBuffer.WriteString("###\n")
	gffBuffer.WriteString("##FASTA\n")
	gffBuffer.WriteString(">" + annotatedSequence.Meta.Name + "\n")

	for letterIndex, letter := range annotatedSequence.Sequence.Sequence {
		letterIndex++
		if letterIndex%70 == 0 && letterIndex != 0 {
			gffBuffer.WriteRune(letter)
			gffBuffer.WriteString("\n")
		} else {
			gffBuffer.WriteRune(letter)
		}
	}
	gffBuffer.WriteString("\n")
	return gffBuffer.Bytes()
}

// ReadGff takes in a filepath for a .gffv3 file and parses it into an Annotated Sequence struct.
func ReadGff(path string) AnnotatedSequence {
	file, err := ioutil.ReadFile(path)
	var annotatedSequence AnnotatedSequence
	if err != nil {
		// return 0, fmt.Errorf("Failed to open file %s for unpack: %s", gzFilePath, err)
	} else {
		annotatedSequence = ParseGff(string(file))
	}
	return annotatedSequence
}

// WriteGff takes an AnnotatedSequence struct and a path string and writes out a gff to that path.
func WriteGff(annotatedSequence AnnotatedSequence, path string) {
	gff := BuildGff(annotatedSequence)
	_ = ioutil.WriteFile(path, gff, 0644)
}
