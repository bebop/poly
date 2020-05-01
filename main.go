package main

import (
	"bytes"
	"fmt"
	"io/ioutil"
	"strconv"
	"strings"
)

type Meta struct {
	// gff specific ish
	GffVersion  string
	RegionName  string
	RegionStart int
	RegionEnd   int
	// genbank specific
	Locus           string
	Size            int
	Type            string
	GenbankDivision string
	Date            string
	Definition      string
	Accession       string
	Version         string
	Keywords        string
	Source          string
	SourceOrganism  string
	References      []Reference
}

// genbank specific
type Reference struct {
	Authors []string
	Title   string
	Journal string
	PubMed  string
}

// from https://github.com/blachlylab/gff3/blob/master/gff3.go
type Record struct {
	Seqid      string
	Source     string
	Type       string
	Start      int
	End        int
	Score      float64
	Strand     byte
	Phase      int
	Attributes map[string]string
}

type Sequence struct {
	Sequence    string
	Description string
}

type AnnotatedSequence struct {
	Meta        Meta
	Annotations []Record
	Sequence    Sequence
}

func main() {
	file, _ := ioutil.ReadFile("data/ecoli-mg1655.gff")
	splitFile := strings.Split(string(file), "\n")
	metaString := splitFile[0:2]
	versionString := metaString[0]
	regionStringArray := strings.Split(metaString[1], " ")

	meta := Meta{}
	meta.GffVersion = strings.Split(versionString, " ")[1]
	meta.RegionName = regionStringArray[1]
	meta.RegionStart, _ = strconv.Atoi(regionStringArray[2])
	meta.RegionEnd, _ = strconv.Atoi(regionStringArray[3])
	meta.Size = meta.RegionEnd - meta.RegionStart

	records := []Record{}
	sequence := Sequence{}
	var sequenceBuffer bytes.Buffer
	fastaFlag := false
	for _, line := range splitFile {
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
			record := Record{}
			fields := strings.Split(line, "\t")
			record.Seqid = fields[0]
			record.Source = fields[1]
			record.Type = fields[2]
			record.Start, _ = strconv.Atoi(fields[3])
			record.End, _ = strconv.Atoi(fields[4])
			record.Score, _ = strconv.ParseFloat(fields[5], 64)
			record.Strand = fields[6][0]
			record.Phase, _ = strconv.Atoi(fields[7])
			record.Attributes = make(map[string]string)
			attributes := fields[8]
			var eqIndex int
			for i := strings.Index(attributes, ";"); i > 0; i = strings.Index(attributes, ";") {
				eqIndex = strings.Index(attributes[:i], "=")
				record.Attributes[attributes[:i][:eqIndex]] = attributes[:i][eqIndex+1:]
				attributes = attributes[i+1:]
			}
			records = append(records, record)
		}
	}
	sequence.Sequence = sequenceBuffer.String()
	annotatedSequence := AnnotatedSequence{}
	annotatedSequence.Meta = meta
	annotatedSequence.Annotations = records
	annotatedSequence.Sequence = sequence
	fmt.Println(annotatedSequence.Sequence.Sequence)
}
