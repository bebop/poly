package poly

import (
	"bufio"
	"bytes"
	"compress/gzip"
	"encoding/json"
	"io"
	"io/ioutil"
	"log"
	"regexp"
	"sort"
	"strconv"
	"strings"

	"github.com/mitchellh/go-wordwrap"
)

/******************************************************************************

File is structured as so:

Structs:
	Sequence - main struct for sequence handling plus sub structs.

File specific parsers, builders readers, and writers:
	Gff - parser, builder, reader, writer
	JSON- parser, reader, writer
	FASTA - parser, reader, writer
	Gbk/gb/genbank - parser, builder, reader, writer
	Genbank flat / seq.gz - parser, reader


******************************************************************************/

/******************************************************************************

Sequence related structs begin here.

******************************************************************************/

// Meta Holds all the meta information of an Sequence struct.
type Meta struct {
	Name        string            `json:"name"`
	GffVersion  string            `json:"gff_version"`
	RegionStart int               `json:"region_start"`
	RegionEnd   int               `json:"region_end"`
	Size        int               `json:"size"`
	Type        string            `json:"type"`
	Date        string            `json:"date"`
	Definition  string            `json:"definition"`
	Accession   string            `json:"accession"`
	Version     string            `json:"version"`
	Keywords    string            `json:"keywords"`
	Organism    string            `json:"organism"`
	Source      string            `json:"source"`
	Origin      string            `json:"origin"`
	Locus       Locus             `json:"locus"`
	References  []Reference       `json:"references"`
	Other       map[string]string `json:"other"`
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

// Location holds nested location info for sequence region.
type Location struct {
	Start             int        `json:"start"`
	End               int        `json:"end"`
	Complement        bool       `json:"complement"`
	Join              bool       `json:"join"`
	FivePrimePartial  bool       `json:"five_prime_partial"`
	ThreePrimePartial bool       `json:"three_prime_partial"`
	SubLocations      []Location `json:"sub_locations"`
}

// Feature holds a single annotation in a struct. from https://github.com/blachlylab/gff3/blob/master/gff3.go
type Feature struct {
	Name string //Seqid in gff, name in gbk
	//gff specific
	Source               string            `json:"source"`
	Type                 string            `json:"type"`
	Score                string            `json:"score"`
	Strand               string            `json:"strand"`
	Phase                string            `json:"phase"`
	Attributes           map[string]string `json:"attributes"`
	GbkLocationString    string            `json:"gbk_location_string"`
	Sequence             string            `json:"sequence"`
	SequenceLocation     Location          `json:"sequence_location"`
	SequenceHash         string            `json:"sequence_hash"`
	Description          string            `json:"description"`
	SequenceHashFunction string            `json:"hash_function"`
	ParentSequence       *Sequence         `json:"-"`
}

// Sequence holds all sequence information in a single struct.
type Sequence struct {
	Meta                 Meta      `json:"meta"`
	Description          string    `json:"description"`
	SequenceHash         string    `json:"sequence_hash"`
	SequenceHashFunction string    `json:"hash_function"`
	Sequence             string    `json:"sequence"`
	Features             []Feature `json:"features"`
}

// AddFeature is the canonical way to add a Feature into a Sequence struct. Appending a Feature struct directly to Sequence.Feature's will break .GetSequence() method.
func (sequence *Sequence) AddFeature(feature Feature) []Feature {
	feature.ParentSequence = sequence
	sequence.Features = append(sequence.Features, feature)
	return sequence.Features
}

/******************************************************************************

Sequence related structs end here.

******************************************************************************/

/******************************************************************************

GFF specific IO related things begin here.

******************************************************************************/

// ParseGff Takes in a string representing a gffv3 file and parses it into an Sequence object.
func ParseGff(file []byte) Sequence {

	gff := string(file)
	sequence := Sequence{}

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

			// Indexing starts at 1 for gff so we need to shift down for Sequence 0 index.
			record.SequenceLocation.Start, _ = strconv.Atoi(fields[3])
			record.SequenceLocation.Start--
			record.SequenceLocation.End, _ = strconv.Atoi(fields[4])

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
			sequence.AddFeature(record)
		}
	}
	sequence.Sequence = sequenceBuffer.String()
	sequence.Meta = meta

	return sequence
}

// BuildGff takes an Annotated sequence and returns a byte array representing a gff to be written out.
func BuildGff(sequence Sequence) []byte {
	var gffBuffer bytes.Buffer

	var versionString string
	if sequence.Meta.GffVersion != "" {
		versionString = "##gff-version " + sequence.Meta.GffVersion + "\n"
	} else {
		versionString = "##gff-version 3 \n"
	}
	gffBuffer.WriteString(versionString)

	var regionString string
	var name string
	var start string
	var end string

	if sequence.Meta.Name != "" {
		name = sequence.Meta.Name
	} else if sequence.Meta.Locus.Name != "" {
		name = sequence.Meta.Locus.Name
	} else if sequence.Meta.Accession != "" {
		name = sequence.Meta.Accession
	} else {
		name = "unknown"
	}

	if sequence.Meta.RegionStart != 0 {
		start = strconv.Itoa(sequence.Meta.RegionStart)
	} else {
		start = "1"
	}

	if sequence.Meta.RegionEnd != 0 {
		end = strconv.Itoa(sequence.Meta.RegionEnd)
	} else if sequence.Meta.Locus.SequenceLength != "" {
		reg, err := regexp.Compile("[^0-9]+")
		if err != nil {
			log.Fatal(err)
		}
		end = reg.ReplaceAllString(sequence.Meta.Locus.SequenceLength, "")
	} else {
		end = "1"
	}

	regionString = "##sequence-region " + name + " " + start + " " + end + "\n"
	gffBuffer.WriteString(regionString)

	for _, feature := range sequence.Features {
		var featureString string

		var featureName string
		if feature.Name != "" {
			featureName = feature.Name
		} else {
			featureName = sequence.Meta.Locus.Name
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

		// Indexing starts at 1 for gff so we need to shift up from Sequence 0 index.
		featureStart := strconv.Itoa(feature.SequenceLocation.Start + 1)
		featureEnd := strconv.Itoa(feature.SequenceLocation.End)

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

// ReadGff takes in a filepath for a .gffv3 file and parses it into an Annotated Sequence struct.
func ReadGff(path string) Sequence {
	file, _ := ioutil.ReadFile(path)
	var sequence Sequence
	sequence = ParseGff(file)
	return sequence
}

// WriteGff takes an Sequence struct and a path string and writes out a gff to that path.
func WriteGff(sequence Sequence, path string) {
	gff := BuildGff(sequence)
	_ = ioutil.WriteFile(path, gff, 0644)
}

/******************************************************************************

GFF specific IO related things end here.

******************************************************************************/

/******************************************************************************

JSON specific IO related things begin here.

******************************************************************************/

// ParseJSON parses an Sequence JSON file and adds appropriate pointers to struct.
func ParseJSON(file []byte) Sequence {
	var sequence Sequence
	json.Unmarshal([]byte(file), &sequence)
	legacyFeatures := sequence.Features
	sequence.Features = []Feature{}

	for _, feature := range legacyFeatures {
		sequence.AddFeature(feature)
	}
	return sequence
}

// ReadJSON reads an Sequence JSON file.
func ReadJSON(path string) Sequence {
	file, _ := ioutil.ReadFile(path)
	sequence := ParseJSON(file)
	return sequence
}

// WriteJSON writes an Sequence struct out to json.
func WriteJSON(sequence Sequence, path string) {
	file, _ := json.MarshalIndent(sequence, "", " ")
	_ = ioutil.WriteFile(path, file, 0644)
}

/******************************************************************************

JSON specific IO related things end here.

******************************************************************************/

type Fasta struct {
	Name     string `json:"name"`
	Sequence string `json:"sequence"`
}

func ReadFASTAGz(path string, sequences chan<- Fasta) {
	file, _ := ioutil.ReadFile(path)
	rdata := bytes.NewReader(file)
	r, _ := gzip.NewReader(rdata)
	go ParseFASTAGz(r, sequences)
}

func ReadFASTAConcurrent(path string, sequences chan<- Fasta) {
	file, _ := ioutil.ReadFile(path)
	rdata := bytes.NewReader(file)
	go ParseFASTAGz(rdata, sequences)
}

func ParseFASTAGz(r io.Reader, sequences chan<- Fasta) {
	// Initialize necessary variables
	var sequenceLines []string
	var name string
	start := true

	// Start the scanner
	scanner := bufio.NewScanner(r)
	for scanner.Scan() {
		line := scanner.Text()
		// if there's nothing on this line skip this iteration of the loop
		if len(line) == 0 {
			continue
		}
		// if it's a comment skip this line
		if line[0:1] == ";" {
			continue
		}

		// start of a fasta line
		if line[0:1] != ">" {
			sequenceLines = append(sequenceLines, line)
		}
		// Process first line of file
		if line[0:1] == ">" && start == true {
			name = line[1:]
			start = false
		}
		// Process subsequent lines
		if line[0:1] == ">" && start == false {
			sequence := strings.Join(sequenceLines, "")
			newFasta := Fasta{
				Name:     name,
				Sequence: sequence}
			// Reset sequence lines
			sequenceLines = []string{}
			// New name
			name = line[1:]
			sequences <- newFasta
		}
	}
	// Add final sequence in file to channel
	sequence := strings.Join(sequenceLines, "")
	newFasta := Fasta{
		Name:     name,
		Sequence: sequence}
	sequences <- newFasta
	close(sequences)
}

/******************************************************************************

FASTA specific IO related things begin here.

******************************************************************************/

// ParseFASTA parses a Sequence struct from a FASTA file and adds appropriate pointers to the structs.
func ParseFASTA(file []byte) Sequence {
	fasta := string(file)
	var sequence Sequence
	var feature Feature
	var features []Feature
	var sequenceBuffer bytes.Buffer
	var start int
	var end int

	lines := strings.Split(fasta, "\n")
	linesLength := len(lines) - 1

	for lineIndex, line := range lines {

		// if there's nothing on this line skip this iteration of the loop
		if len(line) == 0 {
			continue
		}

		// if it's a comment skip this line
		if line[0:1] == ";" {
			continue
		}

		if line[0:1] == ">" && lineIndex == 0 { // if it's the first description
			feature.Description = line[1:]

		} else if line[0:1] == ">" || lineIndex == linesLength { // if it's a description or the last line

			// if end of file write line to buffer
			if lineIndex == linesLength {
				sequenceBuffer.WriteString(line)
			}

			// setting sequence location
			feature.SequenceLocation.Start = start
			end = len(sequenceBuffer.String())
			feature.SequenceLocation.End = end

			// setting start to end after assigning to location in feature.
			start = end

			// adding new feature to features slice
			features = append(features, feature)

			// resetting feature
			feature = Feature{}

			// if it's the last line
			if lineIndex != linesLength {
				feature.Description = line[1:]
			}

		} else {
			sequenceBuffer.WriteString(line)
		}
	}

	sequence.Sequence = sequenceBuffer.String()

	// add features last so that internal pointer to parent sequence is accurate
	for _, feature := range features {
		sequence.AddFeature(feature)
	}

	return sequence
}

// BuildFASTA builds a FASTA string from a Sequence struct.
func BuildFASTA(sequence Sequence) []byte {
	var fastaBuffer bytes.Buffer
	const maxLineLength = 70

	for featureIndex, feature := range sequence.Features {

		// if there isn't a descriptive comment don't write out feature to fasta file.
		if feature.Description == "" {
			continue
		}
		// write feature comment
		fastaBuffer.WriteString(">" + feature.Description + "\n")

		// range over sequence and add spacing
		for characterIndex, character := range feature.GetSequence() {
			characterIndex++
			if characterIndex%maxLineLength == 0 && characterIndex != 0 {
				fastaBuffer.WriteRune(character)
				fastaBuffer.WriteString("\n")
			} else {
				fastaBuffer.WriteRune(character)
			}
		}

		// if it's the end write new line.
		if featureIndex != len(sequence.Features)-1 {
			fastaBuffer.WriteString("\n\n")
		}
	}

	return fastaBuffer.Bytes()
}

// ReadFASTA reads a Sequence struct from a FASTA file.
func ReadFASTA(path string) Sequence {
	file, _ := ioutil.ReadFile(path)
	sequence := ParseFASTA(file)
	return sequence
}

// WriteFASTA writes a Sequence struct out to FASTA.
func WriteFASTA(sequence Sequence, path string) {
	_ = ioutil.WriteFile(path, BuildFASTA(sequence), 0644)
}

/******************************************************************************

FASTA specific IO related things end here.

******************************************************************************/

/******************************************************************************

GBK specific IO related things begin here.

******************************************************************************/

// ParseGbk takes in a string representing a gbk/gb/genbank file and parses it into an Sequence object.
func ParseGbk(file []byte) Sequence {

	gbk := string(file)
	lines := strings.Split(gbk, "\n")

	// Create meta struct
	meta := Meta{}
	meta.Other = make(map[string]string)

	// Create features struct
	features := []Feature{}

	// Create sequence struct
	sequence := Sequence{}

	for numLine := 0; numLine < len(lines); numLine++ {
		line := lines[numLine]
		splitLine := strings.Split(line, " ")
		subLines := lines[numLine+1:]

		// This is to keep the cursor from scrolling to the bottom another time after GetSequence() is called.
		// Break has to be in scope and can't be called within switch statement.
		// Otherwise it will just break the switch which is redundant.
		sequenceBreakFlag := false
		if sequenceBreakFlag == true {
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
func BuildGbk(sequence Sequence) []byte {
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
func ReadGbk(path string) Sequence {
	file, _ := ioutil.ReadFile(path)
	var sequence Sequence
	sequence = ParseGbk(file)
	return sequence
}

// WriteGbk takes an Sequence struct and a path string and writes out a gff to that path.
func WriteGbk(sequence Sequence, path string) {
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

// used in feature check functions.
var genbankSubLevelFeatures = []string{
	"ORGANISM",
	"AUTHORS",
	"TITLE",
	"JOURNAL",
	"PUBMED",
	"REMARK",
}

// all gene feature types in genbank
var genbankGeneFeatureTypes = []string{
	"assembly_gap",
	"C_region",
	"CDS",
	"centromere",
	"D-loop",
	"D_segment",
	"exon",
	"gap",
	"gene",
	"iDNA",
	"intron",
	"J_segment",
	"mat_peptide",
	"misc_binding",
	"misc_difference",
	"misc_feature",
	"misc_recomb",
	"misc_RNA",
	"misc_structure",
	"mobile_element",
	"modified_base",
	"mRNA",
	"ncRNA",
	"N_region",
	"old_sequence",
	"operon",
	"oriT",
	"polyA_site",
	"precursor_RNA",
	"prim_transcript",
	"primer_bind",
	"propeptide",
	"protein_bind",
	"regulatory",
	"repeat_region",
	"rep_origin",
	"rRNA",
	"S_region",
	"sig_peptide",
	"source",
	"stem_loop",
	"STS",
	"telomere",
	"tmRNA",
	"transit_peptide",
	"tRNA",
	"unsure",
	"V_region",
	"V_segment",
	"variation",
	"3'UTR",
	"5'UTR",
}

// all genbank feature qualifiers
var genbankGeneQualifierTypes = []string{
	"/allele=",
	"/altitude=",
	"/anticodon=",
	"/artificial_location",
	"/bio_material=",
	"/bound_moiety=",
	"/cell_line=",
	"/cell_type=",
	"/chromosome=",
	"/citation=",
	"/clone=",
	"/clone_lib=",
	"/codon_start=",
	"/collected_by=",
	"/collection_date=",
	"/compare=",
	"/country=",
	"/cultivar=",
	"/culture_collection=",
	"/db_xref=",
	"/dev_stage=",
	"/direction=",
	"/EC_number=",
	"/ecotype=",
	"/environmental_sample",
	"/estimated_length=",
	"/exception=",
	"/experiment=",
	"/focus",
	"/frequency=",
	"/function=",
	"/gap_type=",
	"/gene=",
	"/gene_synonym=",
	"/germline",
	"/haplogroup=",
	"/haplotype=",
	"/host=",
	"/identified_by=",
	"/inference=",
	"/isolate=",
	"/isolation_source=",
	"/lab_host=",
	"/lat_lon=",
	"/linkage_evidence=",
	"/locus_tag=",
	"/macronuclear",
	"/map=",
	"/mating_type=",
	"/metagenome_source",
	"/mobile_element_type=",
	"/mod_base=",
	"/mol_type=",
	"/ncRNA_class=",
	"/note=",
	"/number=",
	"/old_locus_tag=",
	"/operon=",
	"/organelle=",
	"/organism=",
	"/partial",
	"/PCR_conditions=",
	"/PCR_primers=",
	"/phenotype=",
	"/plasmid=",
	"/pop_variant=",
	"/product=",
	"/protein_id=",
	"/proviral",
	"/pseudo",
	"/pseudogene=",
	"/rearranged",
	"/replace=",
	"/ribosomal_slippage",
	"/rpt_family=",
	"/rpt_type=",
	"/rpt_unit_range=",
	"/rpt_unit_seq=",
	"/satellite=",
	"/segment=",
	"/serotype=",
	"/serovar=",
	"/sex=",
	"/specimen_voucher=",
	"/standard_name=",
	"/strain=",
	"/sub_clone=",
	"/submitter_seqid=",
	"/sub_species=",
	"/sub_strain=",
	"/tag_peptide=",
	"/tissue_lib=",
	"/tissue_type=",
	"/transgenic",
	"/translation=",
	"/transl_except=",
	"/transl_table=",
	"/trans_splicing",
	"/type_material=",
	"/variety=",
}

// indeces for random points of interests on a gbk line.
const metaIndex = 0
const subMetaIndex = 5
const qualifierIndex = 21

func quickMetaCheck(line string) bool {
	flag := false
	if string(line[metaIndex]) != " " && string(line[0:2]) != "//" {
		flag = true
	}
	return flag
}

func quickSubMetaCheck(line string) bool {
	flag := false

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
		feature.GbkLocationString = strings.TrimSpace(splitLine[len(splitLine)-1])
		feature.SequenceLocation = parseGbkLocation(feature.GbkLocationString)

		// initialize attributes.
		feature.Attributes = make(map[string]string)

		// end of feature declaration line. Bump to next line and begin looking for qualifiers.
		lineIndex++
		line = lines[lineIndex]

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

func parseGbkLocation(locationString string) Location {
	var location Location
	if !(strings.ContainsAny(locationString, "(")) { // Case checks for simple expression of x..x
		// to remove FivePrimePartial and ThreePrimePartial indicators from start and end before converting to int.
		partialRegex, _ := regexp.Compile("<|>")
		startEndSplit := strings.Split(locationString, "..")
		start, _ := strconv.Atoi(partialRegex.ReplaceAllString(startEndSplit[0], ""))
		end, _ := strconv.Atoi(partialRegex.ReplaceAllString(startEndSplit[1], ""))

		location = Location{Start: start - 1, End: end}

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
	if location.Start == 0 && location.End == 0 && location.Join == false && location.Complement == false {
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
func buildGbkLocationString(location Location) string {

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
func buildGbkFeatureString(feature Feature) string {
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
func ParseGbkMulti(file []byte) []Sequence {

	gbk := string(file)
	genbankFiles := strings.SplitAfter(gbk, "//\n")

	//Remove last genbankFile in list. The real file terminates with //, which
	//will be interpreted as an empty genbankFile.
	genbankFiles = genbankFiles[:len(genbankFiles)-1]

	//Iterate through each genbankFile in genbankFiles list and Parse it
	//using the ParseGbk function. Return output.
	var sequences []Sequence
	for _, f := range genbankFiles {
		sequences = append(sequences, ParseGbk([]byte(f)))
	}

	return sequences

}

// ParseGbkFlat specifically takes the output of a Genbank Flat file that from
// the genbank ftp dumps. These files have 10 line headers, which are entirely
// removed
func ParseGbkFlat(file []byte) []Sequence {

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
func ReadGbkMulti(path string) []Sequence {
	file, _ := ioutil.ReadFile(path)
	sequences := ParseGbkMulti(file)
	return sequences
}

// ReadGbkFlat reads flat genbank files, like the ones provided by the NCBI FTP server (after decompression)
func ReadGbkFlat(path string) []Sequence {
	file, _ := ioutil.ReadFile(path)
	sequences := ParseGbkFlat(file)
	return sequences
}

// ReadGbkFlatGz reads flat gzip'd genbank files, like the ones provided by the NCBI FTP server
func ReadGbkFlatGz(path string) []Sequence {
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
