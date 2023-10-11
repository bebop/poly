package genbank

import (
	"encoding/json"
	"errors"
	"fmt"
	"io"
	"os"
	"path/filepath"
	"strings"
	"testing"

	"reflect"

	"github.com/TimothyStiles/poly/transform"
	"github.com/google/go-cmp/cmp"
	"github.com/google/go-cmp/cmp/cmpopts"
	"github.com/stretchr/testify/assert"
)

/******************************************************************************

Gbk/gb/genbank related benchmarks begin here.

******************************************************************************/

var singleGbkPaths = []string{
	"../../data/t4_intron.gb",
	"../../data/puc19.gbk",
	"../../data/puc19_snapgene.gb",
	"../../data/benchling.gb",
	"../../data/phix174.gb",
	"../../data/sample.gbk",
	// "../../data/pichia_chr1_head.gb",
}

func TestGbkIO(t *testing.T) {
	tmpDataDir, err := os.MkdirTemp("", "data-*")
	if err != nil {
		t.Error(err)
	}
	defer os.RemoveAll(tmpDataDir)

	// test single gbk read, write, build, parse
	for _, gbkPath := range singleGbkPaths {
		gbk, _ := Read(gbkPath)

		tmpGbkFilePath := filepath.Join(tmpDataDir, filepath.Base(gbkPath))
		_ = Write(gbk, tmpGbkFilePath)

		writeTestGbk, _ := Read(tmpGbkFilePath)
		if diff := cmp.Diff(gbk, writeTestGbk, []cmp.Option{cmpopts.IgnoreFields(Feature{}, "ParentSequence")}...); diff != "" {
			t.Errorf("Parsing the output of Build() does not produce the same output as parsing the original file, \"%s\", read with Read(). Got this diff:\n%s", filepath.Base(gbkPath), diff)
		}
	} // end test single gbk read, write, build, parse
}

func TestMultiLineFeatureParse(t *testing.T) {
	pichia, _ := Read("../../data/pichia_chr1_head.gb")
	var multilineOutput string
	for _, feature := range pichia.Features {
		multilineOutput = feature.Location.GbkLocationString
	}

	if multilineOutput != "join(<459260..459456,459556..459637,459685..459739,459810..>460126)" {
		t.Errorf("Failed to parse multiline genbank feature string")
	}
}

func TestMultiGenbankIO(t *testing.T) {
	tmpDataDir, err := os.MkdirTemp("", "data-*")
	if err != nil {
		t.Error(err)
	}
	defer os.RemoveAll(tmpDataDir)

	// Test multiline Genbank features
	gbkPath := "../../data/multiGbk_test.seq"
	multiGbk, _ := ReadMulti(gbkPath)
	tmpGbkFilePath := filepath.Join(tmpDataDir, filepath.Base(gbkPath))
	_ = WriteMulti(multiGbk, tmpGbkFilePath)

	writeTestGbk, _ := ReadMulti(tmpGbkFilePath)

	if diff := cmp.Diff(multiGbk, writeTestGbk, []cmp.Option{cmpopts.IgnoreFields(Feature{}, "ParentSequence")}...); diff != "" {
		t.Errorf("Parsing the output of Build() does not produce the same output as parsing the original file, \"%s\", read with Read(). Got this diff:\n%s", filepath.Base(gbkPath), diff)
	}
}

func TestGbkLocationStringBuilder(t *testing.T) {
	tmpDataDir, err := os.MkdirTemp("", "data-*")
	if err != nil {
		t.Error(err)
	}
	defer os.RemoveAll(tmpDataDir)

	scrubbedGbk, err := Read("../../data/sample.gbk")
	if err != nil {
		t.Error(err)
	}

	// removing gbkLocationString from features to allow testing for gbkLocationBuilder
	for featureIndex := range scrubbedGbk.Features {
		scrubbedGbk.Features[featureIndex].Location.GbkLocationString = ""
	}

	tmpGbkFilePath := filepath.Join(tmpDataDir, "sample.gbk")
	_ = Write(scrubbedGbk, tmpGbkFilePath)

	testInputGbk, _ := Read("../../data/sample.gbk")
	testOutputGbk, _ := Read(tmpGbkFilePath)

	if diff := cmp.Diff(testInputGbk, testOutputGbk, []cmp.Option{cmpopts.IgnoreFields(Feature{}, "ParentSequence")}...); diff != "" {
		t.Errorf("Issue with partial location building. Parsing the output of Build() does not produce the same output as parsing the original file read with Read(). Got this diff:\n%s", diff)
	}
}

func TestGbLocationStringBuilder(t *testing.T) {
	tmpDataDir, err := os.MkdirTemp("", "data-*")
	if err != nil {
		t.Error(err)
	}
	defer os.RemoveAll(tmpDataDir)

	scrubbedGb, _ := Read("../../data/t4_intron.gb")

	// removing gbkLocationString from features to allow testing for gbkLocationBuilder
	for featureIndex := range scrubbedGb.Features {
		scrubbedGb.Features[featureIndex].Location.GbkLocationString = ""
	}

	tmpGbFilePath := filepath.Join(tmpDataDir, "t4_intron_test.gb")
	_ = Write(scrubbedGb, tmpGbFilePath)

	testInputGb, _ := Read("../../data/t4_intron.gb")
	testOutputGb, _ := Read(tmpGbFilePath)

	if diff := cmp.Diff(testInputGb, testOutputGb, []cmp.Option{cmpopts.IgnoreFields(Feature{}, "ParentSequence")}...); diff != "" {
		t.Errorf("Issue with either Join or complement location building. Parsing the output of Build() does not produce the same output as parsing the original file read with Read(). Got this diff:\n%s", diff)
	}
}

func TestPartialLocationParseRegression(t *testing.T) {
	gbk, _ := Read("../../data/sample.gbk")

	for _, feature := range gbk.Features {
		if feature.Location.GbkLocationString == "687..3158>" && (feature.Location.Start != 686 || feature.Location.End != 3158) {
			t.Errorf("Partial location for three prime location parsing has failed. Parsing the output of Build() does not produce the same output as parsing the original file read with Read()")
		}
	}
	gbk, err := Read("../../data/sample.gbk")
	if err != nil {
		t.Errorf("Failed to read sample.gbk. Got err: %s", err)
	}

	for _, feature := range gbk.Features {
		if feature.Location.GbkLocationString == "687..3158>" && (feature.Location.Start != 686 || feature.Location.End != 3158) {
			t.Errorf("Partial location for three prime location parsing has failed. Parsing the output of Build() does not produce the same output as parsing the original file read with Read(). Got location start %d and location end %d. Expected 687..3158>.", feature.Location.Start, feature.Location.End)
		} else if feature.Location.GbkLocationString == "<1..206" && (feature.Location.Start != 0 || feature.Location.End != 206) {
			t.Errorf("Partial location for five prime location parsing has failed. Parsing the output of Build() does not produce the same output as parsing the original file read with Read().")
		}
	}
}

func TestSubLocationStringParseRegression(t *testing.T) {
	location := "join(complement(5306942..5307394),complement(5304401..5305029),complement(5303328..5303393),complement(5301928..5302004))"
	parsedLocation, err := parseLocation(location)
	if err != nil {
		t.Errorf("Failed to parse location string. Got err: %s", err)
	}
	jsonFile, err := os.Open("../../data/parseLocationRegressionTest.json")
	// if we os.Open returns an error then handle it
	if err != nil {
		fmt.Println(err)
	}
	defer jsonFile.Close()

	byteValue, _ := io.ReadAll(jsonFile)
	var testParsedLocation Location
	err = json.Unmarshal(byteValue, &testParsedLocation)
	if err != nil {
		t.Errorf("Failed to unmarshal json. Got err: %s", err)
	}

	if diff := cmp.Diff(parsedLocation, testParsedLocation); diff != "" {
		t.Errorf("Failed to parse sublocation string. Got this diff:\n%s", diff)
	}
}

func TestSnapgeneGenbankRegression(t *testing.T) {
	snapgene, err := Read("../../data/puc19_snapgene.gb")

	if snapgene.Sequence == "" {
		t.Errorf("Parsing snapgene returned an empty string. Got error: %s", err)
	}
}

func TestGetSequenceMethod(t *testing.T) {
	gbk, _ := Read("../../data/t4_intron.gb")

	// Check to see if GetSequence method works on Features struct
	feature, _ := gbk.Features[1].GetSequence()
	seq := "atgagattacaacgccagagcatcaaagattcagaagttagaggtaaatggtattttaatatcatcggtaaagattctgaacttgttgaaaaagctgaacatcttttacgtgatatgggatgggaagatgaatgcgatggatgtcctctttatgaagacggagaaagcgcaggattttggatttaccattctgacgtcgagcagtttaaagctgattggaaaattgtgaaaaagtctgtttga"
	if feature != seq {
		t.Errorf("Feature GetSequence method has failed. Got this:\n%s instead of \n%s", feature, seq)
	}
}

func TestLocationParser(t *testing.T) {
	gbk, _ := Read("../../data/t4_intron.gb")

	// Read 1..243
	feature, _ := gbk.Features[1].GetSequence()
	seq := "atgagattacaacgccagagcatcaaagattcagaagttagaggtaaatggtattttaatatcatcggtaaagattctgaacttgttgaaaaagctgaacatcttttacgtgatatgggatgggaagatgaatgcgatggatgtcctctttatgaagacggagaaagcgcaggattttggatttaccattctgacgtcgagcagtttaaagctgattggaaaattgtgaaaaagtctgtttga"
	if feature != seq {
		t.Errorf("Feature sequence parser has changed on test '1..243'. Got this:\n%s instead of \n%s", feature, seq)
	}

	// Read join(893..1441,2459..2770)
	featureJoin, _ := gbk.Features[6].GetSequence()
	seqJoin := "atgaaacaataccaagatttaattaaagacatttttgaaaatggttatgaaaccgatgatcgtacaggcacaggaacaattgctctgttcggatctaaattacgctgggatttaactaaaggttttcctgcggtaacaactaagaagctcgcctggaaagcttgcattgctgagctaatatggtttttatcaggaagcacaaatgtcaatgatttacgattaattcaacacgattcgttaatccaaggcaaaacagtctgggatgaaaattacgaaaatcaagcaaaagatttaggataccatagcggtgaacttggtccaatttatggaaaacagtggcgtgattttggtggtgtagaccaaattatagaagttattgatcgtattaaaaaactgccaaatgataggcgtcaaattgtttctgcatggaatccagctgaacttaaatatatggcattaccgccttgtcatatgttctatcagtttaatgtgcgtaatggctatttggatttgcagtggtatcaacgctcagtagatgttttcttgggtctaccgtttaatattgcgtcatatgctacgttagttcatattgtagctaagatgtgtaatcttattccaggggatttgatattttctggtggtaatactcatatctatatgaatcacgtagaacaatgtaaagaaattttgaggcgtgaacctaaagagctttgtgagctggtaataagtggtctaccttataaattccgatatctttctactaaagaacaattaaaatatgttcttaaacttaggcctaaagatttcgttcttaacaactatgtatcacaccctcctattaaaggaaagatggcggtgtaa"
	if featureJoin != seqJoin {
		t.Errorf("Feature sequence parser has changed on test 'join(893..1441,2459..2770)'. Got this:\n%s instead of \n%s", featureJoin, seqJoin)
	}

	// Read complement(2791..3054)
	featureComplement, _ := gbk.Features[10].GetSequence()
	seqComplement := "ttattcactacccggcatagacggcccacgctggaataattcgtcatattgtttttccgttaaaacagtaatatcgtagtaacagtcagaagaagttttaactgtggaaattttattatcaaaatactcacgagtcattttatgagtatagtattttttaccataaatggtaataggctgttctggtcctggaacttctaactcgcttgggttaggaagtgtaaaaagaactacaccagaagtatctttaaatcgtaaaatcat"
	if featureComplement != seqComplement {
		t.Errorf("Feature sequence parser has changed on test 'complement(2791..3054)'. Got this:\n%s instead of \n%s", featureComplement, seqComplement)
	}

	// Read join(complement(315..330),complement(339..896))
	// Note: it is known that some software, like Snapgene, assumes that since both strands are in the reverse direction
	// that the first sequence should be appended to the reverse sequence, instead of the second sequence
	// getting appended to the first. Biopython appends the second sequence to the first, and that is logically
	// the most obvious thing to do, so we are implementing it that way.
	featureJoinComplement, _ := gbk.Features[3].GetSequence()
	seqJoinComplement := "ataccaatttaatcattcatttatatactgattccgtaagggttgttacttcatctattttataccaatgcgtttcaaccatttcacgcttgcttatatcatcaagaaaacttgcgtctaattgaactgttgaattaacacgatgccttttaacgatgcgagaaacaactacttcatctgcataaggtaatgcagcatataacagagcaggcccgccaattacacttactttagaattctgatcaagcatagtttcgaatggtgcattagggcttgacacttgaatttcgccgccagaaatgtaagttatatattgctcccaagtaatatagaaatgtgctaaatcgccgtctttagttacaggataatcacgcgcaaggtcacacaccacaatatggctacgaccaggaagtaatgtaggcaatgactggaacgttttagcacccataatcataattgtgccttcagtacgagctttaaaattctggaggtcctttttaactcgtccccatggtaaaccatcacctaaaccgaatgctaattcattaaagccgtcgaccgttttagttggaga"
	if featureJoinComplement != seqJoinComplement {
		t.Errorf("Feature sequence parser has changed on test 'join(complement(315..330),complement(339..896))'. Got this:\n%s instead of \n%s", featureJoinComplement, seqJoinComplement)
	}

	// Read complement(join(893..1098,1101..2770))
	featureComplementJoin, _ := gbk.Features[5].GetSequence()
	seqComplementJoin := "ttacaccgccatctttcctttaataggagggtgtgatacatagttgttaagaacgaaatctttaggcctaagtttaagaacatattttaattgttctttagtagaaagatatcggaatttataaggtagaccacttattaccagctcacaaagctctttaggttcacgcctcaaaatttctttacattgttctacgtgattcatatagatatgagtattaccaccagaaaatatcaaatcccctggaataagattacacatcttagctacaatatgaactaacgtagcatatgacgcaatattaaacggtagcattatgttcagataaggtcgttaatcttaccccggaattatatccagctgcatgtcaccatgcagagcagactatatctccaacttgttaaagcaagttgtctatcgtttcgagtcacttgaccctactccccaaagggatagtcgttaggcatttatgtagaaccaattccatttatcagattttacacgataagtaactaatccagacgaaattttaaaatgtctagctgcatctgctgcacaatcaaaaataaccccatcacatgaaatctttttaatattactaggctttttacctttcatcttttctgatattttagatttagttatgtctgaatgcttatgattaaagaatgaattattttcacctgaacgatttctgcatttactacaagtataagcagaagtttgtatgcgaacaccgcacttacaaaacttatgggtttctggattccaacgcccgtttttacttccgggtttactgtaaagagctttccgaccatcaggtccaagtttaagcatcttagctttaacagtttcagaacgtttcttaataatttcttcttttaatggatgcgtagaacatgtatcaccaaacgttgcatcagcaatattgtatccattaattttagaattaagctctttaatccaaaaattttctcgttcaataatcaaatctttctcatatggaatttcttccaaaatagaacattcaaacacattaccatgtttgttaaaagacctctgaagttttatagaagaatggcatcctttttctaaatctttaaaatgcctcttccatctcttttcaaaatctttagcacttcctacatatactttattgtttaaagtatttttaatctgataaattccgcttttcataaatacctctttaaatatagaagtatttattaaagggcaagtcctacaatttagcacgggattgtctactagagaggttccccgtttagatagattacaagtataagtcaccttatactcaggcctcaattaacccaagaaaacatctactgagcgttgataccactgcaaatccaaatagccattacgcacattaaactgatagaacatatgacaaggcggtaatgccatatatttaagttcagctggattccatgcagaaacaatttgacgcctatcatttggcagttttttaatacgatcaataacttctataatttggtctacaccaccaaaatcacgccactgttttccataaattggaccaagttcaccgctatggtatcctaaatcttttgcttgattttcgtaattttcatcccagactgttttgccttggattaacgaatcgtgttgaattaatcgtaaatcatacatttgtgcttcctgataaaaaccatattagctcagcaatgcaagctttccaggcgagcttcttagttgttaccgcaggaaaacctttagttaaatcccagcgtaatttagatccgaacagagcaattgttcctgtgcctgtacgatcatcggtttcataaccattttcaaaaatgtctttaattaaatcttggtattgtttcat"
	if featureComplementJoin != seqComplementJoin {
		t.Errorf("Feature sequence parser has changed on test 'complement(join(893..1098,1101..2770))'. Got this:\n%s instead of \n%s", featureComplementJoin, seqComplementJoin)
	}
}

func TestGenbankNewlineParsingRegression(t *testing.T) {
	gbk, _ := Read("../../data/puc19.gbk")

	for _, feature := range gbk.Features {
		if feature.Location.Start == 410 && feature.Location.End == 1750 && feature.Type == "CDS" {
			if feature.Attributes["product"] != "chromosomal replication initiator informational ATPase" {
				t.Errorf("Newline parsing has failed.")
			}
			break
		}
	}
}

func BenchmarkRead(b *testing.B) {
	for i := 0; i < b.N; i++ {
		_, _ = Read("../../data/bsub.gbk")
	}
}

func BenchmarkRead1(b *testing.B)     { BenchmarkRead(b) }
func BenchmarkRead10(b *testing.B)    { BenchmarkRead(b) }
func BenchmarkRead100(b *testing.B)   { BenchmarkRead(b) }
func BenchmarkRead1000(b *testing.B)  { BenchmarkRead(b) }
func BenchmarkRead10000(b *testing.B) { BenchmarkRead(b) }

/******************************************************************************

Gbk/gb/genbank related benchmarks end here.

******************************************************************************/

func TestBenchlingGenbank(t *testing.T) {
	sequence, _ := Read("../../data/benchling.gb")

	if len(sequence.Features) != 17 {
		t.Errorf("Parsing benchling genbank file not returned the correct quantity of features")
	}
}

func TestParse(t *testing.T) {
	type args struct {
		r io.Reader
	}
	tests := []struct {
		name    string
		args    args
		want    Genbank
		wantErr bool
	}{
		// TODO: Add test cases.
		// empty line in genbank meta data
		// {

		// 	name:    "empty line in genbank meta data",
		// 	args:    args{r: strings.NewReader("LOCUS       puc19.gbk               2686 bp    DNA     circular  22-OCT-2019")},
		// 	wantErr: true,
		// },
	}
	for _, tt := range tests {
		t.Run(tt.name, func(t *testing.T) {
			got, err := Parse(tt.args.r)
			if (err != nil) != tt.wantErr {
				t.Errorf("Parse() error = %v, wantErr %v", err, tt.wantErr)
				return
			}
			if !reflect.DeepEqual(got, tt.want) {
				t.Errorf("Parse() = %v, want %v", got, tt.want)
			}
		})
	}
}

func TestParseMulti(t *testing.T) {
	type args struct {
		r io.Reader
	}
	tests := []struct {
		name    string
		args    args
		want    []Genbank
		wantErr bool
	}{
		// TODO: Add test cases.
	}
	for _, tt := range tests {
		t.Run(tt.name, func(t *testing.T) {
			got, err := ParseMulti(tt.args.r)
			if (err != nil) != tt.wantErr {
				t.Errorf("ParseMulti() error = %v, wantErr %v", err, tt.wantErr)
				return
			}
			if !reflect.DeepEqual(got, tt.want) {
				t.Errorf("ParseMulti() = %v, want %v", got, tt.want)
			}
		})
	}
}

// this was hand-written and tests the same as the above suite.
func TestFeature_GetSequence_Legacy(t *testing.T) {
	// This test is a little too complex and contrived for an example function.
	// Essentially, it's testing GetSequence()'s ability to parse and retrieve sequences from complex location structures.
	// This was originally covered in the old package system  it was not covered in the new package system so I decided to include it here.

	// Sequence for greenflourescent protein (GFP) that we're using as test data for this example.
	gfpSequence := "ATGGCTAGCAAAGGAGAAGAACTTTTCACTGGAGTTGTCCCAATTCTTGTTGAATTAGATGGTGATGTTAATGGGCACAAATTTTCTGTCAGTGGAGAGGGTGAAGGTGATGCTACATACGGAAAGCTTACCCTTAAATTTATTTGCACTACTGGAAAACTACCTGTTCCATGGCCAACACTTGTCACTACTTTCTCTTATGGTGTTCAATGCTTTTCCCGTTATCCGGATCATATGAAACGGCATGACTTTTTCAAGAGTGCCATGCCCGAAGGTTATGTACAGGAACGCACTATATCTTTCAAAGATGACGGGAACTACAAGACGCGTGCTGAAGTCAAGTTTGAAGGTGATACCCTTGTTAATCGTATCGAGTTAAAAGGTATTGATTTTAAAGAAGATGGAAACATTCTCGGACACAAACTCGAGTACAACTATAACTCACACAATGTATACATCACGGCAGACAAACAAAAGAATGGAATCAAAGCTAACTTCAAAATTCGCCACAACATTGAAGATGGATCCGTTCAACTAGCAGACCATTATCAACAAAATACTCCAATTGGCGATGGCCCTGTCCTTTTACCAGACAACCATTACCTGTCGACACAATCTGCCCTTTCGAAAGATCCCAACGAAAAGCGTGACCACATGGTCCTTCTTGAGTTTGTAACTGCTGCTGGGATTACACATGGCATGGATGAGCTCTACAAATAA"

	sequenceLength := len(gfpSequence)

	// Splitting the sequence into two parts to make a multi-location feature.
	sequenceFirstHalf := gfpSequence[:sequenceLength/2]
	sequenceSecondHalf := transform.ReverseComplement(gfpSequence[sequenceLength/2:]) // This feature is reverse complemented.

	// rejoining the two halves into a single string where the second half of the sequence is reverse complemented.
	gfpSequenceModified := sequenceFirstHalf + sequenceSecondHalf

	// initialize sequence and feature structs.
	var sequence Genbank
	var feature Feature

	// set the initialized sequence struct's sequence.
	sequence.Sequence = gfpSequenceModified
	// initialize sublocations to be usedin the feature.

	var subLocation Location
	var subLocationReverseComplemented Location

	subLocation.Start = 0
	subLocation.End = sequenceLength / 2

	subLocationReverseComplemented.Start = sequenceLength / 2
	subLocationReverseComplemented.End = sequenceLength
	subLocationReverseComplemented.Complement = true // According to genbank complement means reverse complement. What a country.

	feature.Description = "Green Fluorescent Protein"
	feature.Location.SubLocations = []Location{subLocation, subLocationReverseComplemented}

	// Add the GFP feature to the sequence struct.
	_ = sequence.AddFeature(&feature)

	// get the GFP feature sequence string from the sequence struct.
	featureSequence, _ := feature.GetSequence()

	// check to see if the feature was inserted properly into the sequence.
	if gfpSequence != featureSequence {
		t.Error("Feature sequence was not properly retrieved.")
	}
}

func Test_parseLoopParameters_init(t *testing.T) {
	type fields struct {
		newLocation     bool
		quoteActive     bool
		attribute       string
		attributeValue  string
		sequenceBuilder strings.Builder
		parseStep       string
		genbank         Genbank
		feature         Feature
		features        []Feature
		metadataTag     string
		metadataData    []string
		genbankStarted  bool
	}
	tests := []struct {
		name   string
		fields fields
	}{
		// TODO: Add test cases.
	}
	for _, tt := range tests {
		t.Run(tt.name, func(t *testing.T) {
			params := &parseLoopParameters{
				newLocation:     tt.fields.newLocation,
				quoteActive:     tt.fields.quoteActive,
				attribute:       tt.fields.attribute,
				attributeValue:  tt.fields.attributeValue,
				sequenceBuilder: tt.fields.sequenceBuilder,
				parseStep:       tt.fields.parseStep,
				genbank:         tt.fields.genbank,
				feature:         tt.fields.feature,
				features:        tt.fields.features,
				metadataTag:     tt.fields.metadataTag,
				metadataData:    tt.fields.metadataData,
				genbankStarted:  tt.fields.genbankStarted,
			}
			params.init()
		})
	}
}

func TestParseMultiNth(t *testing.T) {
	type args struct {
		r     io.Reader
		count int
	}
	tests := []struct {
		name    string
		args    args
		want    []Genbank
		wantErr bool
	}{
		// TODO: Add test cases.
	}
	for _, tt := range tests {
		t.Run(tt.name, func(t *testing.T) {
			got, err := ParseMultiNth(tt.args.r, tt.args.count)
			if (err != nil) != tt.wantErr {
				t.Errorf("ParseMultiNth() error = %v, wantErr %v", err, tt.wantErr)
				return
			}
			if !reflect.DeepEqual(got, tt.want) {
				t.Errorf("ParseMultiNth() = %v, want %v", got, tt.want)
			}
		})
	}
}

func Test_parseMetadata(t *testing.T) {
	type args struct {
		metadataData []string
	}
	tests := []struct {
		name string
		args args
		want string
	}{
		// TODO: Add test cases.
	}
	for _, tt := range tests {
		t.Run(tt.name, func(t *testing.T) {
			if got := parseMetadata(tt.args.metadataData); got != tt.want {
				t.Errorf("parseMetadata() = %v, want %v", got, tt.want)
			}
		})
	}
}

func Test_parseReferences(t *testing.T) {
	type args struct {
		metadataData []string
	}
	tests := []struct {
		name    string
		args    args
		want    Reference
		wantErr bool
	}{
		// TODO: Add test cases.
	}
	for _, tt := range tests {
		t.Run(tt.name, func(t *testing.T) {
			got, err := parseReferences(tt.args.metadataData)
			if (err != nil) != tt.wantErr {
				t.Errorf("parseReferences() error = %v, wantErr %v", err, tt.wantErr)
				return
			}
			if !reflect.DeepEqual(got, tt.want) {
				t.Errorf("parseReferences() = %v, want %v", got, tt.want)
			}
		})
	}
}

func TestReference_addKey(t *testing.T) {
	type fields struct {
		Authors string
		Title   string
		Journal string
		PubMed  string
		Remark  string
		Range   string
	}
	type args struct {
		referenceKey   string
		referenceValue string
	}
	tests := []struct {
		name    string
		fields  fields
		args    args
		wantErr bool
	}{
		// TODO: Add test cases.
	}
	for _, tt := range tests {
		t.Run(tt.name, func(t *testing.T) {
			reference := &Reference{
				Authors: tt.fields.Authors,
				Title:   tt.fields.Title,
				Journal: tt.fields.Journal,
				PubMed:  tt.fields.PubMed,
				Remark:  tt.fields.Remark,
				Range:   tt.fields.Range,
			}
			if err := reference.addKey(tt.args.referenceKey, tt.args.referenceValue); (err != nil) != tt.wantErr {
				t.Errorf("Reference.addKey() error = %v, wantErr %v", err, tt.wantErr)
			}
		})
	}
}

func Test_parseLocus(t *testing.T) {
	type args struct {
		locusString string
	}
	tests := []struct {
		name string
		args args
		want Locus
	}{
		// TODO: Add test cases.
	}
	for _, tt := range tests {
		t.Run(tt.name, func(t *testing.T) {
			if got := parseLocus(tt.args.locusString); !reflect.DeepEqual(got, tt.want) {
				t.Errorf("parseLocus() = %v, want %v", got, tt.want)
			}
		})
	}
}

func TestRead(t *testing.T) {
	type args struct {
		path string
	}
	tests := []struct {
		name    string
		args    args
		want    Genbank
		wantErr bool
	}{
		// TODO: Add test cases.
		{
			name: "error on missing file",
			args: args{
				path: "../../afdaljhdfa.txt",
			},
			wantErr: true,
		},
		{
			name: "error on malformed file",
			args: args{
				path: "../../data/malformed_read_test.gbk",
			},
			wantErr: true,
		},
	}
	for _, tt := range tests {
		t.Run(tt.name, func(t *testing.T) {
			got, err := Read(tt.args.path)
			if (err != nil) != tt.wantErr {
				t.Errorf("Read() error = %v, wantErr %v", err, tt.wantErr)
				return
			}
			if !reflect.DeepEqual(got, tt.want) {
				t.Errorf("Read() = %v, want %v", got, tt.want)
			}
		})
	}
}

func TestReadMulti(t *testing.T) {
	type args struct {
		path string
	}
	tests := []struct {
		name    string
		args    args
		want    []Genbank
		wantErr bool
	}{
		// TODO: Add test cases.
	}
	for _, tt := range tests {
		t.Run(tt.name, func(t *testing.T) {
			got, err := ReadMulti(tt.args.path)
			if (err != nil) != tt.wantErr {
				t.Errorf("ReadMulti() error = %v, wantErr %v", err, tt.wantErr)
				return
			}
			if !reflect.DeepEqual(got, tt.want) {
				t.Errorf("ReadMulti() = %v, want %v", got, tt.want)
			}
		})
	}
}

func Test_parseLocation(t *testing.T) {
	type args struct {
		locationString string
	}
	tests := []struct {
		name    string
		args    args
		want    Location
		wantErr bool
	}{
		// TODO: Add test cases.
	}
	for _, tt := range tests {
		t.Run(tt.name, func(t *testing.T) {
			got, err := parseLocation(tt.args.locationString)
			if (err != nil) != tt.wantErr {
				t.Errorf("parseLocation() error = %v, wantErr %v", err, tt.wantErr)
				return
			}
			if !reflect.DeepEqual(got, tt.want) {
				t.Errorf("parseLocation() = %v, want %v", got, tt.want)
			}
		})
	}
}

func Test_buildMetaString(t *testing.T) {
	type args struct {
		name string
		data string
	}
	tests := []struct {
		name string
		args args
		want string
	}{
		// TODO: Add test cases.
	}
	for _, tt := range tests {
		t.Run(tt.name, func(t *testing.T) {
			if got := buildMetaString(tt.args.name, tt.args.data); got != tt.want {
				t.Errorf("buildMetaString() = %v, want %v", got, tt.want)
			}
		})
	}
}

func TestBuildLocationString(t *testing.T) {
	type args struct {
		location Location
	}
	tests := []struct {
		name string
		args args
		want string
	}{
		// TODO: Add test cases.
	}
	for _, tt := range tests {
		t.Run(tt.name, func(t *testing.T) {
			if got := BuildLocationString(tt.args.location); got != tt.want {
				t.Errorf("BuildLocationString() = %v, want %v", got, tt.want)
			}
		})
	}
}

func Test_generateWhiteSpace(t *testing.T) {
	type args struct {
		length int
	}
	tests := []struct {
		name string
		args args
		want string
	}{
		// TODO: Add test cases.
	}
	for _, tt := range tests {
		t.Run(tt.name, func(t *testing.T) {
			if got := generateWhiteSpace(tt.args.length); got != tt.want {
				t.Errorf("generateWhiteSpace() = %v, want %v", got, tt.want)
			}
		})
	}
}

func TestRead_error(t *testing.T) {
	readErr := errors.New("open /tmp/file: no such file or directory")
	oldReadFileFn := readFileFn
	readFileFn = func(filename string) ([]byte, error) {
		return nil, readErr
	}
	defer func() {
		readFileFn = oldReadFileFn
	}()
	_, err := Read("/tmp/file")
	assert.EqualError(t, err, readErr.Error())
}

func TestBuildFeatureString(t *testing.T) {
	feature := Feature{
		Type:        "test type",
		Description: "a description",
		Location: Location{
			GbkLocationString: "gbk location",
		},
	}
	str := BuildFeatureString(feature)
	assert.Equal(t, str, "     test type       gbk location\n")
}

func TestParse_error(t *testing.T) {
	parseMultiErr := errors.New("parse error")
	oldParseMultiNthFn := parseMultiNthFn
	parseMultiNthFn = func(r io.Reader, count int) ([]Genbank, error) {
		return nil, parseMultiErr
	}
	defer func() {
		parseMultiNthFn = oldParseMultiNthFn
	}()
	_, err := Parse(strings.NewReader(""))
	assert.EqualError(t, err, parseMultiErr.Error())

	_, err = ParseMulti(strings.NewReader(""))
	assert.EqualError(t, err, parseMultiErr.Error())
}

func TestParseReferences_error(t *testing.T) {
	parseReferencesErr := errors.New("Failed in parsing reference above line 13. Got error: ")
	oldParseReferencesFn := parseReferencesFn
	parseReferencesFn = func(metadataData []string) (Reference, error) {
		return Reference{}, errors.New("")
	}
	defer func() {
		parseReferencesFn = oldParseReferencesFn
	}()
	file, _ := os.Open("../../data/puc19.gbk")
	_, err := parseMultiNthFn(file, 1)
	assert.EqualError(t, err, parseReferencesErr.Error())
}

func TestIssue303Regression(t *testing.T) {
	seq, _ := Read("../../data/puc19_303_regression.gbk")
	expectedAttribute := "16S rRNA(adenine(1518)-N(6)/adenine(1519)-N(6))-dimethyltransferase"
	for _, feature := range seq.Features {
		if feature.Attributes["locus_tag"] == "JCVISYN3A_0004" && feature.Type == "CDS" {
			if feature.Attributes["product"] != expectedAttribute {
				t.Errorf("Failed to get proper expected attribute. Got: %s Expected: %s", feature.Attributes["product"], expectedAttribute)
			}
		}
		if feature.Attributes["locus_tag"] == "JCVISYN3A_0051" && feature.Type == "CDS" {
			if _, ok := feature.Attributes["pseudo"]; !ok {
				t.Errorf("pseudo should be in attributes")
			}
		}
	}
}

func TestConsortiumRegression(t *testing.T) {
	_, err := Read("../../data/puc19_consrtm.gbk")
	if err != nil {
		t.Errorf("Failed to read consrtm. Got err: %s", err)
	}
}
