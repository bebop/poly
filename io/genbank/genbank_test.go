package genbank

import (
	"io"
	"io/ioutil"
	"os"
	"path/filepath"
	"reflect"
	"strings"
	"testing"

	"github.com/TimothyStiles/poly/transform"
	"github.com/google/go-cmp/cmp"
	"github.com/google/go-cmp/cmp/cmpopts"
)

/******************************************************************************

Gbk/gb/genbank related benchmarks begin here.

******************************************************************************/

func TestGbkIO(t *testing.T) {
	tmpDataDir, err := ioutil.TempDir("", "data-*")
	if err != nil {
		t.Error(err)
	}
	defer os.RemoveAll(tmpDataDir)

	gbk, _ := Read("../../data/puc19.gbk")

	tmpGbkFilePath := filepath.Join(tmpDataDir, "puc19.gbk")
	_ = Write(gbk, tmpGbkFilePath)

	writeTestGbk, _ := Read(tmpGbkFilePath)
	if diff := cmp.Diff(gbk, writeTestGbk, []cmp.Option{cmpopts.IgnoreFields(Feature{}, "ParentSequence")}...); diff != "" {
		t.Errorf("Parsing the output of Build() does not produce the same output as parsing the original file read with Read(). Got this diff:\n%s", diff)
	}

	// Test multiline Genbank features
	pichia, _ := Read("../../data/pichia_chr1_head.gb")
	var multilineOutput string
	for _, feature := range pichia.Features {
		multilineOutput = feature.Location.GbkLocationString
	}

	if multilineOutput != "join(<459260..459456,459556..459637,459685..459739,459810..>460126)" {
		t.Errorf("Failed to parse multiline genbank feature string")
	}
}

func TestGbkLocationStringBuilder(t *testing.T) {
	tmpDataDir, err := ioutil.TempDir("", "data-*")
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
	tmpDataDir, err := ioutil.TempDir("", "data-*")
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

func TestGenbank_AddFeature(t *testing.T) {
	type fields struct {
		Meta     Meta
		Features []Feature
		Sequence string
	}
	type args struct {
		feature *Feature
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
			sequence := &Genbank{
				Meta:     tt.fields.Meta,
				Features: tt.fields.Features,
				Sequence: tt.fields.Sequence,
			}
			if err := sequence.AddFeature(tt.args.feature); (err != nil) != tt.wantErr {
				t.Errorf("Genbank.AddFeature() error = %v, wantErr %v", err, tt.wantErr)
			}
		})
	}
}

func TestFeature_GetSequence(t *testing.T) {
	type fields struct {
		Type                 string
		Description          string
		Attributes           map[string]string
		SequenceHash         string
		SequenceHashFunction string
		Sequence             string
		Location             Location
		ParentSequence       *Genbank
	}
	tests := []struct {
		name    string
		fields  fields
		want    string
		wantErr bool
	}{
		// TODO: Add test cases.
	}
	for _, tt := range tests {
		t.Run(tt.name, func(t *testing.T) {
			feature := Feature{
				Type:                 tt.fields.Type,
				Description:          tt.fields.Description,
				Attributes:           tt.fields.Attributes,
				SequenceHash:         tt.fields.SequenceHash,
				SequenceHashFunction: tt.fields.SequenceHashFunction,
				Sequence:             tt.fields.Sequence,
				Location:             tt.fields.Location,
				ParentSequence:       tt.fields.ParentSequence,
			}
			got, err := feature.GetSequence()
			if (err != nil) != tt.wantErr {
				t.Errorf("Feature.GetSequence() error = %v, wantErr %v", err, tt.wantErr)
				return
			}
			if got != tt.want {
				t.Errorf("Feature.GetSequence() = %v, want %v", got, tt.want)
			}
		})
	}
}

func Test_getFeatureSequence(t *testing.T) {
	type args struct {
		feature  Feature
		location Location
	}
	tests := []struct {
		name    string
		args    args
		want    string
		wantErr bool
	}{
		// TODO: Add test cases.
	}
	for _, tt := range tests {
		t.Run(tt.name, func(t *testing.T) {
			got, err := getFeatureSequence(tt.args.feature, tt.args.location)
			if (err != nil) != tt.wantErr {
				t.Errorf("getFeatureSequence() error = %v, wantErr %v", err, tt.wantErr)
				return
			}
			if got != tt.want {
				t.Errorf("getFeatureSequence() = %v, want %v", got, tt.want)
			}
		})
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

	feature.Description = "Green Flourescent Protein"
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

func TestBuild(t *testing.T) {
	type args struct {
		gbk Genbank
	}
	tests := []struct {
		name    string
		args    args
		want    []byte
		wantErr bool
	}{
		// TODO: Add test cases.
	}
	for _, tt := range tests {
		t.Run(tt.name, func(t *testing.T) {
			got, err := Build(tt.args.gbk)
			if (err != nil) != tt.wantErr {
				t.Errorf("Build() error = %v, wantErr %v", err, tt.wantErr)
				return
			}
			if !reflect.DeepEqual(got, tt.want) {
				t.Errorf("Build() = %v, want %v", got, tt.want)
			}
		})
	}
}

func TestBuildMulti(t *testing.T) {
	type args struct {
		sequence []Genbank
	}
	tests := []struct {
		name    string
		args    args
		want    []byte
		wantErr bool
	}{
		// TODO: Add test cases.
	}
	for _, tt := range tests {
		t.Run(tt.name, func(t *testing.T) {
			got, err := BuildMulti(tt.args.sequence)
			if (err != nil) != tt.wantErr {
				t.Errorf("BuildMulti() error = %v, wantErr %v", err, tt.wantErr)
				return
			}
			if !reflect.DeepEqual(got, tt.want) {
				t.Errorf("BuildMulti() = %v, want %v", got, tt.want)
			}
		})
	}
}

func Test_buildMultiNth(t *testing.T) {
	type args struct {
		sequences []Genbank
		count     int
	}
	tests := []struct {
		name    string
		args    args
		want    []byte
		wantErr bool
	}{
		// TODO: Add test cases.
	}
	for _, tt := range tests {
		t.Run(tt.name, func(t *testing.T) {
			got, err := buildMultiNth(tt.args.sequences, tt.args.count)
			if (err != nil) != tt.wantErr {
				t.Errorf("buildMultiNth() error = %v, wantErr %v", err, tt.wantErr)
				return
			}
			if !reflect.DeepEqual(got, tt.want) {
				t.Errorf("buildMultiNth() = %v, want %v", got, tt.want)
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

func TestReadMultiNth(t *testing.T) {
	type args struct {
		path  string
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
			got, err := ReadMultiNth(tt.args.path, tt.args.count)
			if (err != nil) != tt.wantErr {
				t.Errorf("ReadMultiNth() error = %v, wantErr %v", err, tt.wantErr)
				return
			}
			if !reflect.DeepEqual(got, tt.want) {
				t.Errorf("ReadMultiNth() = %v, want %v", got, tt.want)
			}
		})
	}
}

func TestWrite(t *testing.T) {
	type args struct {
		sequences Genbank
		path      string
	}
	tests := []struct {
		name    string
		args    args
		wantErr bool
	}{
		// TODO: Add test cases.
	}
	for _, tt := range tests {
		t.Run(tt.name, func(t *testing.T) {
			if err := Write(tt.args.sequences, tt.args.path); (err != nil) != tt.wantErr {
				t.Errorf("Write() error = %v, wantErr %v", err, tt.wantErr)
			}
		})
	}
}

func TestWriteMulti(t *testing.T) {
	type args struct {
		sequences []Genbank
		path      string
	}
	tests := []struct {
		name    string
		args    args
		wantErr bool
	}{
		// TODO: Add test cases.
	}
	for _, tt := range tests {
		t.Run(tt.name, func(t *testing.T) {
			if err := WriteMulti(tt.args.sequences, tt.args.path); (err != nil) != tt.wantErr {
				t.Errorf("WriteMulti() error = %v, wantErr %v", err, tt.wantErr)
			}
		})
	}
}

func Test_getSourceOrganism(t *testing.T) {
	type args struct {
		metadataData []string
	}
	tests := []struct {
		name  string
		args  args
		want  string
		want1 string
		want2 []string
	}{
		// TODO: Add test cases.
	}
	for _, tt := range tests {
		t.Run(tt.name, func(t *testing.T) {
			got, got1, got2 := getSourceOrganism(tt.args.metadataData)
			if got != tt.want {
				t.Errorf("getSourceOrganism() got = %v, want %v", got, tt.want)
			}
			if got1 != tt.want1 {
				t.Errorf("getSourceOrganism() got1 = %v, want %v", got1, tt.want1)
			}
			if !reflect.DeepEqual(got2, tt.want2) {
				t.Errorf("getSourceOrganism() got2 = %v, want %v", got2, tt.want2)
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

func TestBuildFeatureString(t *testing.T) {
	type args struct {
		feature Feature
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
			if got := BuildFeatureString(tt.args.feature); got != tt.want {
				t.Errorf("BuildFeatureString() = %v, want %v", got, tt.want)
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
