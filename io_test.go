package poly

import (
	"bytes"
	"fmt"
	"io/ioutil"
	"os"
	"path/filepath"
	"strings"
	"testing"

	"github.com/google/go-cmp/cmp"
	"github.com/google/go-cmp/cmp/cmpopts"
	"github.com/pmezard/go-difflib/difflib"
)

/******************************************************************************

File is structured as so:

Gff - io tests, and benchmarks.
Gbk/gb/genbank - benchmarks.
JSON - io tests.
FASTA - fasta tests.

******************************************************************************/

/******************************************************************************

Gff related tests and benchmarks begin here.

******************************************************************************/

func ExampleReadGff() {

	sequence := ReadGff("data/ecoli-mg1655-short.gff")
	fmt.Println(sequence.Meta.Name)
	// Output: U00096.3
}

func ExampleParseGff() {
	file, _ := ioutil.ReadFile("data/ecoli-mg1655-short.gff")
	sequence := ParseGff(file)

	fmt.Println(sequence.Meta.Name)
	// Output: U00096.3
}

func ExampleBuildGff() {

	sequence := ReadGff("data/ecoli-mg1655-short.gff")
	gffBytes := BuildGff(sequence)
	reparsedSequence := ParseGff(gffBytes)

	fmt.Println(reparsedSequence.Meta.Name)
	// Output: U00096.3

}

func ExampleWriteGff() {
	tmpDataDir, err := ioutil.TempDir("", "data-*")
	defer os.RemoveAll(tmpDataDir)
	if err != nil {
		fmt.Println(err.Error())
	}

	sequence := ReadGff("data/ecoli-mg1655-short.gff")

	tmpGffFilePath := filepath.Join(tmpDataDir, "ecoli-mg1655-short.gff")
	WriteGff(sequence, tmpGffFilePath)

	testSequence := ReadGff(tmpGffFilePath)

	fmt.Println(testSequence.Meta.Name)
	// Output: U00096.3
}

// TODO should delete output files.

func TestGffIO(t *testing.T) {
	tmpDataDir, err := ioutil.TempDir("", "data-*")
	defer os.RemoveAll(tmpDataDir)
	if err != nil {
		t.Error(err)
	}

	testInputPath := "data/ecoli-mg1655-short.gff"
	tmpGffFilePath := filepath.Join(tmpDataDir, "ecoli-mg1655-short.gff")

	testSequence := ReadGff(testInputPath)
	WriteGff(testSequence, tmpGffFilePath)

	readTestSequence := ReadGff(tmpGffFilePath)

	if diff := cmp.Diff(testSequence, readTestSequence, cmpopts.IgnoreFields(Feature{}, "ParentSequence")); diff != "" {
		t.Errorf("Parsing the output of BuildGff() does not produce the same output as parsing the original file read with ReadGff(). Got this diff:\n%s", diff)
	}

	original, _ := ioutil.ReadFile(testInputPath)
	builtOutput, _ := ioutil.ReadFile(tmpGffFilePath)
	gffDiff := difflib.UnifiedDiff{
		A:        difflib.SplitLines(string(original)),
		B:        difflib.SplitLines(string(builtOutput)),
		FromFile: testInputPath,
		ToFile:   tmpGffFilePath,
		Context:  3,
	}

	gffDiffText, _ := difflib.GetUnifiedDiffString(gffDiff)

	if gffDiffText != "" {
		t.Errorf("BuildGff() does not output the same file as was input through ReadGff(). Got this diff:\n%s", gffDiffText)
	}

}

func BenchmarkReadGff(b *testing.B) {
	for i := 0; i < b.N; i++ {
		ReadGff("data/ecoli-mg1655-short.gff")
	}
}

func BenchmarkReadGff1(b *testing.B)     { BenchmarkReadGff(b) }
func BenchmarkReadGff10(b *testing.B)    { BenchmarkReadGff(b) }
func BenchmarkReadGff100(b *testing.B)   { BenchmarkReadGff(b) }
func BenchmarkReadGff1000(b *testing.B)  { BenchmarkReadGff(b) }
func BenchmarkReadGff10000(b *testing.B) { BenchmarkReadGff(b) }

/******************************************************************************

Gff related tests and benchmarks end here.

******************************************************************************/

/******************************************************************************

Gbk/gb/genbank related benchmarks begin here.

******************************************************************************/

func ExampleReadGbk() {
	sequence := ReadGbk("data/puc19.gbk")
	fmt.Println(sequence.Meta.Locus.ModificationDate)
	// Output: 22-OCT-2019
}

func ExampleParseGbk() {
	file, _ := ioutil.ReadFile("data/puc19.gbk")
	sequence := ParseGbk(file)

	fmt.Println(sequence.Meta.Locus.ModificationDate)
	// Output: 22-OCT-2019
}

func ExampleBuildGbk() {
	sequence := ReadGbk("data/puc19.gbk")
	gbkBytes := BuildGbk(sequence)
	testSequence := ParseGbk(gbkBytes)

	fmt.Println(testSequence.Meta.Locus.ModificationDate)
	// Output: 22-OCT-2019
}

func ExampleWriteGbk() {
	tmpDataDir, err := ioutil.TempDir("", "data-*")
	defer os.RemoveAll(tmpDataDir)
	if err != nil {
		fmt.Println(err.Error())
	}

	sequence := ReadGbk("data/puc19.gbk")

	tmpGbkFilePath := filepath.Join(tmpDataDir, "puc19.gbk")
	WriteGbk(sequence, tmpGbkFilePath)

	testSequence := ReadGbk(tmpGbkFilePath)

	fmt.Println(testSequence.Meta.Locus.ModificationDate)
	// Output: 22-OCT-2019
}

func TestGbkIO(t *testing.T) {
	tmpDataDir, err := ioutil.TempDir("", "data-*")
	defer os.RemoveAll(tmpDataDir)
	if err != nil {
		t.Error(err)
	}

	gbk := ReadGbk("data/puc19.gbk")

	tmpGbkFilePath := filepath.Join(tmpDataDir, "puc19.gbk")
	WriteGbk(gbk, tmpGbkFilePath)

	writeTestGbk := ReadGbk(tmpGbkFilePath)
	if diff := cmp.Diff(gbk, writeTestGbk, cmpopts.IgnoreFields(Feature{}, "ParentSequence")); diff != "" {
		t.Errorf("Parsing the output of BuildGbk() does not produce the same output as parsing the original file read with ReadGbk(). Got this diff:\n%s", diff)
	}
}

func TestGbkLocationStringBuilder(t *testing.T) {
	tmpDataDir, err := ioutil.TempDir("", "data-*")
	defer os.RemoveAll(tmpDataDir)
	if err != nil {
		t.Error(err)
	}

	scrubbedGbk := ReadGbk("data/sample.gbk")

	// removing gbkLocationString from features to allow testing for gbkLocationBuilder
	for featureIndex := range scrubbedGbk.Features {
		scrubbedGbk.Features[featureIndex].GbkLocationString = ""
	}

	tmpGbkFilePath := filepath.Join(tmpDataDir, "sample.gbk")
	WriteGbk(scrubbedGbk, tmpGbkFilePath)

	testInputGbk := ReadGbk("data/sample.gbk")
	testOutputGbk := ReadGbk(tmpGbkFilePath)

	if diff := cmp.Diff(testInputGbk, testOutputGbk, cmpopts.IgnoreFields(Feature{}, "ParentSequence")); diff != "" {
		t.Errorf("Issue with partial location building. Parsing the output of BuildGbk() does not produce the same output as parsing the original file read with ReadGbk(). Got this diff:\n%s", diff)
	}
}

func TestGbLocationStringBuilder(t *testing.T) {
	tmpDataDir, err := ioutil.TempDir("", "data-*")
	defer os.RemoveAll(tmpDataDir)
	if err != nil {
		t.Error(err)
	}
	scrubbedGb := ReadGbk("data/t4_intron.gb")

	// removing gbkLocationString from features to allow testing for gbkLocationBuilder
	for featureIndex := range scrubbedGb.Features {
		scrubbedGb.Features[featureIndex].GbkLocationString = ""
	}

	tmpGbFilePath := filepath.Join(tmpDataDir, "t4_intron_test.gb")
	WriteGbk(scrubbedGb, tmpGbFilePath)

	testInputGb := ReadGbk("data/t4_intron.gb")
	testOutputGb := ReadGbk(tmpGbFilePath)

	if diff := cmp.Diff(testInputGb, testOutputGb, cmpopts.IgnoreFields(Feature{}, "ParentSequence")); diff != "" {
		t.Errorf("Issue with either Join or complement location building. Parsing the output of BuildGbk() does not produce the same output as parsing the original file read with ReadGbk(). Got this diff:\n%s", diff)
	}
}

func TestPartialLocationParseRegression(t *testing.T) {
	gbk := ReadGbk("data/sample.gbk")

	for _, feature := range gbk.Features {
		if feature.GbkLocationString == "687..3158>" && (feature.SequenceLocation.Start != 686 || feature.SequenceLocation.End != 3158) {
			t.Errorf("Partial location for three prime location parsing has failed. Parsing the output of BuildGbk() does not produce the same output as parsing the original file read with ReadGbk()")
		} else if feature.GbkLocationString == "<1..206" && (feature.SequenceLocation.Start != 0 || feature.SequenceLocation.End != 206) {
			t.Errorf("Partial location for five prime location parsing has failed. Parsing the output of BuildGbk() does not produce the same output as parsing the original file read with ReadGbk().")
		}
	}
}

func TestLocusParseRegression(t *testing.T) {
	gbk := ReadGbk("data/puc19.gbk").Meta.Locus
	json := ReadJSON("data/puc19static.json").Meta.Locus

	if diff := cmp.Diff(gbk, json, cmpopts.IgnoreFields(Feature{}, "ParentSequence")); diff != "" {
		t.Errorf("The meta parser has changed behaviour. Got this diff:\n%s", diff)
	}
}

func TestSnapgeneGenbankRegression(t *testing.T) {
	snapgene := ReadGbk("data/puc19_snapgene.gb")

	if snapgene.Sequence == "" {
		t.Errorf("Parsing snapgene returned an empty string")
	}
}

func TestGenbankNewlineParsingRegression(t *testing.T) {
	gbk := ReadGbk("data/puc19.gbk")

	for _, feature := range gbk.Features {
		if feature.SequenceLocation.Start == 410 && feature.SequenceLocation.End == 1750 && feature.Type == "CDS" {
			if feature.Attributes["product"] != "chromosomal replication initiator informational ATPase" {
				t.Errorf("Newline parsing has failed.")
			}
			break
		}
	}
}

func BenchmarkReadGbk(b *testing.B) {
	for i := 0; i < b.N; i++ {
		ReadGbk("data/bsub.gbk")
	}
}

func BenchmarkReadGbk1(b *testing.B)     { BenchmarkReadGbk(b) }
func BenchmarkReadGbk10(b *testing.B)    { BenchmarkReadGbk(b) }
func BenchmarkReadGbk100(b *testing.B)   { BenchmarkReadGbk(b) }
func BenchmarkReadGbk1000(b *testing.B)  { BenchmarkReadGbk(b) }
func BenchmarkReadGbk10000(b *testing.B) { BenchmarkReadGbk(b) }

/******************************************************************************

Gbk/gb/genbank related benchmarks end here.

******************************************************************************/

/******************************************************************************

JSON related tests begin here.

******************************************************************************/

func ExampleReadJSON() {
	sequence := ReadJSON("data/sample.json")

	fmt.Println(sequence.Meta.Source)
	//output: Saccharomyces cerevisiae (baker's yeast)
}

func ExampleParseJSON() {
	file, _ := ioutil.ReadFile("data/sample.json")
	sequence := ParseJSON(file)

	fmt.Println(sequence.Meta.Source)
	//output: Saccharomyces cerevisiae (baker's yeast)
}

func ExampleWriteJSON() {
	tmpDataDir, err := ioutil.TempDir("", "data-*")
	defer os.RemoveAll(tmpDataDir)
	if err != nil {
		fmt.Println(err.Error())
	}

	sequence := ReadJSON("data/sample.json")

	tmpJSONFilePath := filepath.Join(tmpDataDir, "sample.json")
	WriteJSON(sequence, tmpJSONFilePath)

	testSequence := ReadJSON(tmpJSONFilePath)

	fmt.Println(testSequence.Meta.Source)
	//output: Saccharomyces cerevisiae (baker's yeast)
}

func TestGbkToJSON(t *testing.T) {
	tmpDataDir, err := ioutil.TempDir("", "data-*")
	defer os.RemoveAll(tmpDataDir)
	if err != nil {
		t.Error(err)
	}

	testSequence := ReadGbk("data/puc19.gbk")

	tmpJSONFilePath := filepath.Join(tmpDataDir, "puc19.json")
	WriteJSON(testSequence, tmpJSONFilePath)

	readTestSequence := ReadJSON(tmpJSONFilePath)

	if diff := cmp.Diff(testSequence, readTestSequence, cmpopts.IgnoreFields(Feature{}, "ParentSequence")); diff != "" {
		t.Errorf(" mismatch (-want +got):\n%s", diff)
	}
}

func TestGffToJSON(t *testing.T) {
	tmpDataDir, err := ioutil.TempDir("", "data-*")
	defer os.RemoveAll(tmpDataDir)
	if err != nil {
		t.Error(err)
	}

	gffTestSequence := ReadGff("data/ecoli-mg1655-short.gff")

	tmpJSONFilePath := filepath.Join(tmpDataDir, "ecoli-mg1655-short.json")
	WriteJSON(gffTestSequence, tmpJSONFilePath)

	gffReadTestSequence := ReadJSON(tmpJSONFilePath)

	if diff := cmp.Diff(gffTestSequence, gffReadTestSequence, cmpopts.IgnoreFields(Feature{}, "ParentSequence")); diff != "" {
		t.Errorf(" mismatch (-want +got):\n%s", diff)
	}

}

/******************************************************************************

JSON related tests end here.

******************************************************************************/

/******************************************************************************

FASTA related tests begin here.

******************************************************************************/

// ExampleReadFASTA shows basic usage for ReadFASTA
func ExampleReadFASTA() {
	sequence := ReadFASTA("data/base.fasta")
	fmt.Println(sequence.Features[0].Description)
	// Output: gi|5524211|gb|AAD44166.1| cytochrome b [Elephas maximus maximus]
}

func ExampleParseFASTA() {
	file, _ := ioutil.ReadFile("data/base.fasta")
	sequence := ParseFASTA(file)

	fmt.Println(sequence.Features[0].Description)
	// Output: gi|5524211|gb|AAD44166.1| cytochrome b [Elephas maximus maximus]
}

func ExampleBuildFASTA() {
	sequence := ReadFASTA("data/base.fasta") // get example data
	fasta := BuildFASTA(sequence)            // build a fasta byte array
	firstLine := string(bytes.Split(fasta, []byte("\n"))[0])

	fmt.Println(firstLine)
	// Output: >gi|5524211|gb|AAD44166.1| cytochrome b [Elephas maximus maximus]
}

func ExampleWriteFASTA() {
	tmpDataDir, err := ioutil.TempDir("", "data-*")
	defer os.RemoveAll(tmpDataDir)
	if err != nil {
		fmt.Println(err.Error())
	}

	sequence := ReadFASTA("data/base.fasta") // get example data

	tmpFASTAFilePath := filepath.Join(tmpDataDir, "base.fasta")
	WriteFASTA(sequence, tmpFASTAFilePath) // write it out again

	testSequence := ReadFASTA(tmpFASTAFilePath) // read it in again

	fmt.Println(testSequence.Features[0].Description)
	// Output: gi|5524211|gb|AAD44166.1| cytochrome b [Elephas maximus maximus]
}

func TestFASTAIO(t *testing.T) {
	tmpDataDir, err := ioutil.TempDir("", "data-*")
	defer os.RemoveAll(tmpDataDir)
	if err != nil {
		t.Error(err)
	}

	inputFilename := "data/base.fasta"
	tmpFASTAFilePath := filepath.Join(tmpDataDir, "test.fasta")

	// read FASTA file
	testSequence := ReadFASTA(inputFilename)

	// write FASTA file
	WriteFASTA(testSequence, tmpFASTAFilePath)

	// read back and diff
	readTestSequence := ReadFASTA(tmpFASTAFilePath)

	if diff := cmp.Diff(testSequence, readTestSequence, cmpopts.IgnoreFields(Feature{}, "ParentSequence")); diff != "" {
		t.Errorf(" mismatch (-want +got):\n%s", diff)
	}
}

/******************************************************************************

FASTA related tests end here.

******************************************************************************/

/******************************************************************************

GbkMulti/GbkFlat related tests begin here.

******************************************************************************/

func ExampleReadGbkMulti() {
	sequences := ReadGbkMulti("data/multiGbk_test.seq")
	var locus []string
	for _, sequence := range sequences {
		locus = append(locus, sequence.Meta.Locus.Name)
	}

	fmt.Println(strings.Join(locus, ", "))
	// Output: AB000100, AB000106
}

func ExampleReadGbkFlat() {
	sequences := ReadGbkFlat("data/long_comment.seq")
	var locus []string
	for _, sequence := range sequences {
		locus = append(locus, sequence.Meta.Locus.Name)
	}

	fmt.Println(strings.Join(locus, ", "))
	// Output: AB000100, AB000106
}

func ExampleReadGbkFlatGz() {
	sequences := ReadGbkFlatGz("data/flatGbk_test.seq.gz")
	//sequences := ReadGbkFlatGz("data/gbbct358.seq.gz")
	var locus []string
	for _, sequence := range sequences {
		locus = append(locus, sequence.Meta.Locus.Name)
	}
	fmt.Println(strings.Join(locus, ", "))
	// Output: AB000100, AB000106
}

func ExampleParseGbkMulti() {
	file, _ := ioutil.ReadFile("data/multiGbk_test.seq")
	sequences := ParseGbkMulti(file)
	var locus []string
	for _, sequence := range sequences {
		locus = append(locus, sequence.Meta.Locus.Name)
	}

	fmt.Println(strings.Join(locus, ", "))
	// Output: AB000100, AB000106
}

func ExampleParseGbkFlat() {
	file, _ := ioutil.ReadFile("data/flatGbk_test.seq")
	sequences := ParseGbkFlat(file)
	var locus []string
	for _, sequence := range sequences {
		locus = append(locus, sequence.Meta.Locus.Name)
	}

	fmt.Println(strings.Join(locus, ", "))
	// Output: AB000100, AB000106
}

/******************************************************************************

GbkMulti/GbkFlat related tests end here.

******************************************************************************/
