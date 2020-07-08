package main

import (
	"io/ioutil"
	"os"
	"testing"

	"github.com/google/go-cmp/cmp"
	"github.com/pmezard/go-difflib/difflib"
)

/******************************************************************************

File is structured as so:

Gff - io tests, and benchmarks.
Gbk/gb/genbank - benchmarks.
JSON - io tests.

******************************************************************************/

/******************************************************************************

Gff related tests and benchmarks begin here.

******************************************************************************/

// TODO should delete output files.
func TestGffIO(t *testing.T) {
	testInputPath := "data/ecoli-mg1655.gff"
	testOutputPath := "data/test.gff"

	testSequence := ReadGff(testInputPath)
	WriteGff(testSequence, testOutputPath)

	readTestSequence := ReadGff(testOutputPath)

	if diff := cmp.Diff(testSequence, readTestSequence); diff != "" {
		t.Errorf("Parsing the output of BuildGff() does not produce the same output as parsing the original file read with ReadGff(). Got this diff:\n%s", diff)
	}

	original, _ := ioutil.ReadFile(testInputPath)
	builtOutput, _ := ioutil.ReadFile(testOutputPath)
	gffDiff := difflib.UnifiedDiff{
		A:        difflib.SplitLines(string(original)),
		B:        difflib.SplitLines(string(builtOutput)),
		FromFile: testInputPath,
		ToFile:   testOutputPath,
		Context:  3,
	}

	gffDiffText, _ := difflib.GetUnifiedDiffString(gffDiff)

	// cleaning up test data.
	os.Remove(testOutputPath)

	if gffDiffText != "" {
		t.Errorf("BuildGff() does not output the same file as was input through ReadGff(). Got this diff:\n%s", gffDiffText)
	}

}

func BenchmarkReadGff(b *testing.B) {
	for i := 0; i < b.N; i++ {
		ParseGff("data/ecoli-mg1655.gff")
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

func TestLocusParseRegression(t *testing.T) {
	gbk := ReadGbk("data/puc19.gbk").Meta.Locus
	json := ReadJSON("data/puc19static.json").Meta.Locus

	if diff := cmp.Diff(gbk, json); diff != "" {
		t.Errorf("The meta parser has changed behaviour. Got this diff:\n%s", diff)
	}
}

func TestSnapgeneGenbankRegression(t *testing.T) {
	snapgene := ReadGbk("data/puc19_snapgene.gb")

	if snapgene.Sequence.Sequence == "" {
		t.Errorf("Parsing snapgene returned an empty string")
	}
}

func TestGenbankNewlineParsingRegression(t *testing.T) {
	gbk := ReadGbk("data/bsub.gbk")

	for _, feature := range gbk.Features {
		if feature.Start == 410 && feature.End == 1750 && feature.Type == "CDS" {
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

// Does not work due to bug https://github.com/TimothyStiles/poly/issues/24
//func TestLocationParser(t *testing.T) {
//	gbk := ReadGbk("data/t4_intron.gb")
//	feature := gbk.Features[7]
//	s := "ATGAAACAATACCAAGATTTAATTAAAGACATTTTTGAAAATGGTTATGAAACCGATGATCGTACAGGCACAGGAACAATTGCTCTGTTCGGATCTAAATTACGCTGGGATTTAACTAAAGGTTTTCCTGCGGTAACAACTAAGAAGCTCGCCTGGAAAGCTTGCATTGCTGAGCTAATATGGTTTTTATCAGGAAGCACAAATGTCAATGATTTACGATTAATTCAACACGATTCGTTAATCCAAGGCAAAACAGTCTGGGATGAAAATTACGAAAATCAAGCAAAAGATTTAGGATACCATAGCGGTGAACTTGGTCCAATTTATGGAAAACAGTGGCGTGATTTTGGTGGTGTAGACCAAATTATAGAAGTTATTGATCGTATTAAAAAACTGCCAAATGATAGGCGTCAAATTGTTTCTGCATGGAATCCAGCTGAACTTAAATATATGGCATTACCGCCTTGTCATATGTTCTATCAGTTTAATGTGCGTAATGGCTATTTGGATTTGCAGTGGTATCAACGCTCAGTAGATGTTTTCTTGGGTCTACCGTTTAATATTGCGTCATATGCTACGTTAGTTCATATTGTAGCTAAGATGTGTAATCTTATTCCAGGGGATTTGATATTTTCTGGTGGTAATACTCATATCTATATGAATCACGTAGAACAATGTAAAGAAATTTTGAGGCGTGAACCTAAAGAGCTTTGTGAGCTGGTAATAAGTGGTCTACCTTATAAATTCCGATATCTTTCTACTAAAGAACAATTAAAATATGTTCTTAAACTTAGGCCTAAAGATTTCGTTCTTAACAACTATGTATCACACCCTCCTATTAAAGGAAAGATGGCGGTGTAA"
//	fmt.Println(gbk.Sequence.Sequence)
//	if featureSeq := feature.Sequence(gbk); featureSeq != s {
//		t.Errorf("Feature parser has changed behaviour. Got this:\n%s\n instead of \n%s", featureSeq, s)
//	}
//}

/******************************************************************************

Gbk/gb/genbank related benchmarks end here.

******************************************************************************/

/******************************************************************************

JSON related tests begin here.

******************************************************************************/

func TestJSONIO(t *testing.T) {
	testSequence := ReadGbk("data/bsub.gbk")
	WriteJSON(testSequence, "data/test.json")
	readTestSequence := ReadJSON("data/test.json")

	// cleaning up test data
	os.Remove("data/test.json")

	if diff := cmp.Diff(testSequence, readTestSequence); diff != "" {
		t.Errorf(" mismatch (-want +got):\n%s", diff)
	}
}

/******************************************************************************

JSON related tests end here.

******************************************************************************/
