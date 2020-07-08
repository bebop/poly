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

func TestLocationParser(t *testing.T) {
	gbk := ReadGbk("data/t4_intron.gb")

	// Read 1..243
	feature := gbk.Features[1].Sequence(gbk)
	seq := "atgagattacaacgccagagcatcaaagattcagaagttagaggtaaatggtattttaatatcatcggtaaagattctgaacttgttgaaaaagctgaacatcttttacgtgatatgggatgggaagatgaatgcgatggatgtcctctttatgaagacggagaaagcgcaggattttggatttaccattctgacgtcgagcagtttaaagctgattggaaaattgtgaaaaagtctgtttga"
	if feature != seq {
		t.Errorf("Feature sequence parser has changed on test '1..243'. Got this:\n%s instead of \n%s", feature, seq)
	}

	// Read join(893..1441,2459..2770)
	featureJoin := gbk.Features[6].Sequence(gbk)
	seqJoin := "atgaaacaataccaagatttaattaaagacatttttgaaaatggttatgaaaccgatgatcgtacaggcacaggaacaattgctctgttcggatctaaattacgctgggatttaactaaaggttttcctgcggtaacaactaagaagctcgcctggaaagcttgcattgctgagctaatatggtttttatcaggaagcacaaatgtcaatgatttacgattaattcaacacgattcgttaatccaaggcaaaacagtctgggatgaaaattacgaaaatcaagcaaaagatttaggataccatagcggtgaacttggtccaatttatggaaaacagtggcgtgattttggtggtgtagaccaaattatagaagttattgatcgtattaaaaaactgccaaatgataggcgtcaaattgtttctgcatggaatccagctgaacttaaatatatggcattaccgccttgtcatatgttctatcagtttaatgtgcgtaatggctatttggatttgcagtggtatcaacgctcagtagatgttttcttgggtctaccgtttaatattgcgtcatatgctacgttagttcatattgtagctaagatgtgtaatcttattccaggggatttgatattttctggtggtaatactcatatctatatgaatcacgtagaacaatgtaaagaaattttgaggcgtgaacctaaagagctttgtgagctggtaataagtggtctaccttataaattccgatatctttctactaaagaacaattaaaatatgttcttaaacttaggcctaaagatttcgttcttaacaactatgtatcacaccctcctattaaaggaaagatggcggtgtaa"
	if featureJoin != seqJoin {
		t.Errorf("Feature sequence parser has changed on test 'join(893..1441,2459..2770)'. Got this:\n%s instead of \n%s", featureJoin, seqJoin)
	}

	// Read complement(2791..3054)
	featureComplement := gbk.Features[10].Sequence(gbk)
	seqComplement := "ttattcactacccggcatagacggcccacgctggaataattcgtcatattgtttttccgttaaaacagtaatatcgtagtaacagtcagaagaagttttaactgtggaaattttattatcaaaatactcacgagtcattttatgagtatagtattttttaccataaatggtaataggctgttctggtcctggaacttctaactcgcttgggttaggaagtgtaaaaagaactacaccagaagtatctttaaatcgtaaaatcat"
	if featureComplement != seqComplement {
		t.Errorf("Feature sequence parser has changed on test 'complement(2791..3054)'. Got this:\n%s instead of \n%s", featureComplement, seqComplement)
	}

	// Read join(complement(315..330),complement(339..896))
	// Note: it is known that some software, like Snapgene, assumes that since both strands are in the reverse direction
	// that the first sequence should be appended to the reverse sequence, instead of the second sequence
	// getting appended to the first. Biopython appends the second sequence to the first, and that is logically
	// the most obvious thing to do, so we are implementing it that way.
	featureJoinComplement := gbk.Features[3].Sequence(gbk)
	seqJoinComplement := "ataccaatttaatcattcatttatatactgattccgtaagggttgttacttcatctattttataccaatgcgtttcaaccatttcacgcttgcttatatcatcaagaaaacttgcgtctaattgaactgttgaattaacacgatgccttttaacgatgcgagaaacaactacttcatctgcataaggtaatgcagcatataacagagcaggcccgccaattacacttactttagaattctgatcaagcatagtttcgaatggtgcattagggcttgacacttgaatttcgccgccagaaatgtaagttatatattgctcccaagtaatatagaaatgtgctaaatcgccgtctttagttacaggataatcacgcgcaaggtcacacaccacaatatggctacgaccaggaagtaatgtaggcaatgactggaacgttttagcacccataatcataattgtgccttcagtacgagctttaaaattctggaggtcctttttaactcgtccccatggtaaaccatcacctaaaccgaatgctaattcattaaagccgtcgaccgttttagttggaga"
	if featureJoinComplement != seqJoinComplement {
		t.Errorf("Feature sequence parser has changed on test 'join(complement(315..330),complement(339..896))'. Got this:\n%s instead of \n%s", featureJoinComplement, seqJoinComplement)
	}

	// Read complement(join(893..1098,1101..2770))
	featureComplementJoin := gbk.Features[5].Sequence(gbk)
	seqComplementJoin := "ttacaccgccatctttcctttaataggagggtgtgatacatagttgttaagaacgaaatctttaggcctaagtttaagaacatattttaattgttctttagtagaaagatatcggaatttataaggtagaccacttattaccagctcacaaagctctttaggttcacgcctcaaaatttctttacattgttctacgtgattcatatagatatgagtattaccaccagaaaatatcaaatcccctggaataagattacacatcttagctacaatatgaactaacgtagcatatgacgcaatattaaacggtagcattatgttcagataaggtcgttaatcttaccccggaattatatccagctgcatgtcaccatgcagagcagactatatctccaacttgttaaagcaagttgtctatcgtttcgagtcacttgaccctactccccaaagggatagtcgttaggcatttatgtagaaccaattccatttatcagattttacacgataagtaactaatccagacgaaattttaaaatgtctagctgcatctgctgcacaatcaaaaataaccccatcacatgaaatctttttaatattactaggctttttacctttcatcttttctgatattttagatttagttatgtctgaatgcttatgattaaagaatgaattattttcacctgaacgatttctgcatttactacaagtataagcagaagtttgtatgcgaacaccgcacttacaaaacttatgggtttctggattccaacgcccgtttttacttccgggtttactgtaaagagctttccgaccatcaggtccaagtttaagcatcttagctttaacagtttcagaacgtttcttaataatttcttcttttaatggatgcgtagaacatgtatcaccaaacgttgcatcagcaatattgtatccattaattttagaattaagctctttaatccaaaaattttctcgttcaataatcaaatctttctcatatggaatttcttccaaaatagaacattcaaacacattaccatgtttgttaaaagacctctgaagttttatagaagaatggcatcctttttctaaatctttaaaatgcctcttccatctcttttcaaaatctttagcacttcctacatatactttattgtttaaagtatttttaatctgataaattccgcttttcataaatacctctttaaatatagaagtatttattaaagggcaagtcctacaatttagcacgggattgtctactagagaggttccccgtttagatagattacaagtataagtcaccttatactcaggcctcaattaacccaagaaaacatctactgagcgttgataccactgcaaatccaaatagccattacgcacattaaactgatagaacatatgacaaggcggtaatgccatatatttaagttcagctggattccatgcagaaacaatttgacgcctatcatttggcagttttttaatacgatcaataacttctataatttggtctacaccaccaaaatcacgccactgttttccataaattggaccaagttcaccgctatggtatcctaaatcttttgcttgattttcgtaattttcatcccagactgttttgccttggattaacgaatcgtgttgaattaatcgtaaatcatacatttgtgcttcctgataaaaaccatattagctcagcaatgcaagctttccaggcgagcttcttagttgttaccgcaggaaaacctttagttaaatcccagcgtaatttagatccgaacagagcaattgttcctgtgcctgtacgatcatcggtttcataaccattttcaaaaatgtctttaattaaatcttggtattgtttcat"
	if featureComplementJoin != seqComplementJoin {
		t.Errorf("Feature sequence parser has changed on test 'complement(join(893..1098,1101..2770))'. Got this:\n%s instead of \n%s", featureComplementJoin, seqComplementJoin)
	}
}

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
