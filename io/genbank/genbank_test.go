package genbank

import (
	"fmt"
	"io/ioutil"
	"os"
	"path/filepath"
	"strings"
	"testing"

	"github.com/TimothyStiles/poly"
	"github.com/google/go-cmp/cmp"
	"github.com/google/go-cmp/cmp/cmpopts"
)

/******************************************************************************

Gbk/gb/genbank related benchmarks begin here.

******************************************************************************/

func ExampleRead() {
	sequence := Read("../../data/puc19.gbk")
	fmt.Println(sequence.Meta.Locus.ModificationDate)
	// Output: 22-OCT-2019
}

func ExampleParse() {
	file, _ := ioutil.ReadFile("../../data/puc19.gbk")
	sequence := Parse(file)

	fmt.Println(sequence.Meta.Locus.ModificationDate)
	// Output: 22-OCT-2019
}

func ExampleBuild() {
	sequence := Read("../../data/puc19.gbk")
	gbkBytes := Build(sequence)
	testSequence := Parse(gbkBytes)

	fmt.Println(testSequence.Meta.Locus.ModificationDate)
	// Output: 22-OCT-2019
}

func ExampleWrite() {
	tmpDataDir, err := ioutil.TempDir("", "data-*")
	if err != nil {
		fmt.Println(err.Error())
	}
	defer os.RemoveAll(tmpDataDir)

	sequence := Read("../../data/puc19.gbk")

	tmpGbkFilePath := filepath.Join(tmpDataDir, "puc19.gbk")
	Write(sequence, tmpGbkFilePath)

	testSequence := Read(tmpGbkFilePath)

	fmt.Println(testSequence.Meta.Locus.ModificationDate)
	// Output: 22-OCT-2019
}

func TestGbkIO(t *testing.T) {
	tmpDataDir, err := ioutil.TempDir("", "data-*")
	if err != nil {
		t.Error(err)
	}
	defer os.RemoveAll(tmpDataDir)

	gbk := Read("../../data/puc19.gbk")

	tmpGbkFilePath := filepath.Join(tmpDataDir, "puc19.gbk")
	Write(gbk, tmpGbkFilePath)

	writeTestGbk := Read(tmpGbkFilePath)
	if diff := cmp.Diff(gbk, writeTestGbk, cmpopts.IgnoreFields(poly.Feature{}, "ParentSequence")); diff != "" {
		t.Errorf("Parsing the output of Build() does not produce the same output as parsing the original file read with Read(). Got this diff:\n%s", diff)
	}

	// Test multiline Genbank features
	pichia := Read("../../data/pichia_chr1_head.gb")
	var multilineOutput string
	for _, feature := range pichia.Features {
		multilineOutput = feature.GbkLocationString
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

	scrubbedGbk := Read("../../data/sample.gbk")

	// removing gbkLocationString from features to allow testing for gbkLocationBuilder
	for featureIndex := range scrubbedGbk.Features {
		scrubbedGbk.Features[featureIndex].GbkLocationString = ""
	}

	tmpGbkFilePath := filepath.Join(tmpDataDir, "sample.gbk")
	Write(scrubbedGbk, tmpGbkFilePath)

	testInputGbk := Read("../../data/sample.gbk")
	testOutputGbk := Read(tmpGbkFilePath)

	if diff := cmp.Diff(testInputGbk, testOutputGbk, cmpopts.IgnoreFields(poly.Feature{}, "ParentSequence")); diff != "" {
		t.Errorf("Issue with partial location building. Parsing the output of Build() does not produce the same output as parsing the original file read with Read(). Got this diff:\n%s", diff)
	}
}

func TestGbLocationStringBuilder(t *testing.T) {
	tmpDataDir, err := ioutil.TempDir("", "data-*")
	if err != nil {
		t.Error(err)
	}
	defer os.RemoveAll(tmpDataDir)

	scrubbedGb := Read("../../data/t4_intron.gb")

	// removing gbkLocationString from features to allow testing for gbkLocationBuilder
	for featureIndex := range scrubbedGb.Features {
		scrubbedGb.Features[featureIndex].GbkLocationString = ""
	}

	tmpGbFilePath := filepath.Join(tmpDataDir, "t4_intron_test.gb")
	Write(scrubbedGb, tmpGbFilePath)

	testInputGb := Read("../../data/t4_intron.gb")
	testOutputGb := Read(tmpGbFilePath)

	if diff := cmp.Diff(testInputGb, testOutputGb, cmpopts.IgnoreFields(poly.Feature{}, "ParentSequence")); diff != "" {
		t.Errorf("Issue with either Join or complement location building. Parsing the output of Build() does not produce the same output as parsing the original file read with Read(). Got this diff:\n%s", diff)
	}
}

func TestPartialLocationParseRegression(t *testing.T) {
	gbk := Read("../../data/sample.gbk")

	for _, feature := range gbk.Features {
		if feature.GbkLocationString == "687..3158>" && (feature.SequenceLocation.Start != 686 || feature.SequenceLocation.End != 3158) {
			t.Errorf("Partial location for three prime location parsing has failed. Parsing the output of Build() does not produce the same output as parsing the original file read with Read()")
		} else if feature.GbkLocationString == "<1..206" && (feature.SequenceLocation.Start != 0 || feature.SequenceLocation.End != 206) {
			t.Errorf("Partial location for five prime location parsing has failed. Parsing the output of Build() does not produce the same output as parsing the original file read with Read().")
		}
	}
}

func TestSnapgeneGenbankRegression(t *testing.T) {
	snapgene := Read("../../data/puc19_snapgene.gb")

	if snapgene.Sequence == "" {
		t.Errorf("Parsing snapgene returned an empty string")
	}
}

func TestGetSequenceMethod(t *testing.T) {

	gbk := Read("../../data/t4_intron.gb")

	// Check to see if GetSequence method works on Features struct
	feature := gbk.Features[1].GetSequence()
	seq := "atgagattacaacgccagagcatcaaagattcagaagttagaggtaaatggtattttaatatcatcggtaaagattctgaacttgttgaaaaagctgaacatcttttacgtgatatgggatgggaagatgaatgcgatggatgtcctctttatgaagacggagaaagcgcaggattttggatttaccattctgacgtcgagcagtttaaagctgattggaaaattgtgaaaaagtctgtttga"
	if feature != seq {
		t.Errorf("Feature GetSequence method has failed. Got this:\n%s instead of \n%s", feature, seq)
	}

}

func TestLocationParser(t *testing.T) {
	gbk := Read("../../data/t4_intron.gb")

	// Read 1..243
	feature := gbk.Features[1].GetSequence()
	seq := "atgagattacaacgccagagcatcaaagattcagaagttagaggtaaatggtattttaatatcatcggtaaagattctgaacttgttgaaaaagctgaacatcttttacgtgatatgggatgggaagatgaatgcgatggatgtcctctttatgaagacggagaaagcgcaggattttggatttaccattctgacgtcgagcagtttaaagctgattggaaaattgtgaaaaagtctgtttga"
	if feature != seq {
		t.Errorf("Feature sequence parser has changed on test '1..243'. Got this:\n%s instead of \n%s", feature, seq)
	}

	// Read join(893..1441,2459..2770)
	featureJoin := gbk.Features[6].GetSequence()
	seqJoin := "atgaaacaataccaagatttaattaaagacatttttgaaaatggttatgaaaccgatgatcgtacaggcacaggaacaattgctctgttcggatctaaattacgctgggatttaactaaaggttttcctgcggtaacaactaagaagctcgcctggaaagcttgcattgctgagctaatatggtttttatcaggaagcacaaatgtcaatgatttacgattaattcaacacgattcgttaatccaaggcaaaacagtctgggatgaaaattacgaaaatcaagcaaaagatttaggataccatagcggtgaacttggtccaatttatggaaaacagtggcgtgattttggtggtgtagaccaaattatagaagttattgatcgtattaaaaaactgccaaatgataggcgtcaaattgtttctgcatggaatccagctgaacttaaatatatggcattaccgccttgtcatatgttctatcagtttaatgtgcgtaatggctatttggatttgcagtggtatcaacgctcagtagatgttttcttgggtctaccgtttaatattgcgtcatatgctacgttagttcatattgtagctaagatgtgtaatcttattccaggggatttgatattttctggtggtaatactcatatctatatgaatcacgtagaacaatgtaaagaaattttgaggcgtgaacctaaagagctttgtgagctggtaataagtggtctaccttataaattccgatatctttctactaaagaacaattaaaatatgttcttaaacttaggcctaaagatttcgttcttaacaactatgtatcacaccctcctattaaaggaaagatggcggtgtaa"
	if featureJoin != seqJoin {
		t.Errorf("Feature sequence parser has changed on test 'join(893..1441,2459..2770)'. Got this:\n%s instead of \n%s", featureJoin, seqJoin)
	}

	// Read complement(2791..3054)
	featureComplement := gbk.Features[10].GetSequence()
	seqComplement := "ttattcactacccggcatagacggcccacgctggaataattcgtcatattgtttttccgttaaaacagtaatatcgtagtaacagtcagaagaagttttaactgtggaaattttattatcaaaatactcacgagtcattttatgagtatagtattttttaccataaatggtaataggctgttctggtcctggaacttctaactcgcttgggttaggaagtgtaaaaagaactacaccagaagtatctttaaatcgtaaaatcat"
	if featureComplement != seqComplement {
		t.Errorf("Feature sequence parser has changed on test 'complement(2791..3054)'. Got this:\n%s instead of \n%s", featureComplement, seqComplement)
	}

	// Read join(complement(315..330),complement(339..896))
	// Note: it is known that some software, like Snapgene, assumes that since both strands are in the reverse direction
	// that the first sequence should be appended to the reverse sequence, instead of the second sequence
	// getting appended to the first. Biopython appends the second sequence to the first, and that is logically
	// the most obvious thing to do, so we are implementing it that way.
	featureJoinComplement := gbk.Features[3].GetSequence()
	seqJoinComplement := "ataccaatttaatcattcatttatatactgattccgtaagggttgttacttcatctattttataccaatgcgtttcaaccatttcacgcttgcttatatcatcaagaaaacttgcgtctaattgaactgttgaattaacacgatgccttttaacgatgcgagaaacaactacttcatctgcataaggtaatgcagcatataacagagcaggcccgccaattacacttactttagaattctgatcaagcatagtttcgaatggtgcattagggcttgacacttgaatttcgccgccagaaatgtaagttatatattgctcccaagtaatatagaaatgtgctaaatcgccgtctttagttacaggataatcacgcgcaaggtcacacaccacaatatggctacgaccaggaagtaatgtaggcaatgactggaacgttttagcacccataatcataattgtgccttcagtacgagctttaaaattctggaggtcctttttaactcgtccccatggtaaaccatcacctaaaccgaatgctaattcattaaagccgtcgaccgttttagttggaga"
	if featureJoinComplement != seqJoinComplement {
		t.Errorf("Feature sequence parser has changed on test 'join(complement(315..330),complement(339..896))'. Got this:\n%s instead of \n%s", featureJoinComplement, seqJoinComplement)
	}

	// Read complement(join(893..1098,1101..2770))
	featureComplementJoin := gbk.Features[5].GetSequence()
	seqComplementJoin := "ttacaccgccatctttcctttaataggagggtgtgatacatagttgttaagaacgaaatctttaggcctaagtttaagaacatattttaattgttctttagtagaaagatatcggaatttataaggtagaccacttattaccagctcacaaagctctttaggttcacgcctcaaaatttctttacattgttctacgtgattcatatagatatgagtattaccaccagaaaatatcaaatcccctggaataagattacacatcttagctacaatatgaactaacgtagcatatgacgcaatattaaacggtagcattatgttcagataaggtcgttaatcttaccccggaattatatccagctgcatgtcaccatgcagagcagactatatctccaacttgttaaagcaagttgtctatcgtttcgagtcacttgaccctactccccaaagggatagtcgttaggcatttatgtagaaccaattccatttatcagattttacacgataagtaactaatccagacgaaattttaaaatgtctagctgcatctgctgcacaatcaaaaataaccccatcacatgaaatctttttaatattactaggctttttacctttcatcttttctgatattttagatttagttatgtctgaatgcttatgattaaagaatgaattattttcacctgaacgatttctgcatttactacaagtataagcagaagtttgtatgcgaacaccgcacttacaaaacttatgggtttctggattccaacgcccgtttttacttccgggtttactgtaaagagctttccgaccatcaggtccaagtttaagcatcttagctttaacagtttcagaacgtttcttaataatttcttcttttaatggatgcgtagaacatgtatcaccaaacgttgcatcagcaatattgtatccattaattttagaattaagctctttaatccaaaaattttctcgttcaataatcaaatctttctcatatggaatttcttccaaaatagaacattcaaacacattaccatgtttgttaaaagacctctgaagttttatagaagaatggcatcctttttctaaatctttaaaatgcctcttccatctcttttcaaaatctttagcacttcctacatatactttattgtttaaagtatttttaatctgataaattccgcttttcataaatacctctttaaatatagaagtatttattaaagggcaagtcctacaatttagcacgggattgtctactagagaggttccccgtttagatagattacaagtataagtcaccttatactcaggcctcaattaacccaagaaaacatctactgagcgttgataccactgcaaatccaaatagccattacgcacattaaactgatagaacatatgacaaggcggtaatgccatatatttaagttcagctggattccatgcagaaacaatttgacgcctatcatttggcagttttttaatacgatcaataacttctataatttggtctacaccaccaaaatcacgccactgttttccataaattggaccaagttcaccgctatggtatcctaaatcttttgcttgattttcgtaattttcatcccagactgttttgccttggattaacgaatcgtgttgaattaatcgtaaatcatacatttgtgcttcctgataaaaaccatattagctcagcaatgcaagctttccaggcgagcttcttagttgttaccgcaggaaaacctttagttaaatcccagcgtaatttagatccgaacagagcaattgttcctgtgcctgtacgatcatcggtttcataaccattttcaaaaatgtctttaattaaatcttggtattgtttcat"
	if featureComplementJoin != seqComplementJoin {
		t.Errorf("Feature sequence parser has changed on test 'complement(join(893..1098,1101..2770))'. Got this:\n%s instead of \n%s", featureComplementJoin, seqComplementJoin)
	}
}

func TestGenbankNewlineParsingRegression(t *testing.T) {
	gbk := Read("../../data/puc19.gbk")

	for _, feature := range gbk.Features {
		if feature.SequenceLocation.Start == 410 && feature.SequenceLocation.End == 1750 && feature.Type == "CDS" {
			if feature.Attributes["product"] != "chromosomal replication initiator informational ATPase" {
				t.Errorf("Newline parsing has failed.")
			}
			break
		}
	}
}

func BenchmarkRead(b *testing.B) {
	for i := 0; i < b.N; i++ {
		Read("../../data/bsub.gbk")
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

/******************************************************************************

GbkMulti/GbkFlat related tests begin here.

******************************************************************************/

func ExampleReadMulti() {
	sequences := ReadMulti("../../data/multiGbk_test.seq")
	var locus []string
	for _, sequence := range sequences {
		locus = append(locus, sequence.Meta.Locus.Name)
	}

	fmt.Println(strings.Join(locus, ", "))
	// Output: AB000100, AB000106
}

func ExampleReadFlat() {
	sequences := ReadFlat("../../data/long_comment.seq")
	var locus []string
	for _, sequence := range sequences {
		locus = append(locus, sequence.Meta.Locus.Name)
	}

	fmt.Println(strings.Join(locus, ", "))
	// Output: AB000100, AB000106
}

func ExampleReadFlatGz() {
	sequences := ReadFlatGz("../../data/flatGbk_test.seq.gz")
	//sequences := ReadFlatGz("../../data/gbbct358.seq.gz")
	var locus []string
	for _, sequence := range sequences {
		locus = append(locus, sequence.Meta.Locus.Name)
	}
	fmt.Println(strings.Join(locus, ", "))
	// Output: AB000100, AB000106
}

func ExampleParseMulti() {
	file, _ := ioutil.ReadFile("../../data/multiGbk_test.seq")
	sequences := ParseMulti(file)
	var locus []string
	for _, sequence := range sequences {
		locus = append(locus, sequence.Meta.Locus.Name)
	}

	fmt.Println(strings.Join(locus, ", "))
	// Output: AB000100, AB000106
}

func ExampleParseFlat() {
	file, _ := ioutil.ReadFile("../../data/flatGbk_test.seq")
	sequences := ParseFlat(file)
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
