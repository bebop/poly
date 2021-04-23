package poly

import (
	"fmt"
	"io/ioutil"
	"os"
	"strings"
	"sync"
	"testing"

	"github.com/google/go-cmp/cmp"
)

func ExampleTranslate() {

	gfpTranslation := "MASKGEELFTGVVPILVELDGDVNGHKFSVSGEGEGDATYGKLTLKFICTTGKLPVPWPTLVTTFSYGVQCFSRYPDHMKRHDFFKSAMPEGYVQERTISFKDDGNYKTRAEVKFEGDTLVNRIELKGIDFKEDGNILGHKLEYNYNSHNVYITADKQKNGIKANFKIRHNIEDGSVQLADHYQQNTPIGDGPVLLPDNHYLSTQSALSKDPNEKRDHMVLLEFVTAAGITHGMDELYK*"
	gfpDnaSequence := "ATGGCTAGCAAAGGAGAAGAACTTTTCACTGGAGTTGTCCCAATTCTTGTTGAATTAGATGGTGATGTTAATGGGCACAAATTTTCTGTCAGTGGAGAGGGTGAAGGTGATGCTACATACGGAAAGCTTACCCTTAAATTTATTTGCACTACTGGAAAACTACCTGTTCCATGGCCAACACTTGTCACTACTTTCTCTTATGGTGTTCAATGCTTTTCCCGTTATCCGGATCATATGAAACGGCATGACTTTTTCAAGAGTGCCATGCCCGAAGGTTATGTACAGGAACGCACTATATCTTTCAAAGATGACGGGAACTACAAGACGCGTGCTGAAGTCAAGTTTGAAGGTGATACCCTTGTTAATCGTATCGAGTTAAAAGGTATTGATTTTAAAGAAGATGGAAACATTCTCGGACACAAACTCGAGTACAACTATAACTCACACAATGTATACATCACGGCAGACAAACAAAAGAATGGAATCAAAGCTAACTTCAAAATTCGCCACAACATTGAAGATGGATCCGTTCAACTAGCAGACCATTATCAACAAAATACTCCAATTGGCGATGGCCCTGTCCTTTTACCAGACAACCATTACCTGTCGACACAATCTGCCCTTTCGAAAGATCCCAACGAAAAGCGTGACCACATGGTCCTTCTTGAGTTTGTAACTGCTGCTGGGATTACACATGGCATGGATGAGCTCTACAAATAA"
	testTranslation := Translate(gfpDnaSequence, GetCodonTable(11)) // need to specify which codons map to which amino acids per NCBI table

	fmt.Println(gfpTranslation == testTranslation)
	// output: true
}

func TestTranslation(t *testing.T) {
	gfpTranslation := "MASKGEELFTGVVPILVELDGDVNGHKFSVSGEGEGDATYGKLTLKFICTTGKLPVPWPTLVTTFSYGVQCFSRYPDHMKRHDFFKSAMPEGYVQERTISFKDDGNYKTRAEVKFEGDTLVNRIELKGIDFKEDGNILGHKLEYNYNSHNVYITADKQKNGIKANFKIRHNIEDGSVQLADHYQQNTPIGDGPVLLPDNHYLSTQSALSKDPNEKRDHMVLLEFVTAAGITHGMDELYK*"
	gfpDnaSequence := "ATGGCTAGCAAAGGAGAAGAACTTTTCACTGGAGTTGTCCCAATTCTTGTTGAATTAGATGGTGATGTTAATGGGCACAAATTTTCTGTCAGTGGAGAGGGTGAAGGTGATGCTACATACGGAAAGCTTACCCTTAAATTTATTTGCACTACTGGAAAACTACCTGTTCCATGGCCAACACTTGTCACTACTTTCTCTTATGGTGTTCAATGCTTTTCCCGTTATCCGGATCATATGAAACGGCATGACTTTTTCAAGAGTGCCATGCCCGAAGGTTATGTACAGGAACGCACTATATCTTTCAAAGATGACGGGAACTACAAGACGCGTGCTGAAGTCAAGTTTGAAGGTGATACCCTTGTTAATCGTATCGAGTTAAAAGGTATTGATTTTAAAGAAGATGGAAACATTCTCGGACACAAACTCGAGTACAACTATAACTCACACAATGTATACATCACGGCAGACAAACAAAAGAATGGAATCAAAGCTAACTTCAAAATTCGCCACAACATTGAAGATGGATCCGTTCAACTAGCAGACCATTATCAACAAAATACTCCAATTGGCGATGGCCCTGTCCTTTTACCAGACAACCATTACCTGTCGACACAATCTGCCCTTTCGAAAGATCCCAACGAAAAGCGTGACCACATGGTCCTTCTTGAGTTTGTAACTGCTGCTGGGATTACACATGGCATGGATGAGCTCTACAAATAA"
	if got := Translate(gfpDnaSequence, GetCodonTable(11)); got != gfpTranslation {
		t.Errorf("TestTranslation has failed. Translate has returned %q, want %q", got, gfpTranslation)
	}

}

func ExampleOptimize() {

	gfpTranslation := "MASKGEELFTGVVPILVELDGDVNGHKFSVSGEGEGDATYGKLTLKFICTTGKLPVPWPTLVTTFSYGVQCFSRYPDHMKRHDFFKSAMPEGYVQERTISFKDDGNYKTRAEVKFEGDTLVNRIELKGIDFKEDGNILGHKLEYNYNSHNVYITADKQKNGIKANFKIRHNIEDGSVQLADHYQQNTPIGDGPVLLPDNHYLSTQSALSKDPNEKRDHMVLLEFVTAAGITHGMDELYK*"

	sequence := ReadGbk("data/puc19.gbk")
	codonTable := GetCodonTable(11)

	optimizationTable := sequence.GetOptimizationTable(codonTable)

	optimizedSequence := Optimize(gfpTranslation, optimizationTable)
	optimizedSequenceTranslation := Translate(optimizedSequence, optimizationTable)

	fmt.Println(optimizedSequenceTranslation == gfpTranslation)
	// output: true
}

func TestOptimize(t *testing.T) {
	gfpTranslation := "MASKGEELFTGVVPILVELDGDVNGHKFSVSGEGEGDATYGKLTLKFICTTGKLPVPWPTLVTTFSYGVQCFSRYPDHMKRHDFFKSAMPEGYVQERTISFKDDGNYKTRAEVKFEGDTLVNRIELKGIDFKEDGNILGHKLEYNYNSHNVYITADKQKNGIKANFKIRHNIEDGSVQLADHYQQNTPIGDGPVLLPDNHYLSTQSALSKDPNEKRDHMVLLEFVTAAGITHGMDELYK*"

	sequence := ReadGbk("data/puc19.gbk")
	codonTable := GetCodonTable(11)

	optimizationTable := sequence.GetOptimizationTable(codonTable)

	optimizedSequence := Optimize(gfpTranslation, optimizationTable)
	optimizedSequenceTranslation := Translate(optimizedSequence, optimizationTable)

	if optimizedSequenceTranslation != gfpTranslation {
		t.Errorf("TestOptimize has failed. Translate has returned %q, want %q", optimizedSequenceTranslation, gfpTranslation)
	}

}

func TestGetCodonFrequency(t *testing.T) {

	translationTable := GetCodonTable(11).generateTranslationTable()

	var codons strings.Builder

	for codon := range translationTable {
		codons.WriteString(codon)
	}

	// converting to string as saved variable for easier debugging.
	codonString := codons.String()

	// getting length as string for easier debugging.
	codonStringlength := len(codonString)

	if codonStringlength != (64 * 3) {
		t.Errorf("TestGetCodonFrequency has failed. aggregrated codon string is not the correct length.")
	}

	codonFrequencyHashMap := getCodonFrequency(codonString)

	if len(codonFrequencyHashMap) != 64 {
		t.Errorf("TestGetCodonFrequency has failed. codonFrequencyHashMap does not contain every codon.")
	}

	for codon, frequency := range codonFrequencyHashMap {
		if frequency != 1 {
			t.Errorf("TestGetCodonFrequency has failed. The codon \"%q\" appears %q times when it should only appear once.", codon, frequency)
		}
	}

	doubleCodonFrequencyHashMap := getCodonFrequency(codonString + codonString)

	if len(doubleCodonFrequencyHashMap) != 64 {
		t.Errorf("TestGetCodonFrequency has failed. doubleCodonFrequencyHashMap does not contain every codon.")
	}

	for codon, frequency := range doubleCodonFrequencyHashMap {
		if frequency != 2 {
			t.Errorf("TestGetCodonFrequency has failed. The codon \"%q\" appears %q times when it should only appear twice.", codon, frequency)
		}
	}

}

/******************************************************************************

JSON related tests begin here.

******************************************************************************/

func ExampleReadCodonJSON() {
	codontable := ReadCodonJSON("data/bsub_codon_test.json")

	fmt.Println(codontable.AminoAcids[0].Codons[0].Weight)
	//output: 28327
}

func ExampleParseCodonJSON() {
	file, _ := ioutil.ReadFile("data/bsub_codon_test.json")
	codontable := ParseCodonJSON(file)

	fmt.Println(codontable.AminoAcids[0].Codons[0].Weight)
	//output: 28327
}

func ExampleWriteCodonJSON() {
	codontable := ReadCodonJSON("data/bsub_codon_test.json")
	WriteCodonJSON(codontable, "data/codon_test.json")
	testCodonTable := ReadCodonJSON("data/codon_test.json")

	// cleaning up test data
	os.Remove("data/codon_test.json")

	fmt.Println(testCodonTable.AminoAcids[0].Codons[0].Weight)
	//output: 28327
}

func TestWriteCodonJSON(t *testing.T) {
	testCodonTable := ReadCodonJSON("data/bsub_codon_test.json")
	WriteCodonJSON(testCodonTable, "data/codon_test1.json")
	readTestCodonTable := ReadCodonJSON("data/codon_test1.json")

	// cleaning up test data
	os.Remove("data/codon_test1.json")

	if diff := cmp.Diff(testCodonTable, readTestCodonTable); diff != "" {
		t.Errorf(" mismatch (-want +got):\n%s", diff)
	}

}

/******************************************************************************

Synthesis Fixer tests begin here

Generally, we are taking tests previously used in DnaChisel along with a few of
our own tests.

https://github.com/Edinburgh-Genome-Foundry/DnaChisel/tree/master/tests/builtin_specifications

The following tests are built into DnaChisel:

- test_AllowPrimer
- test_AvoidBlastMatches
- test_AvoidChanges
- test_AvoidHairpins
- test_AvoidMatches
- test_AvoidRareCodons
- test_AvoidStopCodons
- test_CodonOptimize
- test_EnforceChanges
- test_EnforceChoice
- test_EnforceGCContent
- test_EnforcePatternOccurence
- test_EnforceRegionsCompatibility
- test_EnforceSequence
- test_EnforceTerminalGCContent
- test_EnforceTranslation
- test_SequenceLengthBounds
- test_UniquifyAllKmers




******************************************************************************/

func ExampleFixCds() {
	bla := "ATGAGTATTCAACATTTCCGTGTCGCCCTTATTCCCTTTTTTGCGGCATTTTGCCTTCCTGTTTTTGCTCACCCAGAAACGCTGGTGAAAGTAAAAGATGCTGAAGATCAGTTGGGTGCACGAGTGGGTTACATCGAACTGGATCTCAACAGCGGTAAGATCCTTGAGAGTTTTCGCCCCGAAGAACGTTTTCCAATGATGAGCACTTTTAAAGTTCTGCTATGTGGCGCGGTATTATCCCGTATTGACGCCGGGCAAGAGCAACTCGGTCGCCGCATACACTATTCTCAGAATGACTTGGTTGAGTACTCACCAGTCACAGAAAAGCATCTTACGGATGGCATGACAGTAAGAGAATTATGCAGTGCTGCCATAACCATGAGTGATAACACTGCGGCCAACTTACTTCTGACAACGATCGGAGGACCGAAGGAGCTAACCGCTTTTTTGCACAACATGGGGGATCATGTAACTCGCCTTGATCGTTGGGAACCGGAGCTGAATGAAGCCATACCAAACGACGAGCGTGACACCACGATGCCTGTAGCAATGGCAACAACGTTGCGCAAACTATTAACTGGCGAACTACTTACTCTAGCTTCCCGGCAACAATTAATAGACTGGATGGAGGCGGATAAAGTTGCAGGACCACTTCTGCGCTCGGCCCTTCCGGCTGGCTGGTTTATTGCTGATAAATCTGGAGCCGGTGAGCGTGGGTCTCGCGGTATCATTGCAGCACTGGGGCCAGATGGTAAGCCCTCCCGTATCGTAGTTATCTACACGACGGGGAGTCAGGCAACTATGGATGAACGAAATAGACAGATCGCTGAGATAGGTGCCTCACTGATTAAGCATTGGTAA"
	sequence := ReadGbk("data/ecoli-mg1655.gff")
	codonTable := GetCodonTable(11)
	optimizationTable := sequence.GetOptimizationTable(codonTable)

	var functions []func(string, chan DnaSuggestion, *sync.WaitGroup)
	functions = append(functions, FindBsaI)

	fixedSeq, _ := FixCds(bla, optimizationTable, functions)
	fmt.Println(fixedSeq)

	// Output: ATGAGTATTCAACATTTCCGTGTCGCCCTTATTCCCTTTTTTGCGGCATTTTGCCTTCCTGTTTTTGCTCACCCAGAAACGCTGGTGAAAGTAAAAGATGCTGAAGATCAGTTGGGTGCACGAGTGGGTTACATCGAACTGGATCTCAACAGCGGTAAGATCCTTGAGAGTTTTCGCCCCGAAGAACGTTTTCCAATGATGAGCACTTTTAAAGTTCTGCTATGTGGCGCGGTATTATCCCGTATTGACGCCGGGCAAGAGCAACTCGGTCGCCGCATACACTATTCTCAGAATGACTTGGTTGAGTACTCACCAGTCACAGAAAAGCATCTTACGGATGGCATGACAGTAAGAGAATTATGCAGTGCTGCCATAACCATGAGTGATAACACTGCGGCCAACTTACTTCTGACAACGATCGGAGGACCGAAGGAGCTAACCGCTTTTTTGCACAACATGGGGGATCATGTAACTCGCCTTGATCGTTGGGAACCGGAGCTGAATGAAGCCATACCAAACGACGAGCGTGACACCACGATGCCTGTAGCAATGGCAACAACGTTGCGCAAACTATTAACTGGCGAACTACTTACTCTAGCTTCCCGGCAACAATTAATAGACTGGATGGAGGCGGATAAAGTTGCAGGACCACTTCTGCGCTCGGCCCTTCCGGCTGGCTGGTTTATTGCTGATAAATCTGGAGCCGGTGAGCGTGGCTCTCGCGGTATCATTGCAGCACTGGGGCCAGATGGTAAGCCCTCCCGTATCGTAGTTATCTACACGACGGGGAGTCAGGCAACTATGGATGAACGAAATAGACAGATCGCTGAGATAGGTGCCTCACTGATTAAGCATTGG

}

/******************************************************************************

Codon Compromise + Add related tests begin here.

******************************************************************************/

func ExampleCompromiseCodonTable() {
	sequence := ReadGbk("data/puc19.gbk")
	codonTable := GetCodonTable(11)
	optimizationTable := sequence.GetOptimizationTable(codonTable)

	sequence2 := ReadGbk("data/phix174.gb")
	codonTable2 := GetCodonTable(11)
	optimizationTable2 := sequence2.GetOptimizationTable(codonTable2)

	finalTable, _ := CompromiseCodonTable(optimizationTable, optimizationTable2, 0.1)
	for _, aa := range finalTable.AminoAcids {
		for _, codon := range aa.Codons {
			if codon.Triplet == "TAA" {
				fmt.Println(codon.Weight)
			}
		}
	}
	//output: 2727
}

func TestCompromiseCodonTable(t *testing.T) {
	sequence := ReadGbk("data/puc19.gbk")
	codonTable := GetCodonTable(11)
	optimizationTable := sequence.GetOptimizationTable(codonTable)

	sequence2 := ReadGbk("data/phix174.gb")
	codonTable2 := GetCodonTable(11)
	optimizationTable2 := sequence2.GetOptimizationTable(codonTable2)

	_, err := CompromiseCodonTable(optimizationTable, optimizationTable2, -1.0) // Fails too low
	if err == nil {
		t.Errorf("Compromise table should fail on -1.0")
	}
	_, err = CompromiseCodonTable(optimizationTable, optimizationTable2, 10.0) // Fails too high
	if err == nil {
		t.Errorf("Compromise table should fail on 10.0")
	}
}

func ExampleAddCodonTable() {
	sequence := ReadGbk("data/puc19.gbk")
	codonTable := GetCodonTable(11)
	optimizationTable := sequence.GetOptimizationTable(codonTable)

	sequence2 := ReadGbk("data/phix174.gb")
	codonTable2 := GetCodonTable(11)
	optimizationTable2 := sequence2.GetOptimizationTable(codonTable2)

	finalTable := AddCodonTable(optimizationTable, optimizationTable2)
	for _, aa := range finalTable.AminoAcids {
		for _, codon := range aa.Codons {
			if codon.Triplet == "GGC" {
				fmt.Println(codon.Weight)
			}
		}
	}
	//output: 90
}
