package poly

import (
	"strings"
	"testing"
)

func TestTranslation(t *testing.T) {
	gfpTranslation := "MASKGEELFTGVVPILVELDGDVNGHKFSVSGEGEGDATYGKLTLKFICTTGKLPVPWPTLVTTFSYGVQCFSRYPDHMKRHDFFKSAMPEGYVQERTISFKDDGNYKTRAEVKFEGDTLVNRIELKGIDFKEDGNILGHKLEYNYNSHNVYITADKQKNGIKANFKIRHNIEDGSVQLADHYQQNTPIGDGPVLLPDNHYLSTQSALSKDPNEKRDHMVLLEFVTAAGITHGMDELYK*"
	gfpDnaSequence := "ATGGCTAGCAAAGGAGAAGAACTTTTCACTGGAGTTGTCCCAATTCTTGTTGAATTAGATGGTGATGTTAATGGGCACAAATTTTCTGTCAGTGGAGAGGGTGAAGGTGATGCTACATACGGAAAGCTTACCCTTAAATTTATTTGCACTACTGGAAAACTACCTGTTCCATGGCCAACACTTGTCACTACTTTCTCTTATGGTGTTCAATGCTTTTCCCGTTATCCGGATCATATGAAACGGCATGACTTTTTCAAGAGTGCCATGCCCGAAGGTTATGTACAGGAACGCACTATATCTTTCAAAGATGACGGGAACTACAAGACGCGTGCTGAAGTCAAGTTTGAAGGTGATACCCTTGTTAATCGTATCGAGTTAAAAGGTATTGATTTTAAAGAAGATGGAAACATTCTCGGACACAAACTCGAGTACAACTATAACTCACACAATGTATACATCACGGCAGACAAACAAAAGAATGGAATCAAAGCTAACTTCAAAATTCGCCACAACATTGAAGATGGATCCGTTCAACTAGCAGACCATTATCAACAAAATACTCCAATTGGCGATGGCCCTGTCCTTTTACCAGACAACCATTACCTGTCGACACAATCTGCCCTTTCGAAAGATCCCAACGAAAAGCGTGACCACATGGTCCTTCTTGAGTTTGTAACTGCTGCTGGGATTACACATGGCATGGATGAGCTCTACAAATAA"
	if got := Translate(gfpDnaSequence, defaultCodonTablesByNumber[11]); got != gfpTranslation {
		t.Errorf("TestTranslation has failed. Translate has returned %q, want %q", got, gfpTranslation)
	}

}

func TestOptimize(t *testing.T) {
	gfpTranslation := "MASKGEELFTGVVPILVELDGDVNGHKFSVSGEGEGDATYGKLTLKFICTTGKLPVPWPTLVTTFSYGVQCFSRYPDHMKRHDFFKSAMPEGYVQERTISFKDDGNYKTRAEVKFEGDTLVNRIELKGIDFKEDGNILGHKLEYNYNSHNVYITADKQKNGIKANFKIRHNIEDGSVQLADHYQQNTPIGDGPVLLPDNHYLSTQSALSKDPNEKRDHMVLLEFVTAAGITHGMDELYK*"
	codonTable := defaultCodonTablesByNumber[11]

	sequence := ReadGbk("data/bsub.gbk")
	rawSequence := sequence.Sequence
	codonTable.CreateWeights(rawSequence)

	optimizedSequence := Optimize(gfpTranslation, codonTable)
	optimizedSequenceTranslation := Translate(optimizedSequence, codonTable)

	if optimizedSequenceTranslation != gfpTranslation {
		t.Errorf("TestOptimize has failed. Translate has returned %q, want %q", optimizedSequenceTranslation, gfpTranslation)
	}

}

func TestGetCodonFrequency(t *testing.T) {

	translationTable := defaultCodonTablesByNumber[11].generateTranslationTable()

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
