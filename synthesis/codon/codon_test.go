package codon

import (
	"errors"
	"os"
	"strings"
	"testing"

	"github.com/TimothyStiles/poly/io/genbank"
	"github.com/google/go-cmp/cmp"
	weightedRand "github.com/mroth/weightedrand"
	"github.com/stretchr/testify/assert"
)

func TestTranslation(t *testing.T) {
	gfpTranslation := "MASKGEELFTGVVPILVELDGDVNGHKFSVSGEGEGDATYGKLTLKFICTTGKLPVPWPTLVTTFSYGVQCFSRYPDHMKRHDFFKSAMPEGYVQERTISFKDDGNYKTRAEVKFEGDTLVNRIELKGIDFKEDGNILGHKLEYNYNSHNVYITADKQKNGIKANFKIRHNIEDGSVQLADHYQQNTPIGDGPVLLPDNHYLSTQSALSKDPNEKRDHMVLLEFVTAAGITHGMDELYK*"
	gfpDnaSequence := "ATGGCTAGCAAAGGAGAAGAACTTTTCACTGGAGTTGTCCCAATTCTTGTTGAATTAGATGGTGATGTTAATGGGCACAAATTTTCTGTCAGTGGAGAGGGTGAAGGTGATGCTACATACGGAAAGCTTACCCTTAAATTTATTTGCACTACTGGAAAACTACCTGTTCCATGGCCAACACTTGTCACTACTTTCTCTTATGGTGTTCAATGCTTTTCCCGTTATCCGGATCATATGAAACGGCATGACTTTTTCAAGAGTGCCATGCCCGAAGGTTATGTACAGGAACGCACTATATCTTTCAAAGATGACGGGAACTACAAGACGCGTGCTGAAGTCAAGTTTGAAGGTGATACCCTTGTTAATCGTATCGAGTTAAAAGGTATTGATTTTAAAGAAGATGGAAACATTCTCGGACACAAACTCGAGTACAACTATAACTCACACAATGTATACATCACGGCAGACAAACAAAAGAATGGAATCAAAGCTAACTTCAAAATTCGCCACAACATTGAAGATGGATCCGTTCAACTAGCAGACCATTATCAACAAAATACTCCAATTGGCGATGGCCCTGTCCTTTTACCAGACAACCATTACCTGTCGACACAATCTGCCCTTTCGAAAGATCCCAACGAAAAGCGTGACCACATGGTCCTTCTTGAGTTTGTAACTGCTGCTGGGATTACACATGGCATGGATGAGCTCTACAAATAA"

	if got, _ := Translate(gfpDnaSequence, GetCodonTable(11)); got != gfpTranslation {
		t.Errorf("TestTranslation has failed. Translate has returned %q, want %q", got, gfpTranslation)
	}
}

// Non-Met start codons should still map to Met for our standard codon tables.
// See https://github.com/TimothyStiles/poly/issues/305.
func TestTranslationAlwaysMapsStartCodonToMet(t *testing.T) {
	gfpTranslation := "MASKGEELFTGVVPILVELDGDVNGHKFSVSGEGEGDATYGKLTLKFICTTGKLPVPWPTLVTTFSYGVQCFSRYPDHMKRHDFFKSAMPEGYVQERTISFKDDGNYKTRAEVKFEGDTLVNRIELKGIDFKEDGNILGHKLEYNYNSHNVYITADKQKNGIKANFKIRHNIEDGSVQLADHYQQNTPIGDGPVLLPDNHYLSTQSALSKDPNEKRDHMVLLEFVTAAGITHGMDELYK*"
	gfpDnaSequence := "TTGGCTAGCAAAGGAGAAGAACTTTTCACTGGAGTTGTCCCAATTCTTGTTGAATTAGATGGTGATGTTAATGGGCACAAATTTTCTGTCAGTGGAGAGGGTGAAGGTGATGCTACATACGGAAAGCTTACCCTTAAATTTATTTGCACTACTGGAAAACTACCTGTTCCATGGCCAACACTTGTCACTACTTTCTCTTATGGTGTTCAATGCTTTTCCCGTTATCCGGATCATATGAAACGGCATGACTTTTTCAAGAGTGCCATGCCCGAAGGTTATGTACAGGAACGCACTATATCTTTCAAAGATGACGGGAACTACAAGACGCGTGCTGAAGTCAAGTTTGAAGGTGATACCCTTGTTAATCGTATCGAGTTAAAAGGTATTGATTTTAAAGAAGATGGAAACATTCTCGGACACAAACTCGAGTACAACTATAACTCACACAATGTATACATCACGGCAGACAAACAAAAGAATGGAATCAAAGCTAACTTCAAAATTCGCCACAACATTGAAGATGGATCCGTTCAACTAGCAGACCATTATCAACAAAATACTCCAATTGGCGATGGCCCTGTCCTTTTACCAGACAACCATTACCTGTCGACACAATCTGCCCTTTCGAAAGATCCCAACGAAAAGCGTGACCACATGGTCCTTCTTGAGTTTGTAACTGCTGCTGGGATTACACATGGCATGGATGAGCTCTACAAATAA"

	if got, _ := Translate(gfpDnaSequence, GetCodonTable(11)); got != gfpTranslation {
		t.Errorf("TestTranslation has failed. Translate has returned %q, want %q", got, gfpTranslation)
	}
}

func TestTranslationErrorsOnIncorrectStartCodon(t *testing.T) {
	badSequence := "GGG"

	if _, gotErr := Translate(badSequence, GetCodonTable(11)); gotErr == nil {
		t.Errorf("Translation should return an error if given an incorrect start codon")
	}
}

func TestTranslationErrorsOnEmptyCodonTable(t *testing.T) {
	emtpyCodonTable := codonTable{}
	_, err := Translate("A", emtpyCodonTable)

	if err != errEmptyCodonTable {
		t.Error("Translation should return an error if given an empty codon table")
	}
}

func TestTranslationErrorsOnEmptyAminoAcidString(t *testing.T) {
	nonEmptyCodonTable := GetCodonTable(1)
	_, err := Translate("", nonEmptyCodonTable)

	if err != errEmptySequenceString {
		t.Error("Translation should return an error if given an empty sequence string")
	}
}

func TestTranslationMixedCase(t *testing.T) {
	gfpTranslation := "MASKGEELFTGVVPILVELDGDVNGHKFSVSGEGEGDATYGKLTLKFICTTGKLPVPWPTLVTTFSYGVQCFSRYPDHMKRHDFFKSAMPEGYVQERTISFKDDGNYKTRAEVKFEGDTLVNRIELKGIDFKEDGNILGHKLEYNYNSHNVYITADKQKNGIKANFKIRHNIEDGSVQLADHYQQNTPIGDGPVLLPDNHYLSTQSALSKDPNEKRDHMVLLEFVTAAGITHGMDELYK*"
	gfpDnaSequence := "atggctagcaaaggagaagaacttttcactggagttgtcccaaTTCTTGTTGAATTAGATGGTGATGTTAATGGGCACAAATTTTCTGTCAGTGGAGAGGGTGAAGGTGATGCTACATACGGAAAGCTTACCCTTAAATTTATTTGCACTACTGGAAAACTACCTGTTCCATGGCCAACACTTGTCACTACTTTCTCTTATGGTGTTCAATGCTTTTCCCGTTATCCGGATCATATGAAACGGCATGACTTTTTCAAGAGTGCCATGCCCGAAGGTTATGTACAGGAACGCACTATATCTTTCAAAGATGACGGGAACTACAAGACGCGTGCTGAAGTCAAGTTTGAAGGTGATACCCTTGTTAATCGTATCGAGTTAAAAGGTATTGATTTTAAAGAAGATGGAAACATTCTCGGACACAAACTCGAGTACAACTATAACTCACACAATGTATACATCACGGCAGACAAACAAAAGAATGGAATCAAAGCTAACTTCAAAATTCGCCACAACATTGAAGATGGATCCGTTCAACTAGCAGACCATTATCAACAAAATACTCCAATTGGCGATGGCCCTGTCCTTTTACCAGACAACCATTACCTGTCGACACAATCTGCCCTTTCGAAAGATCCCAACGAAAAGCGTGACCACATGGTCCTTCTTGAGTTTGTAACTGCTGCTGGGATTACACATGGCATGGATGAGCTCTACAAATAA"
	if got, _ := Translate(gfpDnaSequence, GetCodonTable(11)); got != gfpTranslation {
		t.Errorf("TestTranslationMixedCase has failed. Translate has returned %q, want %q", got, gfpTranslation)
	}
}

func TestTranslationLowerCase(t *testing.T) {
	gfpTranslation := "MASKGEELFTGVVPILVELDGDVNGHKFSVSGEGEGDATYGKLTLKFICTTGKLPVPWPTLVTTFSYGVQCFSRYPDHMKRHDFFKSAMPEGYVQERTISFKDDGNYKTRAEVKFEGDTLVNRIELKGIDFKEDGNILGHKLEYNYNSHNVYITADKQKNGIKANFKIRHNIEDGSVQLADHYQQNTPIGDGPVLLPDNHYLSTQSALSKDPNEKRDHMVLLEFVTAAGITHGMDELYK*"
	gfpDnaSequence := "atggctagcaaaggagaagaacttttcactggagttgtcccaattcttgttgaattagatggtgatgttaatgggcacaaattttctgtcagtggagagggtgaaggtgatgctacatacggaaagcttacccttaaatttatttgcactactggaaaactacctgttccatggccaacacttgtcactactttctcttatggtgttcaatgcttttcccgttatccggatcatatgaaacggcatgactttttcaagagtgccatgcccgaaggttatgtacaggaacgcactatatctttcaaagatgacgggaactacaagacgcgtgctgaagtcaagtttgaaggtgatacccttgttaatcgtatcgagttaaaaggtattgattttaaagaagatggaaacattctcggacacaaactcgagtacaactataactcacacaatgtatacatcacggcagacaaacaaaagaatggaatcaaagctaacttcaaaattcgccacaacattgaagatggatccgttcaactagcagaccattatcaacaaaatactccaattggcgatggccctgtccttttaccagacaaccattacctgtcgacacaatctgccctttcgaaagatcccaacgaaaagcgtgaccacatggtccttcttgagtttgtaactgctgctgggattacacatggcatggatgagctctacaaataa"
	if got, _ := Translate(gfpDnaSequence, GetCodonTable(11)); got != gfpTranslation {
		t.Errorf("TestTranslationLowerCase has failed. Translate has returned %q, want %q", got, gfpTranslation)
	}
}

func TestOptimize(t *testing.T) {
	gfpTranslation := "MASKGEELFTGVVPILVELDGDVNGHKFSVSGEGEGDATYGKLTLKFICTTGKLPVPWPTLVTTFSYGVQCFSRYPDHMKRHDFFKSAMPEGYVQERTISFKDDGNYKTRAEVKFEGDTLVNRIELKGIDFKEDGNILGHKLEYNYNSHNVYITADKQKNGIKANFKIRHNIEDGSVQLADHYQQNTPIGDGPVLLPDNHYLSTQSALSKDPNEKRDHMVLLEFVTAAGITHGMDELYK*"

	sequence, _ := genbank.Read("../../data/puc19.gbk")
	codonTable := GetCodonTable(11)

	// a string builder to build a single concatenated string of all coding regions
	var codingRegionsBuilder strings.Builder

	// iterate through the features of the genbank file and if the feature is a coding region, append the sequence to the string builder
	for _, feature := range sequence.Features {
		if feature.Type == "CDS" {
			sequence, _ := feature.GetSequence()
			codingRegionsBuilder.WriteString(sequence)
		}
	}

	// get the concatenated sequence string of the coding regions
	codingRegions := codingRegionsBuilder.String()

	// weight our codon optimization table using the regions we collected from the genbank file above
	optimizationTable := codonTable.OptimizeTable(codingRegions)

	optimizedSequence, _ := Optimize(gfpTranslation, optimizationTable)
	optimizedSequenceTranslation, _ := Translate(optimizedSequence, optimizationTable)

	if optimizedSequenceTranslation != gfpTranslation {
		t.Errorf("TestOptimize has failed. Translate has returned %q, want %q", optimizedSequenceTranslation, gfpTranslation)
	}
}

func TestOptimizeSameSeed(t *testing.T) {
	var gfpTranslation = "MASKGEELFTGVVPILVELDGDVNGHKFSVSGEGEGDATYGKLTLKFICTTGKLPVPWPTLVTTFSYGVQCFSRYPDHMKRHDFFKSAMPEGYVQERTISFKDDGNYKTRAEVKFEGDTLVNRIELKGIDFKEDGNILGHKLEYNYNSHNVYITADKQKNGIKANFKIRHNIEDGSVQLADHYQQNTPIGDGPVLLPDNHYLSTQSALSKDPNEKRDHMVLLEFVTAAGITHGMDELYK*"
	var sequence, _ = genbank.Read("../../data/puc19.gbk")
	var codonTable = GetCodonTable(11)

	// a string builder to build a single concatenated string of all coding regions
	var codingRegionsBuilder strings.Builder

	// iterate through the features of the genbank file and if the feature is a coding region, append the sequence to the string builder
	for _, feature := range sequence.Features {
		if feature.Type == "CDS" {
			sequence, _ := feature.GetSequence()
			codingRegionsBuilder.WriteString(sequence)
		}
	}

	// get the concatenated sequence string of the coding regions
	codingRegions := codingRegionsBuilder.String()

	var optimizationTable = codonTable.OptimizeTable(codingRegions)
	randomSeed := 10

	optimizedSequence, _ := Optimize(gfpTranslation, optimizationTable, randomSeed)
	otherOptimizedSequence, _ := Optimize(gfpTranslation, optimizationTable, randomSeed)

	if optimizedSequence != otherOptimizedSequence {
		t.Error("Optimized sequence with the same random seed are not the same")
	}
}

func TestOptimizeDifferentSeed(t *testing.T) {
	var gfpTranslation = "MASKGEELFTGVVPILVELDGDVNGHKFSVSGEGEGDATYGKLTLKFICTTGKLPVPWPTLVTTFSYGVQCFSRYPDHMKRHDFFKSAMPEGYVQERTISFKDDGNYKTRAEVKFEGDTLVNRIELKGIDFKEDGNILGHKLEYNYNSHNVYITADKQKNGIKANFKIRHNIEDGSVQLADHYQQNTPIGDGPVLLPDNHYLSTQSALSKDPNEKRDHMVLLEFVTAAGITHGMDELYK*"
	var sequence, _ = genbank.Read("../../data/puc19.gbk")
	var codonTable = GetCodonTable(11)

	// a string builder to build a single concatenated string of all coding regions
	var codingRegionsBuilder strings.Builder

	// iterate through the features of the genbank file and if the feature is a coding region, append the sequence to the string builder
	for _, feature := range sequence.Features {
		if feature.Type == "CDS" {
			sequence, _ := feature.GetSequence()
			codingRegionsBuilder.WriteString(sequence)
		}
	}

	// get the concatenated sequence string of the coding regions
	codingRegions := codingRegionsBuilder.String()

	var optimizationTable = codonTable.OptimizeTable(codingRegions)

	optimizedSequence, _ := Optimize(gfpTranslation, optimizationTable)
	otherOptimizedSequence, _ := Optimize(gfpTranslation, optimizationTable)

	if optimizedSequence == otherOptimizedSequence {
		t.Error("Optimized sequence with different random seed have the same result")
	}
}

func TestOptimizeErrorsOnEmptyCodonTable(t *testing.T) {
	emtpyCodonTable := codonTable{}
	_, err := Optimize("A", emtpyCodonTable)

	if err != errEmptyCodonTable {
		t.Error("Optimize should return an error if given an empty codon table")
	}
}

func TestOptimizeErrorsOnEmptyAminoAcidString(t *testing.T) {
	nonEmptyCodonTable := GetCodonTable(1)
	_, err := Optimize("", nonEmptyCodonTable)

	if err != errEmptyAminoAcidString {
		t.Error("Optimize should return an error if given an empty amino acid string")
	}
}
func TestOptimizeErrorsOnInvalidAminoAcid(t *testing.T) {
	aminoAcids := "TOP"
	table := GetCodonTable(1) // does not contain 'O'

	_, optimizeErr := Optimize(aminoAcids, table)
	assert.EqualError(t, optimizeErr, invalidAminoAcidError{'O'}.Error())
}

func TestOptimizeErrorsOnBrokenChooser(t *testing.T) {
	gfpTranslation := "MASKGEELFTGVVPILVELDGDVNGHKFSVSGEGEGDATYGKLTLKFICTTGKLPVPWPTLVTTFSYGVQCFSRYPDHMKRHDFFKSAMPEGYVQERTISFKDDGNYKTRAEVKFEGDTLVNRIELKGIDFKEDGNILGHKLEYNYNSHNVYITADKQKNGIKANFKIRHNIEDGSVQLADHYQQNTPIGDGPVLLPDNHYLSTQSALSKDPNEKRDHMVLLEFVTAAGITHGMDELYK*"

	chooserErr := errors.New("chooser rigged to fail")

	codonTable := &mockTable{
		ChooserFn: func() (map[string]weightedRand.Chooser, error) {
			return nil, chooserErr
		},
		IsEmptyFn: func() bool {
			return false
		},
	}

	_, err := Optimize(gfpTranslation, codonTable)
	assert.EqualError(t, err, chooserErr.Error())
}

func TestGetCodonFrequency(t *testing.T) {
	translationTable := GetCodonTable(11).GenerateTranslationTable()

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

func TestChooserError(t *testing.T) {
	codonTable := GetCodonTable(11)

	oldChooserFn := newChooserFn
	newChooserFn = func(choices ...weightedRand.Choice) (*weightedRand.Chooser, error) {
		return nil, errors.New("new chooser rigged to fail")
	}
	defer func() {
		newChooserFn = oldChooserFn
	}()

	_, err := codonTable.Chooser()
	assert.EqualError(t, err, "weightedRand.NewChooser() error: new chooser rigged to fail")
}

/******************************************************************************

JSON related tests begin here.

******************************************************************************/

func TestWriteCodonJSON(t *testing.T) {
	testCodonTable := ReadCodonJSON("../../data/bsub_codon_test.json")
	WriteCodonJSON(testCodonTable, "../../data/codon_test1.json")
	readTestCodonTable := ReadCodonJSON("../../data/codon_test1.json")

	// cleaning up test data
	os.Remove("../../data/codon_test1.json")

	if diff := cmp.Diff(testCodonTable, readTestCodonTable); diff != "" {
		t.Errorf(" mismatch (-want +got):\n%s", diff)
	}
}

/*
*****************************************************************************

Codon Compromise + Add related tests begin here.

*****************************************************************************
*/
func TestCompromiseCodonTable(t *testing.T) {
	sequence, _ := genbank.Read("../../data/puc19.gbk")
	codonTable := GetCodonTable(11)

	// a string builder to build a single concatenated string of all coding regions
	var codingRegionsBuilder strings.Builder

	// iterate through the features of the genbank file and if the feature is a coding region, append the sequence to the string builder
	for _, feature := range sequence.Features {
		if feature.Type == "CDS" {
			sequence, _ := feature.GetSequence()
			codingRegionsBuilder.WriteString(sequence)
		}
	}

	// get the concatenated sequence string of the coding regions
	codingRegions := codingRegionsBuilder.String()

	// weight our codon optimization table using the regions we collected from the genbank file above
	optimizationTable := codonTable.OptimizeTable(codingRegions)

	sequence2, _ := genbank.Read("../../data/phix174.gb")
	codonTable2 := GetCodonTable(11)

	// a string builder to build a single concatenated string of all coding regions
	var codingRegionsBuilder2 strings.Builder

	// iterate through the features of the genbank file and if the feature is a coding region, append the sequence to the string builder
	for _, feature := range sequence2.Features {
		if feature.Type == "CDS" {
			sequence, _ := feature.GetSequence()
			codingRegionsBuilder2.WriteString(sequence)
		}
	}

	// get the concatenated sequence string of the coding regions
	codingRegions2 := codingRegionsBuilder2.String()

	// weight our codon optimization table using the regions we collected from the genbank file above
	optimizationTable2 := codonTable2.OptimizeTable(codingRegions2)

	_, err := CompromiseCodonTable(optimizationTable, optimizationTable2, -1.0) // Fails too low
	if err == nil {
		t.Errorf("Compromise table should fail on -1.0")
	}
	_, err = CompromiseCodonTable(optimizationTable, optimizationTable2, 10.0) // Fails too high
	if err == nil {
		t.Errorf("Compromise table should fail on 10.0")
	}
}

type mockTable struct {
	codonTable
	ChooserFn func() (map[string]weightedRand.Chooser, error)
	IsEmptyFn func() bool
}

func (t *mockTable) Chooser() (map[string]weightedRand.Chooser, error) {
	return t.ChooserFn()
}

func (t *mockTable) IsEmpty() bool {
	return t.IsEmptyFn()
}

func TestCapitalizationRegression(t *testing.T) {
	// Tests to make sure that amino acids are capitalized
	gfpTranslation := "MaSKGEELFTGVVPILVELDGDVNGHKFSVSGEGEGDATYGKLTLKFICTTGKLPVPWPTLVTTFSYGVQCFSRYPDHMKRHDFFKSAMPEGYVQERTISFKDDGNYKTRAEVKFEGDTLVNRIELKGIDFKEDGNILGHKLEYNYNSHNVYITADKQKNGIKANFKIRHNIEDGSVQLADHYQQNTPIGDGPVLLPDNHYLSTQSALSKDPNEKRDHMVLLEFVTAAGITHGMDELYK*"

	sequence, _ := genbank.Read("../../data/puc19.gbk")
	codonTable := GetCodonTable(11)

	// a string builder to build a single concatenated string of all coding regions
	var codingRegionsBuilder strings.Builder

	// iterate through the features of the genbank file and if the feature is a coding region, append the sequence to the string builder
	for _, feature := range sequence.Features {
		if feature.Type == "CDS" {
			sequence, _ := feature.GetSequence()
			codingRegionsBuilder.WriteString(sequence)
		}
	}

	// get the concatenated sequence string of the coding regions
	codingRegions := codingRegionsBuilder.String()

	// weight our codon optimization table using the regions we collected from the genbank file above
	optimizationTable := codonTable.OptimizeTable(codingRegions)

	optimizedSequence, _ := Optimize(gfpTranslation, optimizationTable)
	optimizedSequenceTranslation, _ := Translate(optimizedSequence, optimizationTable)

	if optimizedSequenceTranslation != strings.ToUpper(gfpTranslation) {
		t.Errorf("TestOptimize has failed. Translate has returned %q, want %q", optimizedSequenceTranslation, gfpTranslation)
	}
}
