package codon

import (
	"errors"
	"os"
	"strings"
	"testing"

	"github.com/TimothyStiles/poly/io/genbank"
	"github.com/google/go-cmp/cmp"
	"github.com/google/go-cmp/cmp/cmpopts"
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

func TestGenerateCodonTable(t *testing.T) {
	t.Parallel()

	sortStrings := cmpopts.SortSlices((func(a, b string) bool { return a < b }))
	sortAminoAcids := cmpopts.SortSlices(func(a, b AminoAcid) bool { return a.Letter < b.Letter })
	sortCodons := cmpopts.SortSlices(func(a, b Codon) bool { return a.Triplet < b.Triplet })

	tests := []struct {
		name string

		aminoAcids string
		starts     string

		wantCodonTable codonTable
	}{
		{
			name: "ok",

			aminoAcids: "FFLLSSSSYY**CC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG",
			starts:     "---M------**--*----M---------------M----------------------------",

			wantCodonTable: codonTable{
				StartCodons: []string{
					"TTG",
					"CTG",
					"ATG",
				},
				StopCodons: []string{
					"TAA",
					"TAG",
					"TGA",
				},
				AminoAcids: []AminoAcid{
					{
						Letter: "Y",
						Codons: []Codon{{Triplet: "TAT", Weight: 1}, {Triplet: "TAC", Weight: 1}},
					},
					{
						Letter: "*",
						Codons: []Codon{
							{Triplet: "TAA", Weight: 1}, {Triplet: "TAG", Weight: 1},
							{Triplet: "TGA", Weight: 1},
						},
					},
					{
						Letter: "Q",
						Codons: []Codon{{Triplet: "CAA", Weight: 1}, {Triplet: "CAG", Weight: 1}},
					},
					{
						Letter: "I",
						Codons: []Codon{
							{Triplet: "ATT", Weight: 1}, {Triplet: "ATC", Weight: 1},
							{Triplet: "ATA", Weight: 1},
						},
					},
					{Letter: "M", Codons: []Codon{{Triplet: "ATG", Weight: 1}}},
					{
						Letter: "A",
						Codons: []Codon{
							{Triplet: "GCT", Weight: 1}, {Triplet: "GCC", Weight: 1},
							{Triplet: "GCA", Weight: 1}, {Triplet: "GCG", Weight: 1},
						},
					},
					{
						Letter: "D",
						Codons: []Codon{{Triplet: "GAT", Weight: 1}, {Triplet: "GAC", Weight: 1}},
					},
					{
						Letter: "S",
						Codons: []Codon{
							{Triplet: "TCT", Weight: 1}, {Triplet: "TCC", Weight: 1},
							{Triplet: "TCA", Weight: 1}, {Triplet: "TCG", Weight: 1},
							{Triplet: "AGT", Weight: 1}, {Triplet: "AGC", Weight: 1},
						},
					},
					{
						Letter: "C",
						Codons: []Codon{{Triplet: "TGT", Weight: 1}, {Triplet: "TGC", Weight: 1}},
					},
					{
						Letter: "R",
						Codons: []Codon{
							{Triplet: "CGT", Weight: 1}, {Triplet: "CGC", Weight: 1},
							{Triplet: "CGA", Weight: 1}, {Triplet: "CGG", Weight: 1},
							{Triplet: "AGA", Weight: 1}, {Triplet: "AGG", Weight: 1},
						},
					},
					{
						Letter: "N",
						Codons: []Codon{{Triplet: "AAT", Weight: 1}, {Triplet: "AAC", Weight: 1}},
					},
					{
						Letter: "V",
						Codons: []Codon{
							{Triplet: "GTT", Weight: 1}, {Triplet: "GTC", Weight: 1},
							{Triplet: "GTA", Weight: 1}, {Triplet: "GTG", Weight: 1},
						},
					},
					{
						Letter: "F",
						Codons: []Codon{{Triplet: "TTT", Weight: 1}, {Triplet: "TTC", Weight: 1}},
					},
					{
						Letter: "L",
						Codons: []Codon{
							{Triplet: "TTA", Weight: 1}, {Triplet: "TTG", Weight: 1},
							{Triplet: "CTT", Weight: 1}, {Triplet: "CTC", Weight: 1},
							{Triplet: "CTA", Weight: 1}, {Triplet: "CTG", Weight: 1},
						},
					},
					{
						Letter: "P",
						Codons: []Codon{
							{Triplet: "CCT", Weight: 1}, {Triplet: "CCC", Weight: 1},
							{Triplet: "CCA", Weight: 1}, {Triplet: "CCG", Weight: 1},
						},
					},
					{Letter: "W", Codons: []Codon{{Triplet: "TGG", Weight: 1}}},
					{
						Letter: "H",
						Codons: []Codon{{Triplet: "CAT", Weight: 1}, {Triplet: "CAC", Weight: 1}},
					},
					{
						Letter: "T",
						Codons: []Codon{
							{Triplet: "ACT", Weight: 1}, {Triplet: "ACC", Weight: 1},
							{Triplet: "ACA", Weight: 1}, {Triplet: "ACG", Weight: 1},
						},
					},
					{
						Letter: "K",
						Codons: []Codon{{Triplet: "AAA", Weight: 1}, {Triplet: "AAG", Weight: 1}},
					},
					{
						Letter: "E",
						Codons: []Codon{{Triplet: "GAA", Weight: 1}, {Triplet: "GAG", Weight: 1}},
					},
					{
						Letter: "G",
						Codons: []Codon{
							{Triplet: "GGT", Weight: 1}, {Triplet: "GGC", Weight: 1},
							{Triplet: "GGA", Weight: 1}, {Triplet: "GGG", Weight: 1},
						},
					},
				},
				Stats: &Stats{
					StartCodonCount: map[string]int{
						"ATG": 1,
						"CTG": 1,
						"TTG": 1,
					},
				},
			},
		},
	}

	for _, tt := range tests {
		var tt = tt
		t.Run(tt.name, func(t *testing.T) {
			t.Parallel()

			got := generateCodonTable(tt.aminoAcids, tt.starts)

			if !cmp.Equal(got, tt.wantCodonTable, sortCodons, sortAminoAcids, sortStrings) {
				t.Errorf("got and tt.wantCodonTable didn't match %s", cmp.Diff(got, tt.wantCodonTable))
			}
		})
	}
}
