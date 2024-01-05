package codon

import (
	"errors"
	"os"
	"strings"
	"testing"

	"github.com/bebop/poly/io/genbank"
	"github.com/google/go-cmp/cmp"
	"github.com/google/go-cmp/cmp/cmpopts"
	weightedRand "github.com/mroth/weightedrand"
	"github.com/stretchr/testify/assert"
)

func TestTranslation(t *testing.T) {
	gfpTranslation := "MASKGEELFTGVVPILVELDGDVNGHKFSVSGEGEGDATYGKLTLKFICTTGKLPVPWPTLVTTFSYGVQCFSRYPDHMKRHDFFKSAMPEGYVQERTISFKDDGNYKTRAEVKFEGDTLVNRIELKGIDFKEDGNILGHKLEYNYNSHNVYITADKQKNGIKANFKIRHNIEDGSVQLADHYQQNTPIGDGPVLLPDNHYLSTQSALSKDPNEKRDHMVLLEFVTAAGITHGMDELYK*"
	gfpDnaSequence := "ATGGCTAGCAAAGGAGAAGAACTTTTCACTGGAGTTGTCCCAATTCTTGTTGAATTAGATGGTGATGTTAATGGGCACAAATTTTCTGTCAGTGGAGAGGGTGAAGGTGATGCTACATACGGAAAGCTTACCCTTAAATTTATTTGCACTACTGGAAAACTACCTGTTCCATGGCCAACACTTGTCACTACTTTCTCTTATGGTGTTCAATGCTTTTCCCGTTATCCGGATCATATGAAACGGCATGACTTTTTCAAGAGTGCCATGCCCGAAGGTTATGTACAGGAACGCACTATATCTTTCAAAGATGACGGGAACTACAAGACGCGTGCTGAAGTCAAGTTTGAAGGTGATACCCTTGTTAATCGTATCGAGTTAAAAGGTATTGATTTTAAAGAAGATGGAAACATTCTCGGACACAAACTCGAGTACAACTATAACTCACACAATGTATACATCACGGCAGACAAACAAAAGAATGGAATCAAAGCTAACTTCAAAATTCGCCACAACATTGAAGATGGATCCGTTCAACTAGCAGACCATTATCAACAAAATACTCCAATTGGCGATGGCCCTGTCCTTTTACCAGACAACCATTACCTGTCGACACAATCTGCCCTTTCGAAAGATCCCAACGAAAAGCGTGACCACATGGTCCTTCTTGAGTTTGTAACTGCTGCTGGGATTACACATGGCATGGATGAGCTCTACAAATAA"

	table, err := NewTranslationTable(11)
	if err != nil {
		t.Fatalf("failed to initialise codon table: %s", err)
	}

	if got, _ := table.Translate(gfpDnaSequence); got != gfpTranslation {
		t.Errorf("TestTranslation has failed. Translate has returned %q, want %q", got, gfpTranslation)
	}
}

func TestTranslationErrorsOnEmptyAminoAcidString(t *testing.T) {
	nonEmptyCodonTable, err := NewTranslationTable(1)
	if err != nil {
		t.Fatalf("failed to initialise codon table: %s", err)
	}

	_, err = nonEmptyCodonTable.Translate("")

	if err != errEmptySequenceString {
		t.Error("Translation should return an error if given an empty sequence string")
	}
}

func TestTranslationMixedCase(t *testing.T) {
	gfpTranslation := "MASKGEELFTGVVPILVELDGDVNGHKFSVSGEGEGDATYGKLTLKFICTTGKLPVPWPTLVTTFSYGVQCFSRYPDHMKRHDFFKSAMPEGYVQERTISFKDDGNYKTRAEVKFEGDTLVNRIELKGIDFKEDGNILGHKLEYNYNSHNVYITADKQKNGIKANFKIRHNIEDGSVQLADHYQQNTPIGDGPVLLPDNHYLSTQSALSKDPNEKRDHMVLLEFVTAAGITHGMDELYK*"
	gfpDnaSequence := "atggctagcaaaggagaagaacttttcactggagttgtcccaaTTCTTGTTGAATTAGATGGTGATGTTAATGGGCACAAATTTTCTGTCAGTGGAGAGGGTGAAGGTGATGCTACATACGGAAAGCTTACCCTTAAATTTATTTGCACTACTGGAAAACTACCTGTTCCATGGCCAACACTTGTCACTACTTTCTCTTATGGTGTTCAATGCTTTTCCCGTTATCCGGATCATATGAAACGGCATGACTTTTTCAAGAGTGCCATGCCCGAAGGTTATGTACAGGAACGCACTATATCTTTCAAAGATGACGGGAACTACAAGACGCGTGCTGAAGTCAAGTTTGAAGGTGATACCCTTGTTAATCGTATCGAGTTAAAAGGTATTGATTTTAAAGAAGATGGAAACATTCTCGGACACAAACTCGAGTACAACTATAACTCACACAATGTATACATCACGGCAGACAAACAAAAGAATGGAATCAAAGCTAACTTCAAAATTCGCCACAACATTGAAGATGGATCCGTTCAACTAGCAGACCATTATCAACAAAATACTCCAATTGGCGATGGCCCTGTCCTTTTACCAGACAACCATTACCTGTCGACACAATCTGCCCTTTCGAAAGATCCCAACGAAAAGCGTGACCACATGGTCCTTCTTGAGTTTGTAACTGCTGCTGGGATTACACATGGCATGGATGAGCTCTACAAATAA"
	table, err := NewTranslationTable(11)
	if err != nil {
		t.Fatalf("failed to initialise codon table: %s", err)
	}

	if got, _ := table.Translate(gfpDnaSequence); got != gfpTranslation {
		t.Errorf("TestTranslationMixedCase has failed. Translate has returned %q, want %q", got, gfpTranslation)
	}
}

func TestTranslationLowerCase(t *testing.T) {
	gfpTranslation := "MASKGEELFTGVVPILVELDGDVNGHKFSVSGEGEGDATYGKLTLKFICTTGKLPVPWPTLVTTFSYGVQCFSRYPDHMKRHDFFKSAMPEGYVQERTISFKDDGNYKTRAEVKFEGDTLVNRIELKGIDFKEDGNILGHKLEYNYNSHNVYITADKQKNGIKANFKIRHNIEDGSVQLADHYQQNTPIGDGPVLLPDNHYLSTQSALSKDPNEKRDHMVLLEFVTAAGITHGMDELYK*"
	gfpDnaSequence := "atggctagcaaaggagaagaacttttcactggagttgtcccaattcttgttgaattagatggtgatgttaatgggcacaaattttctgtcagtggagagggtgaaggtgatgctacatacggaaagcttacccttaaatttatttgcactactggaaaactacctgttccatggccaacacttgtcactactttctcttatggtgttcaatgcttttcccgttatccggatcatatgaaacggcatgactttttcaagagtgccatgcccgaaggttatgtacaggaacgcactatatctttcaaagatgacgggaactacaagacgcgtgctgaagtcaagtttgaaggtgatacccttgttaatcgtatcgagttaaaaggtattgattttaaagaagatggaaacattctcggacacaaactcgagtacaactataactcacacaatgtatacatcacggcagacaaacaaaagaatggaatcaaagctaacttcaaaattcgccacaacattgaagatggatccgttcaactagcagaccattatcaacaaaatactccaattggcgatggccctgtccttttaccagacaaccattacctgtcgacacaatctgccctttcgaaagatcccaacgaaaagcgtgaccacatggtccttcttgagtttgtaactgctgctgggattacacatggcatggatgagctctacaaataa"

	table, err := NewTranslationTable(11)
	if err != nil {
		t.Fatalf("failed to initialise codon table: %s", err)
	}

	if got, _ := table.Translate(gfpDnaSequence); got != gfpTranslation {
		t.Errorf("TestTranslationLowerCase has failed. Translate has returned %q, want %q", got, gfpTranslation)
	}
}

func TestOptimize(t *testing.T) {
	gfpTranslation := "MASKGEELFTGVVPILVELDGDVNGHKFSVSGEGEGDATYGKLTLKFICTTGKLPVPWPTLVTTFSYGVQCFSRYPDHMKRHDFFKSAMPEGYVQERTISFKDDGNYKTRAEVKFEGDTLVNRIELKGIDFKEDGNILGHKLEYNYNSHNVYITADKQKNGIKANFKIRHNIEDGSVQLADHYQQNTPIGDGPVLLPDNHYLSTQSALSKDPNEKRDHMVLLEFVTAAGITHGMDELYK*"

	sequence, _ := genbank.Read("../../data/puc19.gbk")

	table, err := NewTranslationTable(11)
	if err != nil {
		t.Fatalf("failed to initialise codon table: %s", err)
	}

	err = table.UpdateWeightsWithSequence(sequence)
	if err != nil {
		t.Error(err)
	}

	codonTable, err := NewTranslationTable(11)
	if err != nil {
		t.Fatalf("failed to initialise codon table: %s", err)
	}

	optimizedSequence, _ := table.Optimize(gfpTranslation)
	optimizedSequenceTranslation, _ := codonTable.Translate(optimizedSequence)

	if optimizedSequenceTranslation != gfpTranslation {
		t.Errorf("TestOptimize has failed. Translate has returned %q, want %q", optimizedSequenceTranslation, gfpTranslation)
	}
}

func TestOptimizeSameSeed(t *testing.T) {
	var gfpTranslation = "MASKGEELFTGVVPILVELDGDVNGHKFSVSGEGEGDATYGKLTLKFICTTGKLPVPWPTLVTTFSYGVQCFSRYPDHMKRHDFFKSAMPEGYVQERTISFKDDGNYKTRAEVKFEGDTLVNRIELKGIDFKEDGNILGHKLEYNYNSHNVYITADKQKNGIKANFKIRHNIEDGSVQLADHYQQNTPIGDGPVLLPDNHYLSTQSALSKDPNEKRDHMVLLEFVTAAGITHGMDELYK*"
	var sequence, _ = genbank.Read("../../data/puc19.gbk")
	optimizationTable, err := NewTranslationTable(11)
	if err != nil {
		t.Fatalf("failed to initialise codon table: %s", err)
	}

	err = optimizationTable.UpdateWeightsWithSequence(sequence)
	if err != nil {
		t.Error(err)
	}
	if err != nil {
		t.Error(err)
	}

	randomSeed := 10

	optimizedSequence, _ := optimizationTable.Optimize(gfpTranslation, randomSeed)
	otherOptimizedSequence, _ := optimizationTable.Optimize(gfpTranslation, randomSeed)

	if optimizedSequence != otherOptimizedSequence {
		t.Error("Optimized sequence with the same random seed are not the same")
	}
}

func TestOptimizeDifferentSeed(t *testing.T) {
	var gfpTranslation = "MASKGEELFTGVVPILVELDGDVNGHKFSVSGEGEGDATYGKLTLKFICTTGKLPVPWPTLVTTFSYGVQCFSRYPDHMKRHDFFKSAMPEGYVQERTISFKDDGNYKTRAEVKFEGDTLVNRIELKGIDFKEDGNILGHKLEYNYNSHNVYITADKQKNGIKANFKIRHNIEDGSVQLADHYQQNTPIGDGPVLLPDNHYLSTQSALSKDPNEKRDHMVLLEFVTAAGITHGMDELYK*"
	var sequence, _ = genbank.Read("../../data/puc19.gbk")
	optimizationTable, err := NewTranslationTable(11)
	if err != nil {
		t.Fatalf("failed to initialise codon table: %s", err)
	}

	err = optimizationTable.UpdateWeightsWithSequence(sequence)
	if err != nil {
		t.Error(err)
	}

	optimizedSequence, _ := optimizationTable.Optimize(gfpTranslation)
	otherOptimizedSequence, _ := optimizationTable.Optimize(gfpTranslation)

	if optimizedSequence == otherOptimizedSequence {
		t.Error("Optimized sequence with different random seed have the same result")
	}
}

func TestOptimizeErrorsOnEmptyAminoAcidString(t *testing.T) {
	nonEmptyCodonTable, err := NewTranslationTable(1)
	if err != nil {
		t.Fatalf("failed to initialise codon table: %s", err)
	}

	_, err = nonEmptyCodonTable.Optimize("")

	if err != errEmptyAminoAcidString {
		t.Error("Optimize should return an error if given an empty amino acid string")
	}
}
func TestOptimizeErrorsOnInvalidAminoAcid(t *testing.T) {
	aminoAcids := "TOP"
	table, err := NewTranslationTable(1)
	if err != nil {
		t.Fatalf("failed to initialise codon table: %s", err)
	}
	// does not contain 'O'

	_, optimizeErr := table.Optimize(aminoAcids)
	assert.EqualError(t, optimizeErr, invalidAminoAcidError{'O'}.Error())
}

func TestGetCodonFrequency(t *testing.T) {
	table, err := NewTranslationTable(11)
	if err != nil {
		t.Fatalf("failed to initialise codon table: %s", err)
	}
	translationTable := table.TranslationMap

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

	// weight our codon optimization table using the regions we collected from the genbank file above

	optimizationTable, err := NewTranslationTable(11)
	if err != nil {
		t.Fatalf("failed to initialise codon table: %s", err)
	}

	err = optimizationTable.UpdateWeightsWithSequence(sequence)
	if err != nil {
		t.Error(err)
	}

	sequence2, _ := genbank.Read("../../data/phix174.gb")
	optimizationTable2, err := NewTranslationTable(11)
	if err != nil {
		t.Fatalf("failed to initialise codon table: %s", err)
	}

	err = optimizationTable2.UpdateWeightsWithSequence(sequence2)
	if err != nil {
		t.Error(err)
	}

	_, err = CompromiseCodonTable(optimizationTable, optimizationTable2, -1.0) // Fails too low
	if err == nil {
		t.Errorf("Compromise table should fail on -1.0")
	}
	_, err = CompromiseCodonTable(optimizationTable, optimizationTable2, 10.0) // Fails too high
	if err == nil {
		t.Errorf("Compromise table should fail on 10.0")
	}

	// replace chooser fn with test one
	newChooserFn = func(choices ...weightedRand.Choice) (*weightedRand.Chooser, error) {
		return nil, errors.New("new chooser rigged to fail")
	}

	defer func() {
		newChooserFn = weightedRand.NewChooser
	}()

	_, err = CompromiseCodonTable(optimizationTable, optimizationTable2, 0.1)
	if err == nil {
		t.Errorf("Compromise table should fail when new chooser func rigged")
	}
}

func TestAddCodonTable(t *testing.T) {
	sequence, _ := genbank.Read("../../data/puc19.gbk")

	// weight our codon optimization table using the regions we collected from the genbank file above

	optimizationTable, err := NewTranslationTable(11)
	if err != nil {
		t.Fatalf("failed to initialise codon table: %s", err)
	}

	err = optimizationTable.UpdateWeightsWithSequence(sequence)
	if err != nil {
		t.Error(err)
	}

	sequence2, _ := genbank.Read("../../data/phix174.gb")
	optimizationTable2, err := NewTranslationTable(11)
	if err != nil {
		t.Fatalf("failed to initialise codon table: %s", err)
	}

	err = optimizationTable2.UpdateWeightsWithSequence(sequence2)
	if err != nil {
		t.Error(err)
	}

	// replace chooser fn with test one
	newChooserFn = func(choices ...weightedRand.Choice) (*weightedRand.Chooser, error) {
		return nil, errors.New("new chooser rigged to fail")
	}

	defer func() {
		newChooserFn = weightedRand.NewChooser
	}()

	_, err = AddCodonTable(optimizationTable, optimizationTable2)
	if err == nil {
		t.Errorf("Compromise table should fail when new chooser func rigged")
	}
}

func TestCapitalizationRegression(t *testing.T) {
	// Tests to make sure that amino acids are capitalized
	gfpTranslation := "MaSKGEELFTGVVPILVELDGDVNGHKFSVSGEGEGDATYGKLTLKFICTTGKLPVPWPTLVTTFSYGVQCFSRYPDHMKRHDFFKSAMPEGYVQERTISFKDDGNYKTRAEVKFEGDTLVNRIELKGIDFKEDGNILGHKLEYNYNSHNVYITADKQKNGIKANFKIRHNIEDGSVQLADHYQQNTPIGDGPVLLPDNHYLSTQSALSKDPNEKRDHMVLLEFVTAAGITHGMDELYK*"

	sequence, _ := genbank.Read("../../data/puc19.gbk")

	optimizationTable, err := NewTranslationTable(11)
	if err != nil {
		t.Fatalf("failed to initialise codon table: %s", err)
	}

	err = optimizationTable.UpdateWeightsWithSequence(sequence)
	if err != nil {
		t.Error(err)
	}

	optimizedSequence, _ := optimizationTable.Optimize(gfpTranslation, 1)
	optimizedSequenceTranslation, _ := optimizationTable.Translate(optimizedSequence)

	if optimizedSequenceTranslation != strings.ToUpper(gfpTranslation) {
		t.Errorf("TestOptimize has failed. Translate has returned %q, want %q", optimizedSequenceTranslation, gfpTranslation)
	}
}

func TestOptimizeSequence(t *testing.T) {
	t.Parallel()

	var (
		gfpTranslation = "MASKGEELFTGVVPILVELDGDVNGHKFSVSGEGEGDATYGKLTLKFICTTGKLPVPWPTLVTTFSYGVQCFSRYPDHMKRHDFFKSAMPEGYVQERTISFKDDGNYKTRAEVKFEGDTLVNRIELKGIDFKEDGNILGHKLEYNYNSHNVYITADKQKNGIKANFKIRHNIEDGSVQLADHYQQNTPIGDGPVLLPDNHYLSTQSALSKDPNEKRDHMVLLEFVTAAGITHGMDELYK*"
		optimisedGFP   = "ATGGCAAGTAAGGGAGAAGAGCTTTTTACCGGCGTAGTACCAATTCTGGTAGAACTGGATGGTGATGTAAACGGTCACAAATTTAGTGTAAGCGGAGAAGGTGAGGGTGATGCTACCTATGGCAAACTGACCCTAAAGTTTATATGCACGACTGGAAAACTTCCGGTACCGTGGCCAACGTTAGTTACAACGTTTTCTTATGGAGTACAGTGCTTCAGCCGCTACCCAGATCATATGAAACGCCATGATTTCTTTAAGAGCGCCATGCCAGAGGGTTATGTTCAGGAGCGCACGATCTCGTTTAAGGATGATGGTAACTATAAGACTCGTGCTGAGGTGAAGTTCGAAGGCGATACCCTTGTAAATCGTATTGAATTGAAGGGTATAGACTTCAAGGAGGATGGAAATATTCTTGGACATAAGCTGGAATACAATTACAATTCACATAACGTTTATATAACTGCCGACAAGCAAAAAAACGGGATAAAAGCTAATTTTAAAATACGCCACAACATAGAGGACGGGTCGGTGCAACTAGCCGATCATTATCAACAAAACACACCAATCGGCGACGGACCAGTTCTGTTGCCCGATAATCATTACTTATCAACCCAAAGTGCCTTAAGTAAGGATCCGAACGAAAAGCGCGATCATATGGTACTTCTTGAGTTTGTTACCGCTGCAGGCATAACGCATGGCATGGACGAGCTATACAAATAA"
		puc19          = func() genbank.Genbank {
			seq, err := genbank.Read("../../data/puc19.gbk")
			if err != nil {
				t.Fatal(err)
			}

			return seq
		}()
	)

	tests := []struct {
		name string

		sequenceToOptimise string
		updateWeightsWith  genbank.Genbank
		wantOptimised      string

		wantUpdateWeightsErr error
		wantOptimiseErr      error
	}{
		{
			name: "ok",

			sequenceToOptimise: gfpTranslation,
			updateWeightsWith:  puc19,
			wantOptimised:      optimisedGFP,

			wantUpdateWeightsErr: nil,
			wantOptimiseErr:      nil,
		},
		{
			name: "giving no sequence to optimise",

			sequenceToOptimise: "",
			updateWeightsWith:  puc19,
			wantOptimised:      "",

			wantUpdateWeightsErr: nil,
			wantOptimiseErr:      errEmptyAminoAcidString,
		},
		{
			name: "updating weights with a sequence with no CDS",

			sequenceToOptimise: "",
			updateWeightsWith:  genbank.Genbank{},
			wantOptimised:      "",

			wantUpdateWeightsErr: errNoCodingRegions,
			wantOptimiseErr:      errEmptyAminoAcidString,
		},
	}

	for _, tt := range tests {
		var tt = tt
		t.Run(tt.name, func(t *testing.T) {
			t.Parallel()

			optimizationTable, err := NewTranslationTable(11)
			if err != nil {
				t.Fatalf("failed to initialise codon table: %s", err)
			}

			err = optimizationTable.UpdateWeightsWithSequence(tt.updateWeightsWith)
			if !errors.Is(err, tt.wantUpdateWeightsErr) {
				t.Errorf("got %v, want %v", err, tt.wantUpdateWeightsErr)
			}

			got, err := optimizationTable.Optimize(tt.sequenceToOptimise, 1)
			if !errors.Is(err, tt.wantOptimiseErr) {
				t.Errorf("got %v, want %v", err, tt.wantOptimiseErr)
			}

			if !cmp.Equal(got, tt.wantOptimised) {
				t.Errorf("got and tt.wantOptimised didn't match %s", cmp.Diff(got, tt.wantOptimised))
			}
		})
	}
}

func TestNewAminoAcidChooser(t *testing.T) {
	var (
		mockError = errors.New("new chooser rigged to fail")
	)

	tests := []struct {
		name string

		aminoAcids []AminoAcid

		chooserFn func(choices ...weightedRand.Choice) (*weightedRand.Chooser, error)

		wantErr error
	}{
		{
			name: "ok",

			aminoAcids: []AminoAcid{
				{
					Letter: "R",
					Codons: []Codon{
						{
							Triplet: "CGU",
							Weight:  1,
						},
					},
				},
			},

			chooserFn: weightedRand.NewChooser,

			wantErr: nil,
		},
		{
			name: "chooser fn constructor error",

			aminoAcids: []AminoAcid{
				{
					Letter: "R",
					Codons: []Codon{
						{
							Triplet: "CGU",
							Weight:  1,
						},
					},
				},
			},

			chooserFn: func(choices ...weightedRand.Choice) (*weightedRand.Chooser, error) {
				return nil, mockError
			},

			wantErr: mockError,
		},
	}

	for _, tt := range tests {
		t.Run(tt.name, func(t *testing.T) {
			// replace chooser fn with test one
			newChooserFn = tt.chooserFn

			defer func() {
				newChooserFn = weightedRand.NewChooser
			}()

			_, err := newAminoAcidChoosers(tt.aminoAcids)
			if !errors.Is(err, tt.wantErr) {
				t.Errorf("got %v, want %v", err, tt.wantErr)
			}
		})
	}
}

func TestUpdateWeights(t *testing.T) {
	var (
		mockError = errors.New("new chooser rigged to fail")
	)

	tests := []struct {
		name string

		aminoAcids []AminoAcid

		chooserFn func(choices ...weightedRand.Choice) (*weightedRand.Chooser, error)

		wantInitErr error
		wantErr     error
	}{
		{
			name: "ok",

			aminoAcids: []AminoAcid{
				{
					Letter: "R",
					Codons: []Codon{
						{
							Triplet: "CGU",
							Weight:  1,
						},
					},
				},
			},

			chooserFn: weightedRand.NewChooser,

			wantErr: nil,
		},
		{
			name: "chooser fn constructor error",

			aminoAcids: []AminoAcid{
				{
					Letter: "R",
					Codons: []Codon{
						{
							Triplet: "CGU",
							Weight:  1,
						},
					},
				},
			},

			chooserFn: func(choices ...weightedRand.Choice) (*weightedRand.Chooser, error) {
				return nil, mockError
			},

			wantInitErr: mockError,
			wantErr:     mockError,
		},
	}

	for _, tt := range tests {
		t.Run(tt.name, func(t *testing.T) {
			// replace chooser fn with test one
			newChooserFn = tt.chooserFn

			defer func() {
				newChooserFn = weightedRand.NewChooser
			}()

			optimizationTable, err := NewTranslationTable(11)
			if !errors.Is(err, tt.wantInitErr) {
				t.Fatalf("got %v, want %v", err, tt.wantInitErr)
				return
			}

			if tt.wantInitErr != nil {
				return
			}

			err = optimizationTable.UpdateWeights(tt.aminoAcids)
			if !errors.Is(err, tt.wantErr) {
				t.Errorf("got %v, want %v", err, tt.wantErr)
			}
		})
	}
}

func TestCopy(t *testing.T) {
	t.Parallel()

	cmpOptions := []cmp.Option{
		cmpopts.IgnoreUnexported(weightedRand.Chooser{}),
	}

	tests := []struct {
		name string

		wantErr error
	}{
		{
			name: "ok",

			wantErr: nil,
		},
	}

	for _, tt := range tests {
		var tt = tt
		t.Run(tt.name, func(t *testing.T) {
			t.Parallel()

			original, err := NewTranslationTable(11)
			if err != nil {
				t.Fatal(err)
			}

			// perform a deep copy (changing the copy will not change the original)

			deepCopy, err := original.Copy()
			if !errors.Is(err, tt.wantErr) {
				t.Errorf("got %v, want %v", err, tt.wantErr)
			}

			// modify fields

			deepCopy.StartCodons[0] = "🍌"
			deepCopy.StopCodons[0] = "🐗"
			deepCopy.AminoAcids = []AminoAcid{}
			deepCopy.Choosers = map[string]weightedRand.Chooser{}
			deepCopy.Stats = &Stats{}
			deepCopy.TranslationMap = map[string]string{}

			// this compares pointers
			if cmp.Equal(deepCopy, original, cmpOptions...) {
				t.Errorf("deepCopy and original matched, we did not want them to %s", cmp.Diff(deepCopy, original, cmpOptions...))
			}

			// we compare the table's fields

			if cmp.Equal(deepCopy.StartCodonTable, original.StartCodons) {
				t.Errorf("deepCopy and original matched, we did not want them to %s", cmp.Diff(deepCopy.StartCodonTable, original.StartCodons))
			}

			if cmp.Equal(deepCopy.StopCodons, original.StopCodons) {
				t.Errorf("deepCopy and original matched, we did not want them to %s", cmp.Diff(deepCopy.StopCodons, original.StopCodons))
			}

			if cmp.Equal(deepCopy.AminoAcids, original.AminoAcids) {
				t.Errorf("deepCopy and original matched, we did not want them to %s", cmp.Diff(deepCopy.AminoAcids, original.AminoAcids))
			}

			if cmp.Equal(deepCopy.Choosers, original.Choosers) {
				t.Errorf("deepCopy and original matched, we did not want them to %s", cmp.Diff(deepCopy.Choosers, original.Choosers))
			}

			if cmp.Equal(deepCopy.Stats, original.Stats) {
				t.Errorf("deepCopy and original matched, we did not want them to %s", cmp.Diff(deepCopy.Stats, original.Stats))
			}

			if cmp.Equal(deepCopy.TranslationMap, original.TranslationMap) {
				t.Errorf("deepCopy and original matched, we did not want them to %s", cmp.Diff(deepCopy.TranslationMap, original.TranslationMap))
			}
		})
	}
}
