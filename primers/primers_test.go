package primers

import (
	"fmt"
	"math"
	"strings"
	"testing"

	"github.com/TimothyStiles/poly/checks"
	"github.com/TimothyStiles/poly/transform"
)

func ExampleMarmurDoty() {
	sequenceString := "ACGTCCGGACTT"
	meltingTemp := MarmurDoty(sequenceString)

	fmt.Println(meltingTemp)
	// output: 31
}

func TestMarmurDoty(t *testing.T) {
	testSeq := "ACGTCCGGACTT"
	expectedTM := 31.0
	if calcTM := MarmurDoty(testSeq); expectedTM != calcTM {
		t.Errorf("MarmurDoty has changed on test. Got %f instead of %f", calcTM, expectedTM)
	}
}

func ExampleSantaLucia() {

	sequenceString := "ACGATGGCAGTAGCATGC" //"GTAAAACGACGGCCAGT" // M13 fwd
	testCPrimer := 0.1e-6                  // primer concentration
	testCNa := 350e-3                      // salt concentration
	testCMg := 0.0                         // magnesium concentration
	expectedTM := 62.7                     // roughly what we're expecting with a margin of error
	meltingTemp, _, _ := SantaLucia(sequenceString, testCPrimer, testCNa, testCMg)
	withinMargin := math.Abs(expectedTM-meltingTemp)/expectedTM >= 0.02 // checking margin of error

	fmt.Println(withinMargin)
	// output: false
}
func TestSantaLucia(t *testing.T) {
	testSeq := "ACGATGGCAGTAGCATGC" //"GTAAAACGACGGCCAGT" // M13 fwd
	testCPrimer := 0.1e-6
	testCNa := 350e-3
	testCMg := 0.0
	expectedTM := 62.7
	if calcTM, _, _ := SantaLucia(testSeq, testCPrimer, testCNa, testCMg); math.Abs(expectedTM-calcTM)/expectedTM >= 0.02 {
		t.Errorf("SantaLucia has changed on test. Got %f instead of %f", calcTM, expectedTM)
	}
}

func TestSantaLuciaReverseComplement(t *testing.T) {
	testSeq := "ACGTAGATCTACGT" //"GTAAAACGACGGCCAGT" // M13 fwd

	testReverseComplement := transform.ReverseComplement(testSeq)
	if testSeq != testReverseComplement {
		t.Errorf("Input is not a reverse complement of it's. Got %q instead of %q", testSeq, testReverseComplement)
	}
	testCPrimer := 0.1e-6
	testCNa := 350e-3
	testCMg := 0.0
	expectedTM := 47.428514
	if calcTM, _, _ := SantaLucia(testSeq, testCPrimer, testCNa, testCMg); math.Abs(expectedTM-calcTM)/expectedTM >= 0.02 {
		t.Errorf("SantaLucia has changed on test. Got %f instead of %f", calcTM, expectedTM)
	}

}

func ExampleMeltingTemp() {
	sequenceString := "GTAAAACGACGGCCAGT" // M13 fwd
	expectedTM := 52.8
	meltingTemp := MeltingTemp(sequenceString)
	withinMargin := math.Abs(expectedTM-meltingTemp)/expectedTM >= 0.02

	fmt.Println(withinMargin)
	// output: false
}

func TestMeltingTemp(t *testing.T) {
	testSeq := "GTAAAACGACGGCCAGT" // M13 fwd
	expectedTM := 52.8
	if calcTM := MeltingTemp(testSeq); math.Abs(expectedTM-calcTM)/expectedTM >= 0.02 {
		t.Errorf("MeltingTemp has changed on test. Got %f instead of %f", calcTM, expectedTM)
	}
}

func ExampleNucleobaseDeBruijnSequence() {
	a := NucleobaseDeBruijnSequence(4)

	fmt.Println(a)
	// Output: AAAATAAAGAAACAATTAATGAATCAAGTAAGGAAGCAACTAACGAACCATATAGATACATTTATTGATTCATGTATGGATGCATCTATCGATCCAGAGACAGTTAGTGAGTCAGGTAGGGAGGCAGCTAGCGAGCCACACTTACTGACTCACGTACGGACGCACCTACCGACCCTTTTGTTTCTTGGTTGCTTCGTTCCTGTGTCTGGGTGGCTGCGTGCCTCTCGGTCGCTCCGTCCCGGGGCGGCCGCGCCCCAAA
}

func ExampleCreateBarcodesWithBannedSequences() {
	barcodes := CreateBarcodesWithBannedSequences(20, 4, []string{"CTCTCGGTCGCTCC"}, []func(string) bool{})

	fmt.Println(barcodes[0])
	// Output: AAAATAAAGAAACAATTAAT
}

func ExampleCreateBarcodes() {
	barcodes := CreateBarcodes(20, 4)

	fmt.Println(barcodes[0])
	// Output: AAAATAAAGAAACAATTAAT
}

func TestCreateBarcode(t *testing.T) {
	testFunc := func(s string) bool {
		return !strings.Contains(s, "GGCCGCGCCCC")
	}
	barcodes := CreateBarcodesWithBannedSequences(20, 4, []string{}, []func(string) bool{testFunc})
	output := barcodes[len(barcodes)-1]
	if output != "CTCTCGGTCGCTCCGTCCCG" {
		t.Errorf("TestUniqueSequence function should return CTCTCGGTCGCTCCGTCCCG. Got:\n%s", output)
	}

	barcodes = CreateBarcodesWithBannedSequences(20, 4, []string{"GGCCGCGCCCC"}, []func(string) bool{})
	output = barcodes[len(barcodes)-1]
	if output != "CTCTCGGTCGCTCCGTCCCG" {
		t.Errorf("TestUniqueSequence string should return CTCTCGGTCGCTCCGTCCCG. Got:\n%s", output)
	}

	barcodes = CreateBarcodesWithBannedSequences(20, 4, []string{transform.ReverseComplement("GGCCGCGCCCC")}, []func(string) bool{})
	output = barcodes[len(barcodes)-1]
	if output != "CTCTCGGTCGCTCCGTCCCG" {
		t.Errorf("TestUniqueSequence string should return CTCTCGGTCGCTCCGTCCCG. Got:\n%s", output)
	}
}

// Example of creating a large number of DNA Primers with the following properties:
// x consistent GC content
// x roughly the same melting temperature
// - do not dimerize with themselves or other primers in the set
// x don't share any subsequences >4bp with any other primer in the set
// - don't bind to any DNA sequences from a list of bkgnd sequences (ie FreeGenes, and genomes of E.coli,B.subtilis,S.cerevisiae,P.pastoris)
func ExamplePrimerWorkflow() {

	WithinMeltingTempRange := func(dna string, target_temp float64, margin float64) bool {
		if math.Abs(MeltingTemp(dna)-target_temp) > margin {
			return false
		} else {
			return true
		}
	}
	WithinGCRange := func(dna string, target_GC float64, margin float64) bool {
		if math.Abs(checks.GcContent(dna)-target_GC) > margin {
			return false
		} else {
			return true
		}
	}

	target_temp := 40.0
	temp_margin := 1.0
	IsBannedMeltingTemperature := func(dna string) bool {
		return !WithinMeltingTempRange(dna, target_temp, temp_margin)
	}

	target_GC := 40.0
	GC_margin := 1.0
	IsBannedGcContent := func(dna string) bool {
		return !WithinGCRange(dna, target_GC, GC_margin)
	}

	banned_functions := []func(string) bool{
		IsBannedMeltingTemperature, IsBannedGcContent,
	}
	barcode_length := 10
	maxSubSequence := 5
	bannedSequences := []string{}
	generated_primers := CreateBarcodesWithBannedSequences(
		barcode_length, maxSubSequence, bannedSequences, banned_functions,
	)

	fmt.Println(generated_primers[0])
}
