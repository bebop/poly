package poly

import (
	"fmt"
	"math"
	"testing"
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

	testReverseComplement := ReverseComplement(testSeq)
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

func ExampleMakePrimers() {
	sequenceString := "ATTAAGCATTCTGCCGACATGGAAGCCATCACAAACGGCATGATGAACCTGAATCGCCAGCGGCATCAGCACCTTGTCGCCTTGCGTATAATATTTGCCCATGGTGAAAACGGGGGCGAAGAAGTTGTCCATATTGGCCACGTTTAAATCAAAACTGGTGAAACTCACCCAGGGATTGGCTGAGACGAAAAACATATTCTCAATAAACCCTTTAGGGAAATAGGCCAGGTTTTCACCGTAACACGCCACATCTTGCGAATATATGTGTAGAAACTGCCGGAAATCGTCGTGGTATTCACTCCAGAGCGATGAAAACGTTTCAGTTTGCTCATGGAAAACGGTGTAACAAGGGTGAACACTATCCCATATCACCAGCTCACCGTCTTTCATTGCCATACGAAATTCCGGATGAGCATTCATCAGGCGGGCAAGAATGTGAATAAAGGCCGGATAAAACTTGTGCTTATTTTTCTTTACGGTCTTTAAAAAGGCC" // some random sequence from camR
	targetTm := 55.0
	forwardPrimer, reversePrimer := MakePrimers(sequenceString, targetTm)

	fmt.Println(forwardPrimer, " ", reversePrimer)
	// Output: ATTAAGCATTCTGCCGACATGG   GGCCTTTTTAAAGACCGTAAAGAAA
}

func ExamplePcr() {
	sequence := "AAATGCCGCAAAAAAGGGAATAAGGGCGACACGGGTCATAGCTGTTTCCTGGTAAAACGACGGCCAGTATGCAGTTTAAGGTTTACACCTATAAAAGAGAGAGCCGTTATCGTCTGTTTGTGGATGTACAGAGTGATATTATTGACACGCCCGGGCGACGGATGGTGATCCCCCTGGCCAGTGCACGTCTGCTGTCAGATAAAGTCTCCCGTGAACTTTACCCGGTGGTGCATATCGGGGATGAAAGCTGGCGCATGATGACCACCGATATGGCCAGTGTGCCGGTCTCCGTTATCGGGGAAGAAGTGGCTGATCTCAGCCACCGCGAAAATGACATCAAAAACGCCATTAACCTGATGTTCTGGGGAATATAAGTCATAGCTGTTTCCTGGTAAAACGACGGCCAGTAGTCAAAAGCCTCCGACCGGAGGCTTTTGACTGTCATAGCTGTTTCCTGGTAAAACGACGGCCAGTGTAAAACGACGGCCAGTACTCTTCCTTTTTCAATATTATTGAAGCATTTATCAGGGTTATTGTCTCATGAGCGGATACATATTTGAATGTATTTAGAAAAATAAACAAATAGGGGTTCCGCGGTCATAGCTGTTTCCTGGTCATAGCTGTTTCCTGGTAAAACGACGGCCAGT"
	m13for := "GTAAAACGACGGCCAGT"
	m13rev := "CAGGAAACAGCTATGAC"

	fragments := Pcr(sequence, true, 45, m13for, m13rev)

	fmt.Println(len(fragments))
	// Output: 4
}

func TestPcr(t *testing.T) {
	sequence := "AAATGCCGCAAAAAAGGGAATAAGGGCGACACGGGTCATAGCTGTTTCCTGGTAAAACGACGGCCAGTATGCAGTTTAAGGTTTACACCTATAAAAGAGAGAGCCGTTATCGTCTGTTTGTGGATGTACAGAGTGATATTATTGACACGCCCGGGCGACGGATGGTGATCCCCCTGGCCAGTGCACGTCTGCTGTCAGATAAAGTCTCCCGTGAACTTTACCCGGTGGTGCATATCGGGGATGAAAGCTGGCGCATGATGACCACCGATATGGCCAGTGTGCCGGTCTCCGTTATCGGGGAAGAAGTGGCTGATCTCAGCCACCGCGAAAATGACATCAAAAACGCCATTAACCTGATGTTCTGGGGAATATAAGTCATAGCTGTTTCCTGGTAAAACGACGGCCAGTAGTCAAAAGCCTCCGACCGGAGGCTTTTGACTGTCATAGCTGTTTCCTGGTAAAACGACGGCCAGTGTAAAACGACGGCCAGTACTCTTCCTTTTTCAATATTATTGAAGCATTTATCAGGGTTATTGTCTCATGAGCGGATACATATTTGAATGTATTTAGAAAAATAAACAAATAGGGGTTCCGCGGTCATAGCTGTTTCCTGGTCATAGCTGTTTCCTGGTAAAACGACGGCCAGT"
	m13for := "GTAAAACGACGGCCAGT"
	m13rev := "CAGGAAACAGCTATGAC"

	// Test example PCR
	fragments := Pcr(sequence, true, 45, m13for, m13rev)
	if len(fragments) != 4 {
		t.Errorf("Should have 4 output fragments with circular PCR")
	}

	// Test example linear PCR
	fragmentsLinear := Pcr(sequence+"GTCATAGCTGTTTCCTG", false, 45, m13for, m13rev)
	if len(fragmentsLinear) != 4 {
		t.Errorf("Should have 4 output fragments with linear PCR (though these will be different than circular PCR)")
	}

	// This should make reverse primer fail due to too high tm requirements
	fragments = Pcr(sequence, true, 50, m13for, m13rev)
	if len(fragments) != 0 {
		t.Errorf("Tm calculated is too high - m13rev should be <50")
	}

	// This should make forward primer fail due to too high tm requirements
	fragments = Pcr(sequence, true, 60, m13for, m13rev)
	if len(fragments) != 0 {
		t.Errorf("Tm calculated is too high - m13for should be <60")
	}

}
