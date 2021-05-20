package poly

import (
	"fmt"
	"math"
	"strings"
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

func ExampleDeBruijn() {
	a := DeBruijn(4)

	fmt.Println(a)
	// Output: AAAATAAAGAAACAATTAATGAATCAAGTAAGGAAGCAACTAACGAACCATATAGATACATTTATTGATTCATGTATGGATGCATCTATCGATCCAGAGACAGTTAGTGAGTCAGGTAGGGAGGCAGCTAGCGAGCCACACTTACTGACTCACGTACGGACGCACCTACCGACCCTTTTGTTTCTTGGTTGCTTCGTTCCTGTGTCTGGGTGGCTGCGTGCCTCTCGGTCGCTCCGTCCCGGGGCGGCCGCGCCCCAAA
}

func ExampleUniqueSequence() {
	uniqueSequences := make(chan string)
	go UniqueSequence(uniqueSequences, 20, 4, []string{"CTCTCGGTCGCTCC"}, []func(string) bool{})
	var output string
	for uniqueSeq := range uniqueSequences {
		output = uniqueSeq
	}

	// Only we only want 1 output to check against
	fmt.Println(output)
	// Output: TCTCGGTCGCTCCGTCCCGG
}

func TestUniqueSequence(t *testing.T) {
	c := make(chan string)
	testFunc := func(s string) bool {
		if strings.Contains(s, "GGCCGCGCCCC") {
			return false
		} else {
			return true
		}
	}
	go UniqueSequence(c, 20, 4, []string{}, []func(string) bool{testFunc})
	var output string
	for a := range c {
		output = a
	}
	if output != "CTCTCGGTCGCTCCGTCCCG" {
		t.Errorf("TestUniqueSequence function should return CTCTCGGTCGCTCCGTCCCG. Got:\n%s", output)
	}

	b := make(chan string)
	go UniqueSequence(b, 20, 4, []string{"GGCCGCGCCCC"}, []func(string) bool{})
	for d := range b {
		output = d
	}
	if output != "CTCTCGGTCGCTCCGTCCCG" {
		t.Errorf("TestUniqueSequence string should return CTCTCGGTCGCTCCGTCCCG. Got:\n%s", output)
	}
}
