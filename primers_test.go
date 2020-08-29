package main

import (
	"math"
	"testing"
)

func TestMarmurDoty(t *testing.T) {
	testSeq := "ACGTCCGGACTT"
	expectedTM := 31.0
	if calcTM := MarmurDoty(testSeq); expectedTM != calcTM {
		t.Errorf("MarmurDoty has changed on test. Got %f instead of %f", calcTM, expectedTM)
	}
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

func TestMeltingTemp(t *testing.T) {
	testSeq := "GTAAAACGACGGCCAGT" // M13 fwd
	expectedTM := 52.8
	if calcTM := MeltingTemp(testSeq); math.Abs(expectedTM-calcTM)/expectedTM >= 0.02 {
		t.Errorf("MeltingTemp has changed on test. Got %f instead of %f", calcTM, expectedTM)
	}
}

func TestCheckHomopolymericRuns(t *testing.T) {
	testSeq := "GTAAAACGACGGCCAGT" // M13 fwd
	if run := checkHomopolymericRuns(testSeq); run == false {
		t.Errorf("checkHomopolymericRuns has changed on test. Got false instead of true")
	}

	testSeq = "ACGATGGCAGTAGCATGC" //"GTAAAACGACGGCCAGT" // M13 fwd

	if run := checkHomopolymericRuns(testSeq); run == true {
		t.Errorf("checkHomopolymericRuns has changed on test. Got true instead of false")
	}

}

func TestGCPercentage(t *testing.T) {
	testSeq := "GcccGcGG"
	expectedPercentage := 1.000000
	if percentage := GCPercentage(testSeq); percentage != expectedPercentage {
		t.Errorf("GCPercentage has changed on test. Got %f instead of %f", percentage, expectedPercentage)
	}
}

func TestCheckRepeats(t *testing.T) {
	testSeq := "AATGGCATGCGCCATT"
	expectedValue := true
	if flag := checkRepeats(testSeq); flag != expectedValue {
		t.Errorf("checkRepeats has changed on reverse complement repeat test for possible primer hairpins. Got %t instead of %t", flag, expectedValue)
	}

	testSeq = "AATCGAAATCGA"
	expectedValue = true
	if flag := checkRepeats(testSeq); flag != expectedValue {
		t.Errorf("checkRepeats has changed on repeat test. Got %t instead of %t", flag, expectedValue)
	}

}

func TestGCClamp(t *testing.T) {
	testSeq := "AATCGAAATCGG"
	expectedValue := true
	if flag := GCClamp(testSeq); flag != expectedValue {
		t.Errorf("GCClamp has changed on repeat test. Got %t instead of %t", flag, expectedValue)
	}

	testSeq = "AATCGAAATCGC"
	expectedValue = true
	if flag := GCClamp(testSeq); flag != expectedValue {
		t.Errorf("GCClamp has changed on repeat test. Got %t instead of %t", flag, expectedValue)
	}

}
