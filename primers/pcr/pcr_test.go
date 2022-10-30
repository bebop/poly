package pcr

import (
	"fmt"
	"strings"
	"testing"

	"github.com/TimothyStiles/poly/transform"
)

// gene is a gene for testing PCR
var gene string = "aataattacaccgagataacacatcatggataaaccgatactcaaagattctatgaagctatttgaggcacttggtacgatcaagtcgcgctcaatgtttggtggcttcggacttttcgctgatgaaacgatgtttgcactggttgtgaatgatcaacttcacatacgagcagaccagcaaacttcatctaacttcgagaagcaagggctaaaaccgtacgtttataaaaagcgtggttttccagtcgttactaagtactacgcgatttccgacgacttgtgggaatccagtgaacgcttgatagaagtagcgaagaagtcgttagaacaagccaatttggaaaaaaagcaacaggcaagtagtaagcccgacaggttgaaagacctgcctaacttacgactagcgactgaacgaatgcttaagaaagctggtataaaatcagttgaacaacttgaagagaaaggtgcattgaatgcttacaaagcgatacgtgactctcactccgcaaaagtaagtattgagctactctgggctttagaaggagcgataaacggcacgcactggagcgtcgttcctcaatctcgcagagaagagctggaaaatgcgctttcttaa"
var geneUpperCaseComplement string = "TTATTAATGTGGCTCTATTGTGTAGTACCTATTTGGCTATGAGTTTCTAAGATACTTCGATAAACTCCGTGAACCATGCTAGTTCAGCGCGAGTTACAAACCACCGAAGCCTGAAAAGCGACTACTTTGCTACAAACGTGACCAACACTTACTAGTTGAAGTGTATGCTCGTCTGGTCGTTTGAAGTAGATTGAAGCTCTTCGTTCCCGATTTTGGCATGCAAATATTTTTCGCACCAAAAGGTCAGCAATGATTCATGATGCGCTAAAGGCTGCTGAACACCCTTAGGTCACTTGCGAACTATCTTCATCGCTTCTTCAGCAATCTTGTTCGGTTAAACCTTTTTTTCGTTGTCCGTTCATCATTCGGGCTGTCCAACTTTCTGGACGGATTGAATGCTGATCGCTGACTTGCTTACGAATTCTTTCGACCATATTTTAGTCAACTTGTTGAACTTCTCTTTCCACGTAACTTACGAATGTTTCGCTATGCACTGAGAGTGAGGCGTTTTCATTCATAACTCGATGAGACCCGAAATCTTCCTCGCTATTTGCCGTGCGTGACCTCGCAGCAAGGAGTTAGAGCGTCTCTTCTCGACCTTTTACGCGAAAGAATT"

func TestRegression(t *testing.T) {
	uppercaseGeneComplement := strings.ToUpper(transform.Complement(gene))
	expectedFirstFragment := "TTATAGGTCTCATACTAATAATTACACCGAGATAACACATCATGGATAAACCGATACTCAAAGATTCTATGAAGCTATTTGAGGCACTTGGTACGATCAAGTCGCGCTCAATGTTTGGTGGCTTCGGACTTTTCGCTGATGAAACGATGTTTGCACTGGTTGTGAATGATCAACTTCACATACGAGCAGACCAGCAAACTTCATCTAACTTCGAGAAGCAAGGGCTAAAACCGTACGTTTATAAAAAGCGTGGTTTTCCAGTCGTTACTAAGTACTACGCGATTTCCGACGACTTGTGGGAATCCAGTGAACGCTTGATAGAAGTAGCGAAGAAGTCGTTAGAACAAGCCAATTTGGAAAAAAAGCAACAGGCAAGTAGTAAGCCCGACAGGTTGAAAGACCTGCCTAACTTACGACTAGCGACTGAACGAATGCTTAAGAAAGCTGGTATAAAATCAGTTGAACAACTTGAAGAGAAAGGTGCATTGAATGCTTACAAAGCGATACGTGACTCTCACTCCGCAAAAGTAAGTATTGAGCTACTCTGGGCTTTAGAAGGAGCGATAAACGGCACGCACTGGAGCGTCGTTCCTCAATCTCGCAGAGAAGAGCTGGAAAATGCGCTTTCTTAAATGAAGAGACCATATA"
	fmt.Println(uppercaseGeneComplement)
	// These primers are arbitrary, but they should work for this test.
	testPrimers := []string{"TATATGGTCTCTTCATTTAAGAAAGCGCATTTTCCAGC", "TTATAGGTCTCATACTAATAATTACACCGAGATAACACATCATGG", "CTGCAGGTCGACTCTAG"}

	fragments, err := Simulate([]string{gene}, 55.0, false, testPrimers)

	if err != nil {
		t.Errorf("%s", err)
	}

	if expectedFirstFragment != fragments[0] {
		t.Errorf("Expected %s, got %s", expectedFirstFragment, fragments[0])
	}

	// Now let's make sure it only amplified our target.
	fmt.Println(fragments[0])

}

func TestSimulatePrimerRejection(t *testing.T) {
	// CTGCAGGTCGACTCTAG is too low tm for this function, so there is a break in the logic.
	primers := []string{"TATATGGTCTCTTCATTTAAGAAAGCGCATTTTCCAGC", "TTATAGGTCTCATACTAATAATTACACCGAGATAACACATCATGG", "CTGCAGGTCGACTCTAG"}
	fragments, _ := Simulate([]string{gene}, 55.0, false, primers)
	if len(fragments) != 1 {
		t.Errorf("Should only have one fragment")
	}
}

func TestSimulateMoreThanOneForward(t *testing.T) {
	// This tests the first bit of logic in simulate.
	// If this primer isn't last forward binding primer AND there is
	// another reverse primer binding site.

	// gatactcaaagattctatgaagctatttgaggcacttggtacg occurs internally inside of
	// gene
	internalPrimer := "gatactcaaagattctatgaagctatttgaggcacttggtacg"

	// reversePrimer is a different primer from normal that will bind inside
	// of gene.
	reversePrimer := "tatcgctttgtaagcattcaatgcacctttctcttcaagttg"

	// outsideForwardPrimer is a primer that binds out of the range of
	// reversePrimer
	outsideForwardPrimer := "gtcgttcctcaatctcgcagagaagagctggaaaatg"

	primers := []string{internalPrimer, reversePrimer, outsideForwardPrimer}
	fragments, _ := Simulate([]string{gene}, 55.0, false, primers)
	if len(fragments) != 1 {
		t.Errorf("Should only have one fragment")
	}
}

func TestSimulateCircular(t *testing.T) {
	// This tests for circular simulations.

	// forwardPrimer binds near to the end of gene
	forwardPrimer := "actctgggctttagaaggagcgataaacggc"
	// reversePrimer binds to the beginning of gene, in the opposite direction
	reversePrimer := "aagtgcctcaaatagcttcatagaatctttgagtatcgg"

	// targetFragment is what the amplification reaction should result in
	targetFragment := "ACTCTGGGCTTTAGAAGGAGCGATAAACGGCACGCACTGGAGCGTCGTTCCTCAATCTCGCAGAGAAGAGCTGGAAAATGCGCTTTCTTAAAATAATTACACCGAGATAACACATCATGGATAAACCGATACTCAAAGATTCTATGAAGCTATTTGAGGCACTT"

	primers := []string{forwardPrimer, reversePrimer}
	fragments, _ := Simulate([]string{gene}, 55.0, true, primers)
	if fragments[0] != targetFragment {
		t.Errorf("Didn't get target fragment from circular pcr. Expected: %s, got: %s", targetFragment, fragments[0])
	}
}

func TestSimulateConcatemerization(t *testing.T) {
	// This tests the concatermization detector

	forwardPrimer := "AATAATTACACCGAGATAACACATCATGG"

	// This reverse primer will add in a forward primer binding site,
	// allowing concatemerization
	reversePrimer := "CCATGATGTGTTATCTCGGTGTAATTATTTTAAGAAAGCGCATTTTCCAGC"

	_, err := Simulate([]string{gene}, 55.0, false, []string{forwardPrimer, reversePrimer})
	if err == nil {
		t.Errorf("Should have gotten concatemerization")
	}
}
