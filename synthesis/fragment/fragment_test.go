package fragment

import (
	"testing"
)

func TestFragment(t *testing.T) {
	gene := "atgaaaaaatttaactggaagaaaatagtcgcgccaattgcaatgctaattattggcttactaggtggtttacttggtgcctttatcctactaacagcagccggggtatcttttaccaatacaacagatactggagtaaaaacggctaagaccgtctacaccaatataacagatacaactaaggctgttaagaaagtacaaaatgccgttgtttctgtcatcaattatcaagaaggttcatcttcagattctctaaatgacctttatggccgtatctttggcggaggggacagttctgattctagccaagaaaattcaaaagattcagatggtctacaggtcgctggtgaaggttctggagtcatctataaaaaagatggcaaagaagcctacatcgtaaccaataaccatgttgtcgatggggctaaaaaacttgaaatcatgctttcggatggttcgaaaattactggtgaacttgttggtaaagacacttactctgacctagcagttgtcaaagtatcttcagataaaataacaactgttgcagaatttgcagactcaaactcccttactgttggtgaaaaagcaattgctatcggtagcccacttggtaccgaatacgccaactcagtaacagaaggaatcgtttctagccttagccgtactataacgatgcaaaacgataatggtgaaactgtatcaacaaacgctatccaaacagatgcagccattaaccctggtaactctggtggtgccctagtcaatattgaaggacaagttatcggtattaattcaagtaaaatttcatcaacgtctgcagtcgctggtagtgctgttgaaggtatggggtttgccattccatcaaacgatgttgttgaaatcatcaatcaattagaaaaagatggtaaagttacacgaccagcactaggaatctcaatagcagatcttaatagcctttctagcagcgcaacttctaaattagatttaccagatgaggtcaaatccggtgttgttgtcggtagtgttcagaaaggtatgccagctgacggtaaacttcaagaatatgatgttatcactgagattgatggtaagaaaatcagctcaaaaactgatattcaaaccaatctttacagccatagtatcggagatactatcaaggtaaccttctatcgtggtaaagataagaaaactgtagatcttaaattaacaaaatctacagaagacatatctgattaa"

	_, _, err := Fragment(gene, 90, 110)
	if err != nil {
		t.Errorf(err.Error())
	}
}

func TestUnfragmentable(t *testing.T) {
	// One should not be able to fragment this
	polyA := "AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA"
	_, _, err := Fragment(polyA, 40, 80)
	if err == nil {
		t.Errorf("polyA should fail to fragment")
	}
}

func TestFragmentSizes(t *testing.T) {
	// This tests if minSize > maxSize
	lacZ := "ATGACCATGATTACGCCAAGCTTGCATGCCTGCAGGTCGACTCTAGAGGATCCCCGGGTACCGAGCTCGAATTCACTGGCCGTCGTTTTACAACGTCGTGACTGGGAAAACCCTGGCGTTACCCAACTTAATCGCCTTGCAGCACATCCCCCTTTCGCCAGCTGGCGTAATAGCGAAGAGGCCCGCACCGATCGCCCTTCCCAACAGTTGCGCAGCCTGAATGGCGAATGGCGCCTGATGCGGTATTTTCTCCTTACGCATCTGTGCGGTATTTCACACCGCATATGGTGCACTCTCAGTACAATCTGCTCTGATGCCGCATAG"
	_, _, err := Fragment(lacZ, 105, 95)
	if err == nil {
		t.Errorf("Fragment should fail when minFragmentSize > maxFragmentSize")
	}

	_, _, err = Fragment(lacZ, 7, 95)
	if err == nil {
		t.Errorf("Fragment should fail when minFragmentSize < 8")
	}

}

func TestSmallFragmentSize(t *testing.T) {
	// The following should succeed, but will require setting minFragmentSize = 12
	lacZ := "ATGACCATGATTACGCCAAGCTTGCATGCCTGCAGGTCGACTCTAGAGGATCCCCGGGTACCGAGCTCGAATTCACTGGCCGTCGTTTTACAACGTCGTGACTGGGAAAACCCTGGCGTTACCCAACTTAATCGCCTTGCAGCACATCCCCCTTTCGCCAGCTGGCGTAATAGCGAAGAGGCCCGCACCGATCGCCCTTCCCAACAGTTGCGCAGCCTGAATGGCGAATGGCGCCTGATGCGGTATTTTCTCCTTACGCATCTGTGCGGTATTTCACACCGCATATGGTGCACTCTCAGTACAATCTGCTCTGATGCCGCATAG"
	_, _, err := Fragment(lacZ, 12, 20)
	if err != nil {
		t.Errorf("Got error in small fragmentation: %s", err)
	}
}

func TestLongFragment(t *testing.T) {
	// A regression test for a bug that sometimes fragmented a sequence to
	// be longer than its max length
	gene := "GGAGGGTCTCAATGCTGGACGATCGCAAATTCAGCGAACAGGAGCTGGTCCGTCGCAACAAATACAAAACGCTGGTCGAGCAAAACAAAGACCCGTACAAGATTACGAACTGGAAACGCAATACCACCCTGCTGAAACTGAATGAGAAATACAAAGACTATAGCAAGGAGGACCTGTTGAACCTGAATCAAGAACTGGTCGTTGTTGCAGGTCGTATCAAACTGTATCGTGAAGCCGGTAAAAAAGCTGCCTTTGTGAACATTGATGATCAAGACTCCTCTATTCAGTTGTACGTGCGCCTGGATGAGATCGGTGATCAGAGCTTCGAGGATTTCCGCAATTTCGACCTGGGTGACATCATTGGTGTTAAAGGTATCATGATGCGCACCGACCACGGCGAGTTGAGCATCCGTTGTAAGGAAGTCGTGCTGCTGAGCAAGGCCCTGCGTCCGCTGCCGGATAAACACGCGGGCATTCAGGATATTGAGGAAAAGTACCGCCGTCGCTATGTGGACCTGATTATGAATCACGACGTGCGCAAGACGTTCCAGGCGCGTACCAAGATCATTCGTACCTTGCAAAACTTTCTGGATAATAAGGGTTACATGGAGGTCGAAACCCCGATCCTGCATAGCCTGAAGGGTGGCGCGAGCGCGAAACCGTTTATTACCCACTACAATGTGCTGAATACGGATGTGTATCTGCGTATCGCGACCGAGCTGCACCTGAAACGCCTGATTGTTGGCGGTTTCGAGGGTGTGTATGAGATCGGTCGCATCTTTCGCAATGAAGGTATGTCCACGCGTCACAATCCGGAATTCACGTCTATCGAACTGTATGTCGCCTATGAGGACATGTTCTTTTTGATGGATCTGACCGAAGAGATTTTTCGCGTTTGTAATGCCGCAGTCAACAGCTCCAGCATCATTGAGTATAACAACGTGAAAATTGACCTGAGCAAGCCGTTTAAGCGCCTGCATATGGTTGACGGTATTAAACAGGTGACCGGCGTCGACTTCTGGCAGGAGATGACGGTCCAACAGGCTCTGGAGCTGGCCAAAAAGCATAAAGTGCACGTTGAAAAACATCAAGAGTCTGTTGGTCACATTATCAATTTGTTCTATGAGGAGTTCGTGGAGTCCACGATTGTTGAGCCGACGTTCGTGTACGGTCACCCGAAGGAAATCTCTCCGCTGGCTAAGAGCAATCCGTCTGACCCGCGTTTCACGGACCGTTTCGAGCTGTTCATTCTGGGTCGTGAGTATGCGAATGCGTTTAGCGAGCTGAATGACCCGATTGACCAGTACGAACGCTTCAAGGCTCAGATTGAGGAGGAAAGCAAGGGCAACGATGAAGCCAACGACATGGACATTGATTTCATCGAGGCTCTGGAACACGCCATGCCGCCGACCGCGGGTATTGGTATCGGCATTGATCGCTTGGTTATGCTGCTGACGAATAGCGAATCCATCAAAGACGTGCTGTTGTTCCCGCAAATGAAGCCGCGCGAATGAAGAGCTTAGAGACCCGCT"
	frags, _, err := Fragment(gene, 79, 94)
	if err != nil {
		t.Errorf(err.Error())
	}
	for _, frag := range frags {
		if len(frag) > 94 {
			t.Errorf("Fragment too long. Expected <94, got %d", len(frag))
		}
	}
}

func TestCheckLongRegresssion(t *testing.T) {
	overhangs := []string{"AGAC"}
	newOverhangs, _ := NextOverhangs(overhangs)
	foundGTCT := false
	for _, overhang := range newOverhangs {
		if overhang == "GTCT" {
			foundGTCT = true
		}
	}
	if foundGTCT {
		t.Errorf("Should not have found GTCT since it is the reverse complement of AGAC")
	}

}

func TestRegressionTestMatching12(t *testing.T) {
	overhangs := []string{"CGAG", "GTCT", "TACT", "AATG", "ATCC", "CGCT", "AAAA", "AAGT", "ATAG", "ATTA", "AGTC", "AGGA", "ACAA", "ACGC", "TAGA", "TACA", "TTAC", "TTCG", "TGAG", "TCGG", "GAAG", "GTGC", "GCCG", "CAGG", "TATC", "GGGG", "AGCA", "GTTG", "ATGA", "TAAA", "TCTC", "TGAC", "ACCA", "CACG", "TACC", "TTGC", "TCAG", "AGAA", "GAAC", "AACT", "TTTG", "ATAA", "AAGG", "AATA", "TAGC", "TGGA", "TGCG", "TAAG"}
	efficiency := SetEfficiency(overhangs)
	if efficiency < 0.12 || efficiency > 0.13 {
		t.Errorf("Expected efficiency of 0.129 - approximately matches NEB ligase fidelity viewer")
	}
}
