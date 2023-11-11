package clone

import (
	"testing"
)

// pOpen plasmid series (https://stanford.freegenes.org/collections/open-genes/products/open-plasmids#description). I use it for essentially all my cloning. -Keoni
var popen = Part{"TAACTATCGTCTTGAGTCCAACCCGGTAAGACACGACTTATCGCCACTGGCAGCAGCCACTGGTAACAGGATTAGCAGAGCGAGGTATGTAGGCGGTGCTACAGAGTTCTTGAAGTGGTGGCCTAACTACGGCTACACTAGAAGAACAGTATTTGGTATCTGCGCTCTGCTGAAGCCAGTTACCTTCGGAAAAAGAGTTGGTAGCTCTTGATCCGGCAAACAAACCACCGCTGGTAGCGGTGGTTTTTTTGTTTGCAAGCAGCAGATTACGCGCAGAAAAAAAGGATCTCAAGAAGGCCTACTATTAGCAACAACGATCCTTTGATCTTTTCTACGGGGTCTGACGCTCAGTGGAACGAAAACTCACGTTAAGGGATTTTGGTCATGAGATTATCAAAAAGGATCTTCACCTAGATCCTTTTAAATTAAAAATGAAGTTTTAAATCAATCTAAAGTATATATGAGTAAACTTGGTCTGACAGTTACCAATGCTTAATCAGTGAGGCACCTATCTCAGCGATCTGTCTATTTCGTTCATCCATAGTTGCCTGACTCCCCGTCGTGTAGATAACTACGATACGGGAGGGCTTACCATCTGGCCCCAGTGCTGCAATGATACCGCGAGAACCACGCTCACCGGCTCCAGATTTATCAGCAATAAACCAGCCAGCCGGAAGGGCCGAGCGCAGAAGTGGTCCTGCAACTTTATCCGCCTCCATCCAGTCTATTAATTGTTGCCGGGAAGCTAGAGTAAGTAGTTCGCCAGTTAATAGTTTGCGCAACGTTGTTGCCATTGCTACAGGCATCGTGGTGTCACGCTCGTCGTTTGGTATGGCTTCATTCAGCTCCGGTTCCCAACGATCAAGGCGAGTTACATGATCCCCCATGTTGTGCAAAAAAGCGGTTAGCTCCTTCGGTCCTCCGATCGTTGTCAGAAGTAAGTTGGCCGCAGTGTTATCACTCATGGTTATGGCAGCACTGCATAATTCTCTTACTGTCATGCCATCCGTAAGATGCTTTTCTGTGACTGGTGAGTACTCAACCAAGTCATTCTGAGAATAGTGTATGCGGCGACCGAGTTGCTCTTGCCCGGCGTCAATACGGGATAATACCGCGCCACATAGCAGAACTTTAAAAGTGCTCATCATTGGAAAACGTTCTTCGGGGCGAAAACTCTCAAGGATCTTACCGCTGTTGAGATCCAGTTCGATGTAACCCACTCGTGCACCCAACTGATCTTCAGCATCTTTTACTTTCACCAGCGTTTCTGGGTGAGCAAAAACAGGAAGGCAAAATGCCGCAAAAAAGGGAATAAGGGCGACACGGAAATGTTGAATACTCATACTCTTCCTTTTTCAATATTATTGAAGCATTTATCAGGGTTATTGTCTCATGAGCGGATACATATTTGAATGTATTTAGAAAAATAAACAAATAGGGGTTCCGCGCACCTGCACCAGTCAGTAAAACGACGGCCAGTAGTCAAAAGCCTCCGACCGGAGGCTTTTGACTTGGTTCAGGTGGAGTGGGAGTAgtcttcGCcatcgCtACTAAAagccagataacagtatgcgtatttgcgcgctgatttttgcggtataagaatatatactgatatgtatacccgaagtatgtcaaaaagaggtatgctatgaagcagcgtattacagtgacagttgacagcgacagctatcagttgctcaaggcatatatgatgtcaatatctccggtctggtaagcacaaccatgcagaatgaagcccgtcgtctgcgtgccgaacgctggaaagcggaaaatcaggaagggatggctgaggtcgcccggtttattgaaatgaacggctcttttgctgacgagaacagggGCTGGTGAAATGCAGTTTAAGGTTTACACCTATAAAAGAGAGAGCCGTTATCGTCTGTTTGTGGATGTACAGAGTGATATTATTGACACGCCCGGGCGACGGATGGTGATCCCCCTGGCCAGTGCACGTCTGCTGTCAGATAAAGTCTCCCGTGAACTTTACCCGGTGGTGCATATCGGGGATGAAAGCTGGCGCATGATGACCACCGATATGGCCAGTGTGCCGGTCTCCGTTATCGGGGAAGAAGTGGCTGATCTCAGCCACCGCGAAAATGACATCAAAAACGCCATTAACCTGATGTTCTGGGGAATATAAATGTCAGGCTCCCTTATACACAGgcgatgttgaagaccaCGCTGAGGTGTCAATCGTCGGAGCCGCTGAGCAATAACTAGCATAACCCCTTGGGGCCTCTAAACGGGTCTTGAGGGGTTTTTTGCATGGTCATAGCTGTTTCCTGAGAGCTTGGCAGGTGATGACACACATTAACAAATTTCGTGAGGAGTCTCCAGAAGAATGCCATTAATTTCCATAGGCTCCGCCCCCCTGACGAGCATCACAAAAATCGACGCTCAAGTCAGAGGTGGCGAAACCCGACAGGACTATAAAGATACCAGGCGTTTCCCCCTGGAAGCTCCCTCGTGCGCTCTCCTGTTCCGACCCTGCCGCTTACCGGATACCTGTCCGCCTTTCTCCCTTCGGGAAGCGTGGCGCTTTCTCATAGCTCACGCTGTAGGTATCTCAGTTCGGTGTAGGTCGTTCGCTCCAAGCTGGGCTGTGTGCACGAACCCCCCGTTCAGCCCGACCGCTGCGCCTTATCCGG", true}

func TestCutWithEnzymeByName(t *testing.T) {
	enzymeManager := NewEnzymeManager(GetBaseRestrictionEnzymes())
	_, err := enzymeManager.CutWithEnzymeByName(popen, true, "EcoFake")
	if err == nil {
		t.Errorf("CutWithEnzymeByName should have failed when looking for fake restriction enzyme EcoFake")
	}
}

func TestCutWithEnzyme(t *testing.T) {
	enzymeManager := NewEnzymeManager(GetBaseRestrictionEnzymes())
	var sequence Part
	bsai := "GGTCTCAATGC"
	bsaiComplement := "ATGCAGAGACC"

	// test(1)
	// Test case of `<-bsaiComplement bsai-> <-bsaiComplement bsai->` where bsaI cuts off of a linear sequence. This tests the line:
	// if !sequence.Circular && (overhangSet[len(overhangSet)-1].Position+enzyme.EnzymeSkip+enzyme.EnzymeOverhangLen > len(sequence))
	sequence = Part{"ATATATA" + bsaiComplement + bsai + "ATGCATCGATCGACTAGCATG" + bsaiComplement + bsai[:8], false}
	fragment, err := enzymeManager.CutWithEnzymeByName(sequence, true, "BsaI")
	if err != nil {
		t.Errorf("CutWithEnzyme should not have failed on test(1). Got error: %s", err)
	}
	if len(fragment) != 1 {
		t.Errorf("CutWithEnzyme in test(1) should be 1 fragment in length")
	}
	if fragment[0].Sequence != "ATGCATCGATCGACTAGCATG" {
		t.Errorf("CutWithEnzyme in test(1) should give fragment with sequence ATGCATCGATCGACTAGCATG . Got sequence: %s", fragment[0].Sequence)
	}

	// test(2)
	// Now if we take the same sequence and circularize it, we get a different result
	sequence.Circular = true
	fragment, err = enzymeManager.CutWithEnzymeByName(sequence, true, "BsaI")
	if err != nil {
		t.Errorf("CutWithEnzyme should not have failed on test(2). Got error: %s", err)
	}
	if len(fragment) != 2 {
		t.Errorf("CutWithEnzyme in test(2) should be 1 fragment in length")
	}
	if fragment[0].Sequence != "ATGCATCGATCGACTAGCATG" || fragment[1].Sequence != "TATA" {
		t.Errorf("CutWithEnzyme in test(2) should give fragment with sequence ATGCATCGATCGACTAGCATG and TATA. Got sequence: %s and %s", fragment[0].Sequence, fragment[1].Sequence)
	}

	// test(3)
	// Let's test if we have a single cut in our plasmid. This should give
	// different results if we have a linear or circular DNA. Since single cuts
	// will give no fragments if you test for directionality, we set the
	// directionality flag to false. This tests the line:
	// if len(overhangs) == 1 && !directional && !sequence.Circular
	sequence = Part{"ATATATATATATATAT" + bsai + "GCGCGCGCGCGCGCGCGCGC", false}
	fragment, err = enzymeManager.CutWithEnzymeByName(sequence, false, "BsaI")
	if err != nil {
		t.Errorf("CutWithEnzyme should not have failed on test(3). Got error: %s", err)
	}
	if len(fragment) != 2 {
		t.Errorf("Cutting a linear fragment with a single cut site should give 2 fragments")
	}
	if fragment[0].Sequence != "GCGCGCGCGCGCGCGCGCGC" || fragment[1].Sequence != "ATATATATATATATATGGTCTCA" {
		t.Errorf("CutWithEnzyme in test(3) should give fragment with sequence GCGCGCGCGCGCGCGCGCGC and ATATATATATATATATGGTCTCA. Got sequence: %s and %s", fragment[0].Sequence, fragment[1].Sequence)
	}

	// test(4)
	// This tests for the above except with a circular fragment. Specifically, it
	// tests the line:
	// if len(overhangs) == 2 && !directional && sequence.Circular
	sequence.Circular = true
	fragment, err = enzymeManager.CutWithEnzymeByName(sequence, false, "BsaI")
	if err != nil {
		t.Errorf("CutWithEnzyme should not have failed on test(4). Got error: %s", err)
	}
	if len(fragment) != 1 {
		t.Errorf("Cutting a circular fragment with a single cut site should give 1 fragments")
	}
	if fragment[0].Sequence != "GCGCGCGCGCGCGCGCGCGCATATATATATATATATGGTCTCA" {
		t.Errorf("CutWithEnzyme in test(4) should give fragment with sequence ATATATATATATATATGGTCTCA. Got Sequence: %s", fragment[0].Sequence)
	}

	// test(5)
	// This tests if we have a fragment where we do not care about directionality
	// but have more than 1 cut site in our fragment. We can use pOpen for this.
	fragment, err = enzymeManager.CutWithEnzymeByName(popen, false, "BbsI")
	if err != nil {
		t.Errorf("CutWithEnzyme should not have failed on test(5). Got error: %s", err)
	}
	if len(fragment) != 2 {
		t.Errorf("Cutting pOpen without a direction should yield 2 fragments")
	}
}

func TestCircularLigate(t *testing.T) {
	// The following tests for complementing overhangs. Specific, this line:
	// newSeed := Fragment{seedFragment.Sequence + seedFragment.ReverseOverhang + ReverseComplement(newFragment.Sequence), seedFragment.ForwardOverhang, ReverseComplement(newFragment.ForwardOverhang)}
	fragment1 := Fragment{"AAAAAA", "GTTG", "CTAT"}
	fragment2 := Fragment{"AAAAAA", "CAAC", "ATAG"}
	outputConstructs, infiniteLoops := CircularLigate([]Fragment{fragment1, fragment2})
	if len(outputConstructs) != 1 {
		t.Errorf("Circular ligation with complementing overhangs should only output 1 valid rotated sequence.")
	}
	if len(infiniteLoops) != 0 {
		t.Errorf("Circular ligation should have no loops")
	}
}

func TestGoldenGate(t *testing.T) {
	enzymeManager := NewEnzymeManager(GetBaseRestrictionEnzymes())
	// Here we test if the enzyme we want to use in a GoldenGate reaction does not exist in our enzyme pool
	fragment1 := Part{"GAAGTGCCATTCCGCCTGACCTGAAGACCAGGAGAAACACGTGGCAAACATTCCGGTCTCAAATGGAAAAGAGCAACGAAACCAACGGCTACCTTGACAGCGCTCAAGCCGGCCCTGCAGCTGGCCCGGGCGCTCCGGGTACCGCCGCGGGTCGTGCACGTCGTTGCGCGGGCTTCCTGCGGCGCCAAGCGCTGGTGCTGCTCACGGTGTCTGGTGTTCTGGCAGGCGCCGGTTTGGGCGCGGCACTGCGTGGGCTCAGCCTGAGCCGCACCCAGGTCACCTACCTGGCCTTCCCCGGCGAGATGCTGCTCCGCATGCTGCGCATGATCATCCTGCCGCTGGTGGTCTGCAGCCTGGTGTCGGGCGCCGCCTCCCTCGATGCCAGCTGCCTCGGGCGTCTGGGCGGTATCGCTGTCGCCTACTTTGGCCTCACCACACTGAGTGCCTCGGCGCTCGCCGTGGCCTTGGCGTTCATCATCAAGCCAGGATCCGGTGCGCAGACCCTTCAGTCCAGCGACCTGGGGCTGGAGGACTCGGGGCCTCCTCCTGTCCCCAAAGAAACGGTGGACTCTTTCCTCGACCTGGCCAGAAACCTGTTTCCCTCCAATCTTGTGGTTGCAGCTTTCCGTACGTATGCAACCGATTATAAAGTCGTGACCCAGAACAGCAGCTCTGGAAATGTAACCCATGAAAAGATCCCCATAGGCACTGAGATAGAAGGGATGAACATTTTAGGATTGGTCCTGTTTGCTCTGGTGTTAGGAGTGGCCTTAAAGAAACTAGGCTCCGAAGGAGAGGACCTCATCCGTTTCTTCAATTCCCTCAACGAGGCGACGATGGTGCTGGTGTCCTGGATTATGTGGTACGCGTCTTCAGGCTAGGTGGAGGCTCAGTG", false}
	fragment2 := Part{"GAAGTGCCATTCCGCCTGACCTGAAGACCAGTACGTACCTGTGGGCATCATGTTCCTTGTTGGAAGCAAGATCGTGGAAATGAAAGACATCATCGTGCTGGTGACCAGCCTGGGGAAATACATCTTCGCATCTATATTGGGCCACGTCATTCATGGTGGTATCGTCCTGCCGCTGATTTATTTTGTTTTCACACGAAAAAACCCATTCAGATTCCTCCTGGGCCTCCTCGCCCCATTTGCGACAGCATTTGCTACGTGCTCCAGCTCAGCGACCCTTCCCTCTATGATGAAGTGCATTGAAGAGAACAATGGTGTGGACAAGAGGATCTCCAGGTTTATTCTCCCCATCGGGGCCACCGTGAACATGGACGGAGCAGCCATCTTCCAGTGTGTGGCCGCGGTGTTCATTGCGCAACTCAACAACGTAGAGCTCAACGCAGGACAGATTTTCACCATTCTAGTGACTGCCACAGCGTCCAGTGTTGGAGCAGCAGGCGTGCCAGCTGGAGGGGTCCTCACCATTGCCATTATCCTGGAGGCCATTGGGCTGCCTACTCATGATCTGCCTCTGATCCTGGCTGTGGACTGGATTGTGGACCGGACCACCACGGTGGTGAATGTGGAAGGGGATGCCCTGGGTGCAGGCATTCTCCACCACCTGAATCAGAAGGCAACAAAGAAAGGCGAGCAGGAACTTGCTGAGGTGAAAGTGGAAGCCATCCCCAACTGCAAGTCTGAGGAGGAAACCTCGCCCCTGGTGACACACCAGAACCCCGCTGGCCCCGTGGCCAGTGCCCCAGAACTGGAATCCAAGGAGTCGGTTCTGTGAAGAGCTTAGAGACCGACGACTGCCTAAGGACATTCGCTGCGTCTTCAGGCTAGGTGGAGGCTCAGTG", false}

	_, _, err := enzymeManager.GoldenGate([]Part{fragment1, fragment2, popen}, "EcoRFake")
	if err == nil {
		t.Errorf("GoldenGate should fail when using enzyme EcoRFake")
	}
	if err.Error() != "Enzyme EcoRFake not found" {
		t.Errorf("Failure of GoldenGate on incorrect enzyme should follow the exact string `Enzyme EcoRFake not found in enzymeMap`. Got: %s", err.Error())
	}
}

func TestSignalKilledGoldenGate(t *testing.T) {
	enzymeManager := NewEnzymeManager(GetBaseRestrictionEnzymes())
	// This previously would crash from using too much RAM.
	fragment1 := Part{"AAAGCACTCTTAGGCCTCTGGAAGACATGGAGGGTCTCAAGGTGATCAAAGGATCTTCTTGAGATCCTTTTTTTCTGCGCGTAATCTTTTGCCCTGTAAACGAAAAAACCACCTGGGTAGTCTTCGCATTTCTTAATCGGTGCCC", false}
	fragment2 := Part{"AAAGCACTCTTAGGCCTCTGGAAGACATTGGGGAGGTGGTTTGATCGAAGGTTAAGTCAGTTGGGGAACTGCTTAACCGTGGTAACTGGCTTTCGCAGAGCACAGCAACCAAATCTGTTAGTCTTCGCATTTCTTAATCGGTGCCC", false}
	fragment3 := Part{"AAAGCACTCTTAGGCCTCTGGAAGACATCTGTCCTTCCAGTGTAGCCGGACTTTGGCGCACACTTCAAGAGCAACCGCGTGTTTAGCTAAACAAATCCTCTGCGAACTCCCAGTTACCTAGTCTTCGCATTTCTTAATCGGTGCCC", false}
	fragment4 := Part{"AAAGCACTCTTAGGCCTCTGGAAGACATTACCAATGGCTGCTGCCAGTGGCGTTTTACCGTGCTTTTCCGGGTTGGACTCAAGTGAACAGTTACCGGATAAGGCGCAGCAGTCGGGCTTAGTCTTCGCATTTCTTAATCGGTGCCC", false}
	fragment5 := Part{"AAAGCACTCTTAGGCCTCTGGAAGACATGGCTGAACGGGGAGTTCTTGCTTACAGCCCAGCTTGGAGCGAACGACCTACACCGAGCCGAGATACCAGTGTGTGAGCTATGAGAAAGCGTAGTCTTCGCATTTCTTAATCGGTGCCC", false}
	fragment6 := Part{"AAAGCACTCTTAGGCCTCTGGAAGACATAGCGCCACACTTCCCGTAAGGGAGAAAGGCGGAACAGGTATCCGGTAAACGGCAGGGTCGGAACAGGAGAGCGCAAGAGGGAGCGACCCGTAGTCTTCGCATTTCTTAATCGGTGCCC", false}
	fragment7 := Part{"AAAGCACTCTTAGGCCTCTGGAAGACATCCCGCCGGAAACGGTGGGGATCTTTAAGTCCTGTCGGGTTTCGCCCGTACTGTCAGATTCATGGTTGAGCCTCACGGCTCCCACAGATGTAGTCTTCGCATTTCTTAATCGGTGCCC", false}
	fragment8 := Part{"AAAGCACTCTTAGGCCTCTGGAAGACATGATGCACCGGAAAAGCGTCTGTTTATGTGAACTCTGGCAGGAGGGCGGAGCCTATGGAAAAACGCCACCGGCGCGGCCCTGCTGTTTTGCCTCACATGTTAGTCTTCGCATTTCTTAATCGGTGCCC", false}
	fragment9 := Part{"AAAGCACTCTTAGGCCTCTGGAAGACATATGTTAGTCCCCTGCTTATCCACGGAATCTGTGGGTAACTTTGTATGTGTCCGCAGCGCAAAAAGAGACCCGCTTAGTCTTCGCATTTCTTAATCGGTGCCC", false}
	fragments := []Part{popen, fragment1, fragment2, fragment3, fragment4, fragment5, fragment6, fragment7, fragment8, fragment9}

	clones, loopingClones, err := enzymeManager.GoldenGate(fragments, "BbsI")
	if err != nil {
		t.Errorf("GoldenGate should not fail with these fragments. Got error: %s", err)
	}
	if len(clones) != 1 {
		t.Errorf("There should be 1 output  Got: %d", len(clones))
	}
	// This should be changed later when we have a better way of informing user of reused overhangs
	if len(loopingClones) != 4 {
		t.Errorf("Should be only 4 looping sequences. Got: %d", len(loopingClones))
	}
}

func TestPanicGoldenGate(t *testing.T) {
	enzymeManager := NewEnzymeManager(GetBaseRestrictionEnzymes())
	// This used to panic with the message:
	// panic: runtime error: slice bounds out of range [:-2] [recovered]
	// It was from the following sequence: GAAGACATAATGGTCTTC . There are 2 intercepting BbsI sites.
	fragment1 := Part{"AAACCGGAGCCATACAGTACGAAGACATGGAGGGTCTCAAATGAAAAAAATCATCGAAACCCAGCGTGCACCGGGAGCAATCGGACCGTACGTCCAGGGAGTCGACCTAGGATCAATGTAGTCTTCGCACTTGGCTTAGATGCAAC", false}
	fragment2 := Part{"AAACCGGAGCCATACAGTACGAAGACATAATGGTCTTCACCTCAGGACAGATCCCGGTCTGCCCGCAGACCGGAGAAATCCCGGCAGACGTCCAGGACCAGGCACGTCTATCACTAGATAGTCTTCGCACTTGGCTTAGATGCAAC", false}
	fragment3 := Part{"AAACCGGAGCCATACAGTACGAAGACATTAGAAAACGTCAAAGCAATCGTCGTCGCAGCAGGACTATCAGTCGGAGACATCATCAAAATGACCGTCTTCATCACCGACCTAAACGACTTAGTCTTCGCACTTGGCTTAGATGCAAC", false}
	fragment4 := Part{"AAACCGGAGCCATACAGTACGAAGACATGACTTCGCAACCATCAACGAAGTCTACAAACAGTTCTTCGACGAACACCAGGCAACCTACCCGACCCGTTCATGCGTCCAGGTCGCACGTCTACTAGTCTTCGCACTTGGCTTAGATGCAAC", false}
	fragment5 := Part{"AAACCGGAGCCATACAGTACGAAGACATCTACCGAAAGACGTCAAACTAGAAATCGAAGCAATCGCAGTCCGTTCAGCAAGAGCTTAGAGACCCGCTTAGTCTTCGCACTTGGCTTAGATGCAAC", false}
	fragments := []Part{popen, fragment1, fragment2, fragment3, fragment4, fragment5}

	_, _, err := enzymeManager.GoldenGate(fragments, "BbsI")
	if err != nil {
		t.Errorf("GoldenGate should not fail with these fragments. Got error: %s", err)
	}
}

func TestCircularCutRegression(t *testing.T) {
	enzymeManager := NewEnzymeManager(GetBaseRestrictionEnzymes())
	// This used to error with 0 fragments since the BsaI cut site is on the other
	// side of the origin from its recognition site.
	plasmid1 := Part{"AAACTACAAGACCCGCGCCGAGGTGAAGTTCGAGGGCGACACCCTGGTGAACCGCATCGAGCTGAAGGGCATCGACTTCAAGGAGGACGGCAACATCCTGGGGCACAAGCTGGAGTACAACTACAACAGCCACAACGTCTATATCATGGCCGACAAGCAGAAGAACGGCATCAAGGTGAACTTCAAGATCCGCCACAACATCGAGGACGGCAGCCGAGaccaagtcgcggccgcgaggtgtcaatcgtcggagtagggataacagggtaatccgctgagcaataactagcataaccccttggggcctctaaacgggtcttgaggggttttttgcatggtcatagctgtttcctgttacgccccgccctgccactcgtcgcagtactgttgtaattcattaagcattctgccgacatggaagccatcacaaacggcatgatgaacctgaatcgccagcggcatcagcaccttgtcgccttgcgtataatatttgcccatggtgaaaacgggggcgaagaagttgtccatattggccacgtttaaatcaaaactggtgaaactcacccagggattggctgacacgaaaaacatattctcaataaaccctttagggaaataggccaggttttcaccgtaacacgccacatcttgcgaatatatgtgtagaaactgccggaaatcgtcgtggtattcactccagagggatgaaaacgtttcagtttgctcatggaaaacggtgtaacaagggtgaacactatcccatatcaccagctcaccatccttcattgccatacgaaattccggatgagcattcatcaggcgggcaagaatgtgaataaaggccggataaaacttgtgcttatttttctttacggtctttaaaaaggccgtaatatccagctgaacggtctggttataggtacattgagcaactgactgaaatgcctcaaaatgttctttacgatgccattgggatatatcaacggtggtatatccagtgatttttttctccattttagcttccttagctcctgaaaatctcgataactcaaaaaatacgcccggtagtgatcttatttcattatggtgaaagttggaacctcttacgtgccgatcatttccataggctccgcccccctgacgagcatcacaaaaatcgacgctcaagtcagaggtggcgaaacccgacaggactataaagataccaggcgtttccccctggaagctccctcgtgcgctctcctgttccgaccctgccgcttaccggatacctgtccgcctttctcccttcgggaagcgtggcgctttctcatagctcacgctgtaggtatctcagttcggtgtaggtcgttcgctccaagctgggctgtgtgcacgaaccccccgttcagcccgaccgctgcgccttatccggtaactatcgtcttgagtccaacccggtaagacacgacttatcgccactggcagcagccactggtaacaggattagcagagcgaggtatgtaggcggtgctacagagttcttgaagtggtggcctaactacggctacactagaaggacagtatttggtatctgcgctctgctgaagccagttaccttcggaaaaagagttggtagctcttgatccggcaaacaaaccaccgctggtagcggtggtttttttgtttgcaagcagcagattacgcgcagaaaaaaaggatctcaagtaaaacgacggccagtagtcaaaagcctccgaccggaggcttttgacttggttcaggtggagtggcggccgcgacttgGTCTC", true}
	newFragments, err := enzymeManager.CutWithEnzymeByName(plasmid1, true, "BsaI")
	if err != nil {
		t.Errorf("Failed to cut: %s", err)
	}
	if len(newFragments) != 1 {
		t.Errorf("Expected 1 new fragment, got: %d", len(newFragments))
	}
}

func benchmarkGoldenGate(b *testing.B, enzymeManager EnzymeManager, parts []Part) {
	for n := 0; n < b.N; n++ {
		_, _, _ = enzymeManager.GoldenGate(parts, "BbsI")
	}
}

func BenchmarkGoldenGate3Parts(b *testing.B) {
	enzymeManager := NewEnzymeManager(GetBaseRestrictionEnzymes())
	fragment1 := Part{"GAAGTGCCATTCCGCCTGACCTGAAGACCAGGAGAAACACGTGGCAAACATTCCGGTCTCAAATGGAAAAGAGCAACGAAACCAACGGCTACCTTGACAGCGCTCAAGCCGGCCCTGCAGCTGGCCCGGGCGCTCCGGGTACCGCCGCGGGTCGTGCACGTCGTTGCGCGGGCTTCCTGCGGCGCCAAGCGCTGGTGCTGCTCACGGTGTCTGGTGTTCTGGCAGGCGCCGGTTTGGGCGCGGCACTGCGTGGGCTCAGCCTGAGCCGCACCCAGGTCACCTACCTGGCCTTCCCCGGCGAGATGCTGCTCCGCATGCTGCGCATGATCATCCTGCCGCTGGTGGTCTGCAGCCTGGTGTCGGGCGCCGCCTCCCTCGATGCCAGCTGCCTCGGGCGTCTGGGCGGTATCGCTGTCGCCTACTTTGGCCTCACCACACTGAGTGCCTCGGCGCTCGCCGTGGCCTTGGCGTTCATCATCAAGCCAGGATCCGGTGCGCAGACCCTTCAGTCCAGCGACCTGGGGCTGGAGGACTCGGGGCCTCCTCCTGTCCCCAAAGAAACGGTGGACTCTTTCCTCGACCTGGCCAGAAACCTGTTTCCCTCCAATCTTGTGGTTGCAGCTTTCCGTACGTATGCAACCGATTATAAAGTCGTGACCCAGAACAGCAGCTCTGGAAATGTAACCCATGAAAAGATCCCCATAGGCACTGAGATAGAAGGGATGAACATTTTAGGATTGGTCCTGTTTGCTCTGGTGTTAGGAGTGGCCTTAAAGAAACTAGGCTCCGAAGGAGAGGACCTCATCCGTTTCTTCAATTCCCTCAACGAGGCGACGATGGTGCTGGTGTCCTGGATTATGTGGTACGCGTCTTCAGGCTAGGTGGAGGCTCAGTG", false}
	fragment2 := Part{"GAAGTGCCATTCCGCCTGACCTGAAGACCAGTACGTACCTGTGGGCATCATGTTCCTTGTTGGAAGCAAGATCGTGGAAATGAAAGACATCATCGTGCTGGTGACCAGCCTGGGGAAATACATCTTCGCATCTATATTGGGCCACGTCATTCATGGTGGTATCGTCCTGCCGCTGATTTATTTTGTTTTCACACGAAAAAACCCATTCAGATTCCTCCTGGGCCTCCTCGCCCCATTTGCGACAGCATTTGCTACGTGCTCCAGCTCAGCGACCCTTCCCTCTATGATGAAGTGCATTGAAGAGAACAATGGTGTGGACAAGAGGATCTCCAGGTTTATTCTCCCCATCGGGGCCACCGTGAACATGGACGGAGCAGCCATCTTCCAGTGTGTGGCCGCGGTGTTCATTGCGCAACTCAACAACGTAGAGCTCAACGCAGGACAGATTTTCACCATTCTAGTGACTGCCACAGCGTCCAGTGTTGGAGCAGCAGGCGTGCCAGCTGGAGGGGTCCTCACCATTGCCATTATCCTGGAGGCCATTGGGCTGCCTACTCATGATCTGCCTCTGATCCTGGCTGTGGACTGGATTGTGGACCGGACCACCACGGTGGTGAATGTGGAAGGGGATGCCCTGGGTGCAGGCATTCTCCACCACCTGAATCAGAAGGCAACAAAGAAAGGCGAGCAGGAACTTGCTGAGGTGAAAGTGGAAGCCATCCCCAACTGCAAGTCTGAGGAGGAAACCTCGCCCCTGGTGACACACCAGAACCCCGCTGGCCCCGTGGCCAGTGCCCCAGAACTGGAATCCAAGGAGTCGGTTCTGTGAAGAGCTTAGAGACCGACGACTGCCTAAGGACATTCGCTGCGTCTTCAGGCTAGGTGGAGGCTCAGTG", false}

	benchmarkGoldenGate(b, enzymeManager, []Part{fragment1, fragment2, popen})
}
