package clone

import (
	"fmt"
	"regexp"
	"testing"

	"github.com/TimothyStiles/poly/seqhash"
)

// Eventually, we want to get the data for this map from ftp://ftp.neb.com/pub/rebase
func getBaseRestrictionEnzymes() []Enzyme {
	return []Enzyme{
		{"BsaI", regexp.MustCompile("GGTCTC"), regexp.MustCompile("GAGACC"), 1, 4, "GGTCTC"},
		{"BbsI", regexp.MustCompile("GAAGAC"), regexp.MustCompile("GTCTTC"), 2, 4, "GAAGAC"},
		{"BtgZI", regexp.MustCompile("GCGATG"), regexp.MustCompile("CATCGC"), 10, 4, "GCGATG"},
	}
}

// pOpen plasmid series (https://stanford.freegenes.org/collections/open-genes/products/open-plasmids#description). I use it for essentially all my cloning. -Keoni
var popen = Part{"TAACTATCGTCTTGAGTCCAACCCGGTAAGACACGACTTATCGCCACTGGCAGCAGCCACTGGTAACAGGATTAGCAGAGCGAGGTATGTAGGCGGTGCTACAGAGTTCTTGAAGTGGTGGCCTAACTACGGCTACACTAGAAGAACAGTATTTGGTATCTGCGCTCTGCTGAAGCCAGTTACCTTCGGAAAAAGAGTTGGTAGCTCTTGATCCGGCAAACAAACCACCGCTGGTAGCGGTGGTTTTTTTGTTTGCAAGCAGCAGATTACGCGCAGAAAAAAAGGATCTCAAGAAGGCCTACTATTAGCAACAACGATCCTTTGATCTTTTCTACGGGGTCTGACGCTCAGTGGAACGAAAACTCACGTTAAGGGATTTTGGTCATGAGATTATCAAAAAGGATCTTCACCTAGATCCTTTTAAATTAAAAATGAAGTTTTAAATCAATCTAAAGTATATATGAGTAAACTTGGTCTGACAGTTACCAATGCTTAATCAGTGAGGCACCTATCTCAGCGATCTGTCTATTTCGTTCATCCATAGTTGCCTGACTCCCCGTCGTGTAGATAACTACGATACGGGAGGGCTTACCATCTGGCCCCAGTGCTGCAATGATACCGCGAGAACCACGCTCACCGGCTCCAGATTTATCAGCAATAAACCAGCCAGCCGGAAGGGCCGAGCGCAGAAGTGGTCCTGCAACTTTATCCGCCTCCATCCAGTCTATTAATTGTTGCCGGGAAGCTAGAGTAAGTAGTTCGCCAGTTAATAGTTTGCGCAACGTTGTTGCCATTGCTACAGGCATCGTGGTGTCACGCTCGTCGTTTGGTATGGCTTCATTCAGCTCCGGTTCCCAACGATCAAGGCGAGTTACATGATCCCCCATGTTGTGCAAAAAAGCGGTTAGCTCCTTCGGTCCTCCGATCGTTGTCAGAAGTAAGTTGGCCGCAGTGTTATCACTCATGGTTATGGCAGCACTGCATAATTCTCTTACTGTCATGCCATCCGTAAGATGCTTTTCTGTGACTGGTGAGTACTCAACCAAGTCATTCTGAGAATAGTGTATGCGGCGACCGAGTTGCTCTTGCCCGGCGTCAATACGGGATAATACCGCGCCACATAGCAGAACTTTAAAAGTGCTCATCATTGGAAAACGTTCTTCGGGGCGAAAACTCTCAAGGATCTTACCGCTGTTGAGATCCAGTTCGATGTAACCCACTCGTGCACCCAACTGATCTTCAGCATCTTTTACTTTCACCAGCGTTTCTGGGTGAGCAAAAACAGGAAGGCAAAATGCCGCAAAAAAGGGAATAAGGGCGACACGGAAATGTTGAATACTCATACTCTTCCTTTTTCAATATTATTGAAGCATTTATCAGGGTTATTGTCTCATGAGCGGATACATATTTGAATGTATTTAGAAAAATAAACAAATAGGGGTTCCGCGCACCTGCACCAGTCAGTAAAACGACGGCCAGTAGTCAAAAGCCTCCGACCGGAGGCTTTTGACTTGGTTCAGGTGGAGTGGGAGTAgtcttcGCcatcgCtACTAAAagccagataacagtatgcgtatttgcgcgctgatttttgcggtataagaatatatactgatatgtatacccgaagtatgtcaaaaagaggtatgctatgaagcagcgtattacagtgacagttgacagcgacagctatcagttgctcaaggcatatatgatgtcaatatctccggtctggtaagcacaaccatgcagaatgaagcccgtcgtctgcgtgccgaacgctggaaagcggaaaatcaggaagggatggctgaggtcgcccggtttattgaaatgaacggctcttttgctgacgagaacagggGCTGGTGAAATGCAGTTTAAGGTTTACACCTATAAAAGAGAGAGCCGTTATCGTCTGTTTGTGGATGTACAGAGTGATATTATTGACACGCCCGGGCGACGGATGGTGATCCCCCTGGCCAGTGCACGTCTGCTGTCAGATAAAGTCTCCCGTGAACTTTACCCGGTGGTGCATATCGGGGATGAAAGCTGGCGCATGATGACCACCGATATGGCCAGTGTGCCGGTCTCCGTTATCGGGGAAGAAGTGGCTGATCTCAGCCACCGCGAAAATGACATCAAAAACGCCATTAACCTGATGTTCTGGGGAATATAAATGTCAGGCTCCCTTATACACAGgcgatgttgaagaccaCGCTGAGGTGTCAATCGTCGGAGCCGCTGAGCAATAACTAGCATAACCCCTTGGGGCCTCTAAACGGGTCTTGAGGGGTTTTTTGCATGGTCATAGCTGTTTCCTGAGAGCTTGGCAGGTGATGACACACATTAACAAATTTCGTGAGGAGTCTCCAGAAGAATGCCATTAATTTCCATAGGCTCCGCCCCCCTGACGAGCATCACAAAAATCGACGCTCAAGTCAGAGGTGGCGAAACCCGACAGGACTATAAAGATACCAGGCGTTTCCCCCTGGAAGCTCCCTCGTGCGCTCTCCTGTTCCGACCCTGCCGCTTACCGGATACCTGTCCGCCTTTCTCCCTTCGGGAAGCGTGGCGCTTTCTCATAGCTCACGCTGTAGGTATCTCAGTTCGGTGTAGGTCGTTCGCTCCAAGCTGGGCTGTGTGCACGAACCCCCCGTTCAGCCCGACCGCTGCGCCTTATCCGG", true}

func TestCutWithEnzymeByName(t *testing.T) {
	em := NewEnzymeManager(getBaseRestrictionEnzymes())
	_, err := em.CutWithEnzymeByName(popen, true, "EcoFake")
	if err == nil {
		t.Errorf("CutWithEnzymeByName should have failed when looking for fake restriction enzyme EcoFake")
	}
}

func TestCutWithEnzyme(t *testing.T) {
	em := NewEnzymeManager(getBaseRestrictionEnzymes())
	var seq Part
	bsai := "GGTCTCAATGC"
	bsaiComplement := "ATGCAGAGACC"

	// test(1)
	// Test case of `<-bsaiComplement bsai-> <-bsaiComplement bsai->` where bsaI cuts off of a linear sequence. This tests the line:
	// if !seq.Circular && (overhangSet[len(overhangSet)-1].Position+enzyme.EnzymeSkip+enzyme.EnzymeOverhangLen > len(sequence))
	seq = Part{"ATATATA" + bsaiComplement + bsai + "ATGCATCGATCGACTAGCATG" + bsaiComplement + bsai[:8], false}
	frag, err := em.CutWithEnzymeByName(seq, true, "BsaI")
	if err != nil {
		t.Errorf("CutWithEnzyme should not have failed on test(1). Got error: %s", err)
	}
	if len(frag) != 1 {
		t.Errorf("CutWithEnzyme in test(1) should be 1 fragment in length")
	}
	if frag[0].Sequence != "ATGCATCGATCGACTAGCATG" {
		t.Errorf("CutWithEnzyme in test(1) should give fragment with sequence ATGCATCGATCGACTAGCATG . Got sequence: %s", frag[0].Sequence)
	}

	// test(2)
	// Now if we take the same sequence and circularize it, we get a different result
	seq.Circular = true
	frag, err = em.CutWithEnzymeByName(seq, true, "BsaI")
	if err != nil {
		t.Errorf("CutWithEnzyme should not have failed on test(2). Got error: %s", err)
	}
	if len(frag) != 2 {
		t.Errorf("CutWithEnzyme in test(2) should be 1 fragment in length")
	}
	if frag[0].Sequence != "ATGCATCGATCGACTAGCATG" || frag[1].Sequence != "TATA" {
		t.Errorf("CutWithEnzyme in test(2) should give fragment with sequence ATGCATCGATCGACTAGCATG and TATA. Got sequence: %s and %s", frag[0].Sequence, frag[1].Sequence)
	}

	// test(3)
	// Let's test if we have a single cut in our plasmid. This should give
	// different results if we have a linear or circular DNA. Since single cuts
	// will give no fragments if you test for directionality, we set the
	// directionality flag to false. This tests the line:
	// if len(overhangs) == 1 && !directional && !seq.Circular
	seq = Part{"ATATATATATATATAT" + bsai + "GCGCGCGCGCGCGCGCGCGC", false}
	frag, err = em.CutWithEnzymeByName(seq, false, "BsaI")
	if err != nil {
		t.Errorf("CutWithEnzyme should not have failed on test(3). Got error: %s", err)
	}
	if len(frag) != 2 {
		t.Errorf("Cutting a linear fragment with a single cut site should give 2 fragments")
	}
	if frag[0].Sequence != "GCGCGCGCGCGCGCGCGCGC" || frag[1].Sequence != "ATATATATATATATATGGTCTCA" {
		t.Errorf("CutWithEnzyme in test(3) should give fragment with sequence GCGCGCGCGCGCGCGCGCGC and ATATATATATATATATGGTCTCA. Got sequence: %s and %s", frag[0].Sequence, frag[1].Sequence)
	}

	// test(4)
	// This tests for the above except with a circular fragment. Specifically, it
	// tests the line:
	// if len(overhangs) == 2 && !directional && seq.Circular
	seq.Circular = true
	frag, err = em.CutWithEnzymeByName(seq, false, "BsaI")
	if err != nil {
		t.Errorf("CutWithEnzyme should not have failed on test(4). Got error: %s", err)
	}
	if len(frag) != 1 {
		t.Errorf("Cutting a circular fragment with a single cut site should give 1 fragments")
	}
	if frag[0].Sequence != "GCGCGCGCGCGCGCGCGCGCATATATATATATATATGGTCTCA" {
		t.Errorf("CutWithEnzyme in test(4) should give fragment with sequence ATATATATATATATATGGTCTCA. Got Sequence: %s", frag[0].Sequence)
	}

	// test(5)
	// This tests if we have a fragment where we do not care about directionality
	// but have more than 1 cut site in our fragment. We can use pOpen for this.
	frag, err = em.CutWithEnzymeByName(popen, false, "BbsI")
	if err != nil {
		t.Errorf("CutWithEnzyme should not have failed on test(5). Got error: %s", err)
	}
	if len(frag) != 2 {
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
	em := NewEnzymeManager(getBaseRestrictionEnzymes())
	// Here we test if the enzyme we want to use in a GoldenGate reaction does not exist in our enzyme pool
	fragment1 := Part{"GAAGTGCCATTCCGCCTGACCTGAAGACCAGGAGAAACACGTGGCAAACATTCCGGTCTCAAATGGAAAAGAGCAACGAAACCAACGGCTACCTTGACAGCGCTCAAGCCGGCCCTGCAGCTGGCCCGGGCGCTCCGGGTACCGCCGCGGGTCGTGCACGTCGTTGCGCGGGCTTCCTGCGGCGCCAAGCGCTGGTGCTGCTCACGGTGTCTGGTGTTCTGGCAGGCGCCGGTTTGGGCGCGGCACTGCGTGGGCTCAGCCTGAGCCGCACCCAGGTCACCTACCTGGCCTTCCCCGGCGAGATGCTGCTCCGCATGCTGCGCATGATCATCCTGCCGCTGGTGGTCTGCAGCCTGGTGTCGGGCGCCGCCTCCCTCGATGCCAGCTGCCTCGGGCGTCTGGGCGGTATCGCTGTCGCCTACTTTGGCCTCACCACACTGAGTGCCTCGGCGCTCGCCGTGGCCTTGGCGTTCATCATCAAGCCAGGATCCGGTGCGCAGACCCTTCAGTCCAGCGACCTGGGGCTGGAGGACTCGGGGCCTCCTCCTGTCCCCAAAGAAACGGTGGACTCTTTCCTCGACCTGGCCAGAAACCTGTTTCCCTCCAATCTTGTGGTTGCAGCTTTCCGTACGTATGCAACCGATTATAAAGTCGTGACCCAGAACAGCAGCTCTGGAAATGTAACCCATGAAAAGATCCCCATAGGCACTGAGATAGAAGGGATGAACATTTTAGGATTGGTCCTGTTTGCTCTGGTGTTAGGAGTGGCCTTAAAGAAACTAGGCTCCGAAGGAGAGGACCTCATCCGTTTCTTCAATTCCCTCAACGAGGCGACGATGGTGCTGGTGTCCTGGATTATGTGGTACGCGTCTTCAGGCTAGGTGGAGGCTCAGTG", false}
	fragment2 := Part{"GAAGTGCCATTCCGCCTGACCTGAAGACCAGTACGTACCTGTGGGCATCATGTTCCTTGTTGGAAGCAAGATCGTGGAAATGAAAGACATCATCGTGCTGGTGACCAGCCTGGGGAAATACATCTTCGCATCTATATTGGGCCACGTCATTCATGGTGGTATCGTCCTGCCGCTGATTTATTTTGTTTTCACACGAAAAAACCCATTCAGATTCCTCCTGGGCCTCCTCGCCCCATTTGCGACAGCATTTGCTACGTGCTCCAGCTCAGCGACCCTTCCCTCTATGATGAAGTGCATTGAAGAGAACAATGGTGTGGACAAGAGGATCTCCAGGTTTATTCTCCCCATCGGGGCCACCGTGAACATGGACGGAGCAGCCATCTTCCAGTGTGTGGCCGCGGTGTTCATTGCGCAACTCAACAACGTAGAGCTCAACGCAGGACAGATTTTCACCATTCTAGTGACTGCCACAGCGTCCAGTGTTGGAGCAGCAGGCGTGCCAGCTGGAGGGGTCCTCACCATTGCCATTATCCTGGAGGCCATTGGGCTGCCTACTCATGATCTGCCTCTGATCCTGGCTGTGGACTGGATTGTGGACCGGACCACCACGGTGGTGAATGTGGAAGGGGATGCCCTGGGTGCAGGCATTCTCCACCACCTGAATCAGAAGGCAACAAAGAAAGGCGAGCAGGAACTTGCTGAGGTGAAAGTGGAAGCCATCCCCAACTGCAAGTCTGAGGAGGAAACCTCGCCCCTGGTGACACACCAGAACCCCGCTGGCCCCGTGGCCAGTGCCCCAGAACTGGAATCCAAGGAGTCGGTTCTGTGAAGAGCTTAGAGACCGACGACTGCCTAAGGACATTCGCTGCGTCTTCAGGCTAGGTGGAGGCTCAGTG", false}

	_, _, err := em.GoldenGate([]Part{fragment1, fragment2, popen}, "EcoRFake")
	if err == nil {
		t.Errorf("GoldenGate should fail when using enzyme EcoRFake")
	}
	if err.Error() != "Enzyme EcoRFake not found" {
		t.Errorf("Failure of GoldenGate on incorrect enzyme should follow the exact string `Enzyme EcoRFake not found in enzymeMap`. Got: %s", err.Error())
	}
}

func ExampleGoldenGate() {
	em := NewEnzymeManager(getBaseRestrictionEnzymes())
	// Fragment 1 has a palindrome at its start. This isn't very common but
	// can occur. These two fragments are real DNA fragments used in the
	// FreeGenes Project. They are used because they were on my computer
	// - Keoni
	fragment1 := Part{"GAAGTGCCATTCCGCCTGACCTGAAGACCAGGAGAAACACGTGGCAAACATTCCGGTCTCAAATGGAAAAGAGCAACGAAACCAACGGCTACCTTGACAGCGCTCAAGCCGGCCCTGCAGCTGGCCCGGGCGCTCCGGGTACCGCCGCGGGTCGTGCACGTCGTTGCGCGGGCTTCCTGCGGCGCCAAGCGCTGGTGCTGCTCACGGTGTCTGGTGTTCTGGCAGGCGCCGGTTTGGGCGCGGCACTGCGTGGGCTCAGCCTGAGCCGCACCCAGGTCACCTACCTGGCCTTCCCCGGCGAGATGCTGCTCCGCATGCTGCGCATGATCATCCTGCCGCTGGTGGTCTGCAGCCTGGTGTCGGGCGCCGCCTCCCTCGATGCCAGCTGCCTCGGGCGTCTGGGCGGTATCGCTGTCGCCTACTTTGGCCTCACCACACTGAGTGCCTCGGCGCTCGCCGTGGCCTTGGCGTTCATCATCAAGCCAGGATCCGGTGCGCAGACCCTTCAGTCCAGCGACCTGGGGCTGGAGGACTCGGGGCCTCCTCCTGTCCCCAAAGAAACGGTGGACTCTTTCCTCGACCTGGCCAGAAACCTGTTTCCCTCCAATCTTGTGGTTGCAGCTTTCCGTACGTATGCAACCGATTATAAAGTCGTGACCCAGAACAGCAGCTCTGGAAATGTAACCCATGAAAAGATCCCCATAGGCACTGAGATAGAAGGGATGAACATTTTAGGATTGGTCCTGTTTGCTCTGGTGTTAGGAGTGGCCTTAAAGAAACTAGGCTCCGAAGGAGAGGACCTCATCCGTTTCTTCAATTCCCTCAACGAGGCGACGATGGTGCTGGTGTCCTGGATTATGTGGTACGCGTCTTCAGGCTAGGTGGAGGCTCAGTG", false}
	fragment2 := Part{"GAAGTGCCATTCCGCCTGACCTGAAGACCAGTACGTACCTGTGGGCATCATGTTCCTTGTTGGAAGCAAGATCGTGGAAATGAAAGACATCATCGTGCTGGTGACCAGCCTGGGGAAATACATCTTCGCATCTATATTGGGCCACGTCATTCATGGTGGTATCGTCCTGCCGCTGATTTATTTTGTTTTCACACGAAAAAACCCATTCAGATTCCTCCTGGGCCTCCTCGCCCCATTTGCGACAGCATTTGCTACGTGCTCCAGCTCAGCGACCCTTCCCTCTATGATGAAGTGCATTGAAGAGAACAATGGTGTGGACAAGAGGATCTCCAGGTTTATTCTCCCCATCGGGGCCACCGTGAACATGGACGGAGCAGCCATCTTCCAGTGTGTGGCCGCGGTGTTCATTGCGCAACTCAACAACGTAGAGCTCAACGCAGGACAGATTTTCACCATTCTAGTGACTGCCACAGCGTCCAGTGTTGGAGCAGCAGGCGTGCCAGCTGGAGGGGTCCTCACCATTGCCATTATCCTGGAGGCCATTGGGCTGCCTACTCATGATCTGCCTCTGATCCTGGCTGTGGACTGGATTGTGGACCGGACCACCACGGTGGTGAATGTGGAAGGGGATGCCCTGGGTGCAGGCATTCTCCACCACCTGAATCAGAAGGCAACAAAGAAAGGCGAGCAGGAACTTGCTGAGGTGAAAGTGGAAGCCATCCCCAACTGCAAGTCTGAGGAGGAAACCTCGCCCCTGGTGACACACCAGAACCCCGCTGGCCCCGTGGCCAGTGCCCCAGAACTGGAATCCAAGGAGTCGGTTCTGTGAAGAGCTTAGAGACCGACGACTGCCTAAGGACATTCGCTGCGTCTTCAGGCTAGGTGGAGGCTCAGTG", false}

	Clones, _, _ := em.GoldenGate([]Part{fragment1, fragment2, popen}, "BbsI")

	fmt.Println(seqhash.RotateSequence(Clones[0]))
	// Output: AAAAAAAGGATCTCAAGAAGGCCTACTATTAGCAACAACGATCCTTTGATCTTTTCTACGGGGTCTGACGCTCAGTGGAACGAAAACTCACGTTAAGGGATTTTGGTCATGAGATTATCAAAAAGGATCTTCACCTAGATCCTTTTAAATTAAAAATGAAGTTTTAAATCAATCTAAAGTATATATGAGTAAACTTGGTCTGACAGTTACCAATGCTTAATCAGTGAGGCACCTATCTCAGCGATCTGTCTATTTCGTTCATCCATAGTTGCCTGACTCCCCGTCGTGTAGATAACTACGATACGGGAGGGCTTACCATCTGGCCCCAGTGCTGCAATGATACCGCGAGAACCACGCTCACCGGCTCCAGATTTATCAGCAATAAACCAGCCAGCCGGAAGGGCCGAGCGCAGAAGTGGTCCTGCAACTTTATCCGCCTCCATCCAGTCTATTAATTGTTGCCGGGAAGCTAGAGTAAGTAGTTCGCCAGTTAATAGTTTGCGCAACGTTGTTGCCATTGCTACAGGCATCGTGGTGTCACGCTCGTCGTTTGGTATGGCTTCATTCAGCTCCGGTTCCCAACGATCAAGGCGAGTTACATGATCCCCCATGTTGTGCAAAAAAGCGGTTAGCTCCTTCGGTCCTCCGATCGTTGTCAGAAGTAAGTTGGCCGCAGTGTTATCACTCATGGTTATGGCAGCACTGCATAATTCTCTTACTGTCATGCCATCCGTAAGATGCTTTTCTGTGACTGGTGAGTACTCAACCAAGTCATTCTGAGAATAGTGTATGCGGCGACCGAGTTGCTCTTGCCCGGCGTCAATACGGGATAATACCGCGCCACATAGCAGAACTTTAAAAGTGCTCATCATTGGAAAACGTTCTTCGGGGCGAAAACTCTCAAGGATCTTACCGCTGTTGAGATCCAGTTCGATGTAACCCACTCGTGCACCCAACTGATCTTCAGCATCTTTTACTTTCACCAGCGTTTCTGGGTGAGCAAAAACAGGAAGGCAAAATGCCGCAAAAAAGGGAATAAGGGCGACACGGAAATGTTGAATACTCATACTCTTCCTTTTTCAATATTATTGAAGCATTTATCAGGGTTATTGTCTCATGAGCGGATACATATTTGAATGTATTTAGAAAAATAAACAAATAGGGGTTCCGCGCACCTGCACCAGTCAGTAAAACGACGGCCAGTAGTCAAAAGCCTCCGACCGGAGGCTTTTGACTTGGTTCAGGTGGAGTGGGAGAAACACGTGGCAAACATTCCGGTCTCAAATGGAAAAGAGCAACGAAACCAACGGCTACCTTGACAGCGCTCAAGCCGGCCCTGCAGCTGGCCCGGGCGCTCCGGGTACCGCCGCGGGTCGTGCACGTCGTTGCGCGGGCTTCCTGCGGCGCCAAGCGCTGGTGCTGCTCACGGTGTCTGGTGTTCTGGCAGGCGCCGGTTTGGGCGCGGCACTGCGTGGGCTCAGCCTGAGCCGCACCCAGGTCACCTACCTGGCCTTCCCCGGCGAGATGCTGCTCCGCATGCTGCGCATGATCATCCTGCCGCTGGTGGTCTGCAGCCTGGTGTCGGGCGCCGCCTCCCTCGATGCCAGCTGCCTCGGGCGTCTGGGCGGTATCGCTGTCGCCTACTTTGGCCTCACCACACTGAGTGCCTCGGCGCTCGCCGTGGCCTTGGCGTTCATCATCAAGCCAGGATCCGGTGCGCAGACCCTTCAGTCCAGCGACCTGGGGCTGGAGGACTCGGGGCCTCCTCCTGTCCCCAAAGAAACGGTGGACTCTTTCCTCGACCTGGCCAGAAACCTGTTTCCCTCCAATCTTGTGGTTGCAGCTTTCCGTACGTATGCAACCGATTATAAAGTCGTGACCCAGAACAGCAGCTCTGGAAATGTAACCCATGAAAAGATCCCCATAGGCACTGAGATAGAAGGGATGAACATTTTAGGATTGGTCCTGTTTGCTCTGGTGTTAGGAGTGGCCTTAAAGAAACTAGGCTCCGAAGGAGAGGACCTCATCCGTTTCTTCAATTCCCTCAACGAGGCGACGATGGTGCTGGTGTCCTGGATTATGTGGTACGTACCTGTGGGCATCATGTTCCTTGTTGGAAGCAAGATCGTGGAAATGAAAGACATCATCGTGCTGGTGACCAGCCTGGGGAAATACATCTTCGCATCTATATTGGGCCACGTCATTCATGGTGGTATCGTCCTGCCGCTGATTTATTTTGTTTTCACACGAAAAAACCCATTCAGATTCCTCCTGGGCCTCCTCGCCCCATTTGCGACAGCATTTGCTACGTGCTCCAGCTCAGCGACCCTTCCCTCTATGATGAAGTGCATTGAAGAGAACAATGGTGTGGACAAGAGGATCTCCAGGTTTATTCTCCCCATCGGGGCCACCGTGAACATGGACGGAGCAGCCATCTTCCAGTGTGTGGCCGCGGTGTTCATTGCGCAACTCAACAACGTAGAGCTCAACGCAGGACAGATTTTCACCATTCTAGTGACTGCCACAGCGTCCAGTGTTGGAGCAGCAGGCGTGCCAGCTGGAGGGGTCCTCACCATTGCCATTATCCTGGAGGCCATTGGGCTGCCTACTCATGATCTGCCTCTGATCCTGGCTGTGGACTGGATTGTGGACCGGACCACCACGGTGGTGAATGTGGAAGGGGATGCCCTGGGTGCAGGCATTCTCCACCACCTGAATCAGAAGGCAACAAAGAAAGGCGAGCAGGAACTTGCTGAGGTGAAAGTGGAAGCCATCCCCAACTGCAAGTCTGAGGAGGAAACCTCGCCCCTGGTGACACACCAGAACCCCGCTGGCCCCGTGGCCAGTGCCCCAGAACTGGAATCCAAGGAGTCGGTTCTGTGAAGAGCTTAGAGACCGACGACTGCCTAAGGACATTCGCTGAGGTGTCAATCGTCGGAGCCGCTGAGCAATAACTAGCATAACCCCTTGGGGCCTCTAAACGGGTCTTGAGGGGTTTTTTGCATGGTCATAGCTGTTTCCTGAGAGCTTGGCAGGTGATGACACACATTAACAAATTTCGTGAGGAGTCTCCAGAAGAATGCCATTAATTTCCATAGGCTCCGCCCCCCTGACGAGCATCACAAAAATCGACGCTCAAGTCAGAGGTGGCGAAACCCGACAGGACTATAAAGATACCAGGCGTTTCCCCCTGGAAGCTCCCTCGTGCGCTCTCCTGTTCCGACCCTGCCGCTTACCGGATACCTGTCCGCCTTTCTCCCTTCGGGAAGCGTGGCGCTTTCTCATAGCTCACGCTGTAGGTATCTCAGTTCGGTGTAGGTCGTTCGCTCCAAGCTGGGCTGTGTGCACGAACCCCCCGTTCAGCCCGACCGCTGCGCCTTATCCGGTAACTATCGTCTTGAGTCCAACCCGGTAAGACACGACTTATCGCCACTGGCAGCAGCCACTGGTAACAGGATTAGCAGAGCGAGGTATGTAGGCGGTGCTACAGAGTTCTTGAAGTGGTGGCCTAACTACGGCTACACTAGAAGAACAGTATTTGGTATCTGCGCTCTGCTGAAGCCAGTTACCTTCGGAAAAAGAGTTGGTAGCTCTTGATCCGGCAAACAAACCACCGCTGGTAGCGGTGGTTTTTTTGTTTGCAAGCAGCAGATTACGCGCAG
}

func TestSignalKilledGoldenGate(t *testing.T) {
	em := NewEnzymeManager(getBaseRestrictionEnzymes())
	// This previously would crash from using too much RAM.
	frag1 := Part{"AAAGCACTCTTAGGCCTCTGGAAGACATGGAGGGTCTCAAGGTGATCAAAGGATCTTCTTGAGATCCTTTTTTTCTGCGCGTAATCTTTTGCCCTGTAAACGAAAAAACCACCTGGGTAGTCTTCGCATTTCTTAATCGGTGCCC", false}
	frag2 := Part{"AAAGCACTCTTAGGCCTCTGGAAGACATTGGGGAGGTGGTTTGATCGAAGGTTAAGTCAGTTGGGGAACTGCTTAACCGTGGTAACTGGCTTTCGCAGAGCACAGCAACCAAATCTGTTAGTCTTCGCATTTCTTAATCGGTGCCC", false}
	frag3 := Part{"AAAGCACTCTTAGGCCTCTGGAAGACATCTGTCCTTCCAGTGTAGCCGGACTTTGGCGCACACTTCAAGAGCAACCGCGTGTTTAGCTAAACAAATCCTCTGCGAACTCCCAGTTACCTAGTCTTCGCATTTCTTAATCGGTGCCC", false}
	frag4 := Part{"AAAGCACTCTTAGGCCTCTGGAAGACATTACCAATGGCTGCTGCCAGTGGCGTTTTACCGTGCTTTTCCGGGTTGGACTCAAGTGAACAGTTACCGGATAAGGCGCAGCAGTCGGGCTTAGTCTTCGCATTTCTTAATCGGTGCCC", false}
	frag5 := Part{"AAAGCACTCTTAGGCCTCTGGAAGACATGGCTGAACGGGGAGTTCTTGCTTACAGCCCAGCTTGGAGCGAACGACCTACACCGAGCCGAGATACCAGTGTGTGAGCTATGAGAAAGCGTAGTCTTCGCATTTCTTAATCGGTGCCC", false}
	frag6 := Part{"AAAGCACTCTTAGGCCTCTGGAAGACATAGCGCCACACTTCCCGTAAGGGAGAAAGGCGGAACAGGTATCCGGTAAACGGCAGGGTCGGAACAGGAGAGCGCAAGAGGGAGCGACCCGTAGTCTTCGCATTTCTTAATCGGTGCCC", false}
	frag7 := Part{"AAAGCACTCTTAGGCCTCTGGAAGACATCCCGCCGGAAACGGTGGGGATCTTTAAGTCCTGTCGGGTTTCGCCCGTACTGTCAGATTCATGGTTGAGCCTCACGGCTCCCACAGATGTAGTCTTCGCATTTCTTAATCGGTGCCC", false}
	frag8 := Part{"AAAGCACTCTTAGGCCTCTGGAAGACATGATGCACCGGAAAAGCGTCTGTTTATGTGAACTCTGGCAGGAGGGCGGAGCCTATGGAAAAACGCCACCGGCGCGGCCCTGCTGTTTTGCCTCACATGTTAGTCTTCGCATTTCTTAATCGGTGCCC", false}
	frag9 := Part{"AAAGCACTCTTAGGCCTCTGGAAGACATATGTTAGTCCCCTGCTTATCCACGGAATCTGTGGGTAACTTTGTATGTGTCCGCAGCGCAAAAAGAGACCCGCTTAGTCTTCGCATTTCTTAATCGGTGCCC", false}
	fragments := []Part{popen, frag1, frag2, frag3, frag4, frag5, frag6, frag7, frag8, frag9}

	clones, loopingClones, err := em.GoldenGate(fragments, "BbsI")
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
	em := NewEnzymeManager(getBaseRestrictionEnzymes())
	// This used to panic with the message:
	// panic: runtime error: slice bounds out of range [:-2] [recovered]
	// It was from the following sequence: GAAGACATAATGGTCTTC . There are 2 intercepting BbsI sites.
	frag1 := Part{"AAACCGGAGCCATACAGTACGAAGACATGGAGGGTCTCAAATGAAAAAAATCATCGAAACCCAGCGTGCACCGGGAGCAATCGGACCGTACGTCCAGGGAGTCGACCTAGGATCAATGTAGTCTTCGCACTTGGCTTAGATGCAAC", false}
	frag2 := Part{"AAACCGGAGCCATACAGTACGAAGACATAATGGTCTTCACCTCAGGACAGATCCCGGTCTGCCCGCAGACCGGAGAAATCCCGGCAGACGTCCAGGACCAGGCACGTCTATCACTAGATAGTCTTCGCACTTGGCTTAGATGCAAC", false}
	frag3 := Part{"AAACCGGAGCCATACAGTACGAAGACATTAGAAAACGTCAAAGCAATCGTCGTCGCAGCAGGACTATCAGTCGGAGACATCATCAAAATGACCGTCTTCATCACCGACCTAAACGACTTAGTCTTCGCACTTGGCTTAGATGCAAC", false}
	frag4 := Part{"AAACCGGAGCCATACAGTACGAAGACATGACTTCGCAACCATCAACGAAGTCTACAAACAGTTCTTCGACGAACACCAGGCAACCTACCCGACCCGTTCATGCGTCCAGGTCGCACGTCTACTAGTCTTCGCACTTGGCTTAGATGCAAC", false}
	frag5 := Part{"AAACCGGAGCCATACAGTACGAAGACATCTACCGAAAGACGTCAAACTAGAAATCGAAGCAATCGCAGTCCGTTCAGCAAGAGCTTAGAGACCCGCTTAGTCTTCGCACTTGGCTTAGATGCAAC", false}
	fragments := []Part{popen, frag1, frag2, frag3, frag4, frag5}

	_, _, err := em.GoldenGate(fragments, "BbsI")
	if err != nil {
		t.Errorf("GoldenGate should not fail with these fragments. Got error: %s", err)
	}
}

func TestCircularCutRegression(t *testing.T) {
	em := NewEnzymeManager(getBaseRestrictionEnzymes())
	// This used to error with 0 fragments since the BsaI cut site is on the other
	// side of the origin from its recognition site.
	plasmid1 := Part{"AAACTACAAGACCCGCGCCGAGGTGAAGTTCGAGGGCGACACCCTGGTGAACCGCATCGAGCTGAAGGGCATCGACTTCAAGGAGGACGGCAACATCCTGGGGCACAAGCTGGAGTACAACTACAACAGCCACAACGTCTATATCATGGCCGACAAGCAGAAGAACGGCATCAAGGTGAACTTCAAGATCCGCCACAACATCGAGGACGGCAGCCGAGaccaagtcgcggccgcgaggtgtcaatcgtcggagtagggataacagggtaatccgctgagcaataactagcataaccccttggggcctctaaacgggtcttgaggggttttttgcatggtcatagctgtttcctgttacgccccgccctgccactcgtcgcagtactgttgtaattcattaagcattctgccgacatggaagccatcacaaacggcatgatgaacctgaatcgccagcggcatcagcaccttgtcgccttgcgtataatatttgcccatggtgaaaacgggggcgaagaagttgtccatattggccacgtttaaatcaaaactggtgaaactcacccagggattggctgacacgaaaaacatattctcaataaaccctttagggaaataggccaggttttcaccgtaacacgccacatcttgcgaatatatgtgtagaaactgccggaaatcgtcgtggtattcactccagagggatgaaaacgtttcagtttgctcatggaaaacggtgtaacaagggtgaacactatcccatatcaccagctcaccatccttcattgccatacgaaattccggatgagcattcatcaggcgggcaagaatgtgaataaaggccggataaaacttgtgcttatttttctttacggtctttaaaaaggccgtaatatccagctgaacggtctggttataggtacattgagcaactgactgaaatgcctcaaaatgttctttacgatgccattgggatatatcaacggtggtatatccagtgatttttttctccattttagcttccttagctcctgaaaatctcgataactcaaaaaatacgcccggtagtgatcttatttcattatggtgaaagttggaacctcttacgtgccgatcatttccataggctccgcccccctgacgagcatcacaaaaatcgacgctcaagtcagaggtggcgaaacccgacaggactataaagataccaggcgtttccccctggaagctccctcgtgcgctctcctgttccgaccctgccgcttaccggatacctgtccgcctttctcccttcgggaagcgtggcgctttctcatagctcacgctgtaggtatctcagttcggtgtaggtcgttcgctccaagctgggctgtgtgcacgaaccccccgttcagcccgaccgctgcgccttatccggtaactatcgtcttgagtccaacccggtaagacacgacttatcgccactggcagcagccactggtaacaggattagcagagcgaggtatgtaggcggtgctacagagttcttgaagtggtggcctaactacggctacactagaaggacagtatttggtatctgcgctctgctgaagccagttaccttcggaaaaagagttggtagctcttgatccggcaaacaaaccaccgctggtagcggtggtttttttgtttgcaagcagcagattacgcgcagaaaaaaaggatctcaagtaaaacgacggccagtagtcaaaagcctccgaccggaggcttttgacttggttcaggtggagtggcggccgcgacttgGTCTC", true}
	newFragments, err := em.CutWithEnzymeByName(plasmid1, true, "BsaI")
	if err != nil {
		t.Errorf("Failed to cut: %s", err)
	}
	if len(newFragments) != 1 {
		t.Errorf("Expected 1 new fragment, got: %d", len(newFragments))
	}
}

func benchmarkGoldenGate(b *testing.B, em EnzymeManager, parts []Part) {
	for n := 0; n < b.N; n++ {
		em.GoldenGate(parts, "BbsI")
	}
}

func BenchmarkGoldenGate3Parts(b *testing.B) {
	em := NewEnzymeManager(getBaseRestrictionEnzymes())
	fragment1 := Part{"GAAGTGCCATTCCGCCTGACCTGAAGACCAGGAGAAACACGTGGCAAACATTCCGGTCTCAAATGGAAAAGAGCAACGAAACCAACGGCTACCTTGACAGCGCTCAAGCCGGCCCTGCAGCTGGCCCGGGCGCTCCGGGTACCGCCGCGGGTCGTGCACGTCGTTGCGCGGGCTTCCTGCGGCGCCAAGCGCTGGTGCTGCTCACGGTGTCTGGTGTTCTGGCAGGCGCCGGTTTGGGCGCGGCACTGCGTGGGCTCAGCCTGAGCCGCACCCAGGTCACCTACCTGGCCTTCCCCGGCGAGATGCTGCTCCGCATGCTGCGCATGATCATCCTGCCGCTGGTGGTCTGCAGCCTGGTGTCGGGCGCCGCCTCCCTCGATGCCAGCTGCCTCGGGCGTCTGGGCGGTATCGCTGTCGCCTACTTTGGCCTCACCACACTGAGTGCCTCGGCGCTCGCCGTGGCCTTGGCGTTCATCATCAAGCCAGGATCCGGTGCGCAGACCCTTCAGTCCAGCGACCTGGGGCTGGAGGACTCGGGGCCTCCTCCTGTCCCCAAAGAAACGGTGGACTCTTTCCTCGACCTGGCCAGAAACCTGTTTCCCTCCAATCTTGTGGTTGCAGCTTTCCGTACGTATGCAACCGATTATAAAGTCGTGACCCAGAACAGCAGCTCTGGAAATGTAACCCATGAAAAGATCCCCATAGGCACTGAGATAGAAGGGATGAACATTTTAGGATTGGTCCTGTTTGCTCTGGTGTTAGGAGTGGCCTTAAAGAAACTAGGCTCCGAAGGAGAGGACCTCATCCGTTTCTTCAATTCCCTCAACGAGGCGACGATGGTGCTGGTGTCCTGGATTATGTGGTACGCGTCTTCAGGCTAGGTGGAGGCTCAGTG", false}
	fragment2 := Part{"GAAGTGCCATTCCGCCTGACCTGAAGACCAGTACGTACCTGTGGGCATCATGTTCCTTGTTGGAAGCAAGATCGTGGAAATGAAAGACATCATCGTGCTGGTGACCAGCCTGGGGAAATACATCTTCGCATCTATATTGGGCCACGTCATTCATGGTGGTATCGTCCTGCCGCTGATTTATTTTGTTTTCACACGAAAAAACCCATTCAGATTCCTCCTGGGCCTCCTCGCCCCATTTGCGACAGCATTTGCTACGTGCTCCAGCTCAGCGACCCTTCCCTCTATGATGAAGTGCATTGAAGAGAACAATGGTGTGGACAAGAGGATCTCCAGGTTTATTCTCCCCATCGGGGCCACCGTGAACATGGACGGAGCAGCCATCTTCCAGTGTGTGGCCGCGGTGTTCATTGCGCAACTCAACAACGTAGAGCTCAACGCAGGACAGATTTTCACCATTCTAGTGACTGCCACAGCGTCCAGTGTTGGAGCAGCAGGCGTGCCAGCTGGAGGGGTCCTCACCATTGCCATTATCCTGGAGGCCATTGGGCTGCCTACTCATGATCTGCCTCTGATCCTGGCTGTGGACTGGATTGTGGACCGGACCACCACGGTGGTGAATGTGGAAGGGGATGCCCTGGGTGCAGGCATTCTCCACCACCTGAATCAGAAGGCAACAAAGAAAGGCGAGCAGGAACTTGCTGAGGTGAAAGTGGAAGCCATCCCCAACTGCAAGTCTGAGGAGGAAACCTCGCCCCTGGTGACACACCAGAACCCCGCTGGCCCCGTGGCCAGTGCCCCAGAACTGGAATCCAAGGAGTCGGTTCTGTGAAGAGCTTAGAGACCGACGACTGCCTAAGGACATTCGCTGCGTCTTCAGGCTAGGTGGAGGCTCAGTG", false}

	benchmarkGoldenGate(b, em, []Part{fragment1, fragment2, popen})
}
