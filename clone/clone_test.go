package clone

import (
	"fmt"
	"log"
	"testing"

	"github.com/Open-Science-Global/poly/seqhash"
)

// pOpen plasmid series (https://stanford.freegenes.org/collections/open-genes/products/open-plasmids#description). I use it for essentially all my cloning. -Keoni
var popen = Part{"TAACTATCGTCTTGAGTCCAACCCGGTAAGACACGACTTATCGCCACTGGCAGCAGCCACTGGTAACAGGATTAGCAGAGCGAGGTATGTAGGCGGTGCTACAGAGTTCTTGAAGTGGTGGCCTAACTACGGCTACACTAGAAGAACAGTATTTGGTATCTGCGCTCTGCTGAAGCCAGTTACCTTCGGAAAAAGAGTTGGTAGCTCTTGATCCGGCAAACAAACCACCGCTGGTAGCGGTGGTTTTTTTGTTTGCAAGCAGCAGATTACGCGCAGAAAAAAAGGATCTCAAGAAGGCCTACTATTAGCAACAACGATCCTTTGATCTTTTCTACGGGGTCTGACGCTCAGTGGAACGAAAACTCACGTTAAGGGATTTTGGTCATGAGATTATCAAAAAGGATCTTCACCTAGATCCTTTTAAATTAAAAATGAAGTTTTAAATCAATCTAAAGTATATATGAGTAAACTTGGTCTGACAGTTACCAATGCTTAATCAGTGAGGCACCTATCTCAGCGATCTGTCTATTTCGTTCATCCATAGTTGCCTGACTCCCCGTCGTGTAGATAACTACGATACGGGAGGGCTTACCATCTGGCCCCAGTGCTGCAATGATACCGCGAGAACCACGCTCACCGGCTCCAGATTTATCAGCAATAAACCAGCCAGCCGGAAGGGCCGAGCGCAGAAGTGGTCCTGCAACTTTATCCGCCTCCATCCAGTCTATTAATTGTTGCCGGGAAGCTAGAGTAAGTAGTTCGCCAGTTAATAGTTTGCGCAACGTTGTTGCCATTGCTACAGGCATCGTGGTGTCACGCTCGTCGTTTGGTATGGCTTCATTCAGCTCCGGTTCCCAACGATCAAGGCGAGTTACATGATCCCCCATGTTGTGCAAAAAAGCGGTTAGCTCCTTCGGTCCTCCGATCGTTGTCAGAAGTAAGTTGGCCGCAGTGTTATCACTCATGGTTATGGCAGCACTGCATAATTCTCTTACTGTCATGCCATCCGTAAGATGCTTTTCTGTGACTGGTGAGTACTCAACCAAGTCATTCTGAGAATAGTGTATGCGGCGACCGAGTTGCTCTTGCCCGGCGTCAATACGGGATAATACCGCGCCACATAGCAGAACTTTAAAAGTGCTCATCATTGGAAAACGTTCTTCGGGGCGAAAACTCTCAAGGATCTTACCGCTGTTGAGATCCAGTTCGATGTAACCCACTCGTGCACCCAACTGATCTTCAGCATCTTTTACTTTCACCAGCGTTTCTGGGTGAGCAAAAACAGGAAGGCAAAATGCCGCAAAAAAGGGAATAAGGGCGACACGGAAATGTTGAATACTCATACTCTTCCTTTTTCAATATTATTGAAGCATTTATCAGGGTTATTGTCTCATGAGCGGATACATATTTGAATGTATTTAGAAAAATAAACAAATAGGGGTTCCGCGCACCTGCACCAGTCAGTAAAACGACGGCCAGTAGTCAAAAGCCTCCGACCGGAGGCTTTTGACTTGGTTCAGGTGGAGTGGGAGTAgtcttcGCcatcgCtACTAAAagccagataacagtatgcgtatttgcgcgctgatttttgcggtataagaatatatactgatatgtatacccgaagtatgtcaaaaagaggtatgctatgaagcagcgtattacagtgacagttgacagcgacagctatcagttgctcaaggcatatatgatgtcaatatctccggtctggtaagcacaaccatgcagaatgaagcccgtcgtctgcgtgccgaacgctggaaagcggaaaatcaggaagggatggctgaggtcgcccggtttattgaaatgaacggctcttttgctgacgagaacagggGCTGGTGAAATGCAGTTTAAGGTTTACACCTATAAAAGAGAGAGCCGTTATCGTCTGTTTGTGGATGTACAGAGTGATATTATTGACACGCCCGGGCGACGGATGGTGATCCCCCTGGCCAGTGCACGTCTGCTGTCAGATAAAGTCTCCCGTGAACTTTACCCGGTGGTGCATATCGGGGATGAAAGCTGGCGCATGATGACCACCGATATGGCCAGTGTGCCGGTCTCCGTTATCGGGGAAGAAGTGGCTGATCTCAGCCACCGCGAAAATGACATCAAAAACGCCATTAACCTGATGTTCTGGGGAATATAAATGTCAGGCTCCCTTATACACAGgcgatgttgaagaccaCGCTGAGGTGTCAATCGTCGGAGCCGCTGAGCAATAACTAGCATAACCCCTTGGGGCCTCTAAACGGGTCTTGAGGGGTTTTTTGCATGGTCATAGCTGTTTCCTGAGAGCTTGGCAGGTGATGACACACATTAACAAATTTCGTGAGGAGTCTCCAGAAGAATGCCATTAATTTCCATAGGCTCCGCCCCCCTGACGAGCATCACAAAAATCGACGCTCAAGTCAGAGGTGGCGAAACCCGACAGGACTATAAAGATACCAGGCGTTTCCCCCTGGAAGCTCCCTCGTGCGCTCTCCTGTTCCGACCCTGCCGCTTACCGGATACCTGTCCGCCTTTCTCCCTTCGGGAAGCGTGGCGCTTTCTCATAGCTCACGCTGTAGGTATCTCAGTTCGGTGTAGGTCGTTCGCTCCAAGCTGGGCTGTGTGCACGAACCCCCCGTTCAGCCCGACCGCTGCGCCTTATCCGG", true}

func TestCutWithEnzymeByName(t *testing.T) {
	_, err := CutWithEnzymeByName(popen, true, "EcoFake")
	if err == nil {
		log.Fatalf("CutWithEnzymeByName should have failed when looking for fake restriction enzyme EcoFake")
	}
}

func TestCutWithEnzyme(t *testing.T) {
	var seq Part
	bsai := "GGTCTCAATGC"
	bsaiComplement := "ATGCAGAGACC"

	// test(1)
	// Test case of `<-bsaiComplement bsai-> <-bsaiComplement bsai->` where bsaI cuts off of a linear sequence. This tests the line:
	// if !seq.Circular && (overhangSet[len(overhangSet)-1].Position+enzyme.EnzymeSkip+enzyme.EnzymeOverhangLen > len(sequence))
	seq = Part{"ATATATA" + bsaiComplement + bsai + "ATGCATCGATCGACTAGCATG" + bsaiComplement + bsai[:8], false}
	frag, err := CutWithEnzymeByName(seq, true, "BsaI")
	if err != nil {
		log.Fatalf("CutWithEnzyme should not have failed on test(1). Got error: %s", err)
	}
	if len(frag) != 1 {
		log.Fatalf("CutWithEnzyme in test(1) should be 1 fragment in length")
	}
	if frag[0].Sequence != "ATGCATCGATCGACTAGCATG" {
		log.Fatalf("CutWithEnzyme in test(1) should give fragment with sequence ATGCATCGATCGACTAGCATG . Got sequence: %s", frag[0].Sequence)
	}

	// test(2)
	// Now if we take the same sequence and circularize it, we get a different result
	seq.Circular = true
	frag, err = CutWithEnzymeByName(seq, true, "BsaI")
	if err != nil {
		log.Fatalf("CutWithEnzyme should not have failed on test(2). Got error: %s", err)
	}
	if len(frag) != 2 {
		log.Fatalf("CutWithEnzyme in test(2) should be 1 fragment in length")
	}
	if frag[0].Sequence != "ATGCATCGATCGACTAGCATG" || frag[1].Sequence != "TATA" {
		log.Fatalf("CutWithEnzyme in test(2) should give fragment with sequence ATGCATCGATCGACTAGCATG and TATA. Got sequence: %s and %s", frag[0].Sequence, frag[1].Sequence)
	}

	// test(3)
	// Let's test if we have a single cut in our plasmid. This should give
	// different results if we have a linear or circular DNA. Since single cuts
	// will give no fragments if you test for directionality, we set the
	// directionality flag to false. This tests the line:
	// if len(overhangs) == 1 && !directional && !seq.Circular
	seq = Part{"ATATATATATATATAT" + bsai + "GCGCGCGCGCGCGCGCGCGC", false}
	frag, err = CutWithEnzymeByName(seq, false, "BsaI")
	if err != nil {
		log.Fatalf("CutWithEnzyme should not have failed on test(3). Got error: %s", err)
	}
	if len(frag) != 2 {
		log.Fatalf("Cutting a linear fragment with a single cut site should give 2 fragments")
	}
	if frag[0].Sequence != "GCGCGCGCGCGCGCGCGCGC" || frag[1].Sequence != "ATATATATATATATATGGTCTCA" {
		log.Fatalf("CutWithEnzyme in test(3) should give fragment with sequence GCGCGCGCGCGCGCGCGCGC and ATATATATATATATATGGTCTCA. Got sequence: %s and %s", frag[0].Sequence, frag[1].Sequence)
	}

	// test(4)
	// This tests for the above except with a circular fragment. Specifically, it
	// tests the line:
	// if len(overhangs) == 2 && !directional && seq.Circular
	seq.Circular = true
	frag, err = CutWithEnzymeByName(seq, false, "BsaI")
	if err != nil {
		log.Fatalf("CutWithEnzyme should not have failed on test(4). Got error: %s", err)
	}
	if len(frag) != 1 {
		log.Fatalf("Cutting a circular fragment with a single cut site should give 1 fragments")
	}
	if frag[0].Sequence != "GCGCGCGCGCGCGCGCGCGCATATATATATATATATGGTCTCA" {
		log.Fatalf("CutWithEnzyme in test(4) should give fragment with sequence ATATATATATATATATGGTCTCA. Got Sequence: %s", frag[0].Sequence)
	}

	// test(5)
	// This tests if we have a fragment where we do not care about directionality
	// but have more than 1 cut site in our fragment. We can use pOpen for this.
	frag, err = CutWithEnzymeByName(popen, false, "BbsI")
	if err != nil {
		log.Fatalf("CutWithEnzyme should not have failed on test(5). Got error: %s", err)
	}
	if len(frag) != 2 {
		log.Fatalf("Cutting pOpen without a direction should yield 2 fragments")
	}
}

func TestCircularLigate(t *testing.T) {
	// The following tests for complementing overhangs. Specific, this line:
	// newSeed := Fragment{seedFragment.Sequence + seedFragment.ReverseOverhang + ReverseComplement(newFragment.Sequence), seedFragment.ForwardOverhang, ReverseComplement(newFragment.ForwardOverhang)}
	fragment1 := Fragment{"AAAAAA", "GTTG", "CTAT"}
	fragment2 := Fragment{"AAAAAA", "CAAC", "ATAG"}
	outputConstructs := CircularLigate([]Fragment{fragment1, fragment2})
	if len(outputConstructs) != 1 {
		fmt.Println(outputConstructs)
		log.Fatalf("Circular ligation with complementing overhangs should only output 1 valid rotated sequence.")
	}
}

func TestGoldenGate(t *testing.T) {
	// Here we test if the enzyme we want to use in a GoldenGate reaction does not exist in our enzyme pool
	fragment1 := Part{"GAAGTGCCATTCCGCCTGACCTGAAGACCAGGAGAAACACGTGGCAAACATTCCGGTCTCAAATGGAAAAGAGCAACGAAACCAACGGCTACCTTGACAGCGCTCAAGCCGGCCCTGCAGCTGGCCCGGGCGCTCCGGGTACCGCCGCGGGTCGTGCACGTCGTTGCGCGGGCTTCCTGCGGCGCCAAGCGCTGGTGCTGCTCACGGTGTCTGGTGTTCTGGCAGGCGCCGGTTTGGGCGCGGCACTGCGTGGGCTCAGCCTGAGCCGCACCCAGGTCACCTACCTGGCCTTCCCCGGCGAGATGCTGCTCCGCATGCTGCGCATGATCATCCTGCCGCTGGTGGTCTGCAGCCTGGTGTCGGGCGCCGCCTCCCTCGATGCCAGCTGCCTCGGGCGTCTGGGCGGTATCGCTGTCGCCTACTTTGGCCTCACCACACTGAGTGCCTCGGCGCTCGCCGTGGCCTTGGCGTTCATCATCAAGCCAGGATCCGGTGCGCAGACCCTTCAGTCCAGCGACCTGGGGCTGGAGGACTCGGGGCCTCCTCCTGTCCCCAAAGAAACGGTGGACTCTTTCCTCGACCTGGCCAGAAACCTGTTTCCCTCCAATCTTGTGGTTGCAGCTTTCCGTACGTATGCAACCGATTATAAAGTCGTGACCCAGAACAGCAGCTCTGGAAATGTAACCCATGAAAAGATCCCCATAGGCACTGAGATAGAAGGGATGAACATTTTAGGATTGGTCCTGTTTGCTCTGGTGTTAGGAGTGGCCTTAAAGAAACTAGGCTCCGAAGGAGAGGACCTCATCCGTTTCTTCAATTCCCTCAACGAGGCGACGATGGTGCTGGTGTCCTGGATTATGTGGTACGCGTCTTCAGGCTAGGTGGAGGCTCAGTG", false}
	fragment2 := Part{"GAAGTGCCATTCCGCCTGACCTGAAGACCAGTACGTACCTGTGGGCATCATGTTCCTTGTTGGAAGCAAGATCGTGGAAATGAAAGACATCATCGTGCTGGTGACCAGCCTGGGGAAATACATCTTCGCATCTATATTGGGCCACGTCATTCATGGTGGTATCGTCCTGCCGCTGATTTATTTTGTTTTCACACGAAAAAACCCATTCAGATTCCTCCTGGGCCTCCTCGCCCCATTTGCGACAGCATTTGCTACGTGCTCCAGCTCAGCGACCCTTCCCTCTATGATGAAGTGCATTGAAGAGAACAATGGTGTGGACAAGAGGATCTCCAGGTTTATTCTCCCCATCGGGGCCACCGTGAACATGGACGGAGCAGCCATCTTCCAGTGTGTGGCCGCGGTGTTCATTGCGCAACTCAACAACGTAGAGCTCAACGCAGGACAGATTTTCACCATTCTAGTGACTGCCACAGCGTCCAGTGTTGGAGCAGCAGGCGTGCCAGCTGGAGGGGTCCTCACCATTGCCATTATCCTGGAGGCCATTGGGCTGCCTACTCATGATCTGCCTCTGATCCTGGCTGTGGACTGGATTGTGGACCGGACCACCACGGTGGTGAATGTGGAAGGGGATGCCCTGGGTGCAGGCATTCTCCACCACCTGAATCAGAAGGCAACAAAGAAAGGCGAGCAGGAACTTGCTGAGGTGAAAGTGGAAGCCATCCCCAACTGCAAGTCTGAGGAGGAAACCTCGCCCCTGGTGACACACCAGAACCCCGCTGGCCCCGTGGCCAGTGCCCCAGAACTGGAATCCAAGGAGTCGGTTCTGTGAAGAGCTTAGAGACCGACGACTGCCTAAGGACATTCGCTGCGTCTTCAGGCTAGGTGGAGGCTCAGTG", false}

	_, err := GoldenGate([]Part{fragment1, fragment2, popen}, "EcoRFake")
	if err == nil {
		log.Fatalf("GoldenGate should fail when using enzyme EcoRFake")
	}
	if err.Error() != "Enzyme EcoRFake not found in enzymeMap" {
		log.Fatalf("Failure of GoldenGate on incorrect enzyme should follow the exact string `Enzyme EcoRFake not found in enzymeMap`. Got: %s", err.Error())
	}
}

func ExampleGoldenGate() {
	// Fragment 1 has a palindrome at its start. This isn't very common but
	// can occur. These two fragments are real DNA fragments used in the
	// FreeGenes Project. They are used because they were on my computer
	// - Keoni
	fragment1 := Part{"GAAGTGCCATTCCGCCTGACCTGAAGACCAGGAGAAACACGTGGCAAACATTCCGGTCTCAAATGGAAAAGAGCAACGAAACCAACGGCTACCTTGACAGCGCTCAAGCCGGCCCTGCAGCTGGCCCGGGCGCTCCGGGTACCGCCGCGGGTCGTGCACGTCGTTGCGCGGGCTTCCTGCGGCGCCAAGCGCTGGTGCTGCTCACGGTGTCTGGTGTTCTGGCAGGCGCCGGTTTGGGCGCGGCACTGCGTGGGCTCAGCCTGAGCCGCACCCAGGTCACCTACCTGGCCTTCCCCGGCGAGATGCTGCTCCGCATGCTGCGCATGATCATCCTGCCGCTGGTGGTCTGCAGCCTGGTGTCGGGCGCCGCCTCCCTCGATGCCAGCTGCCTCGGGCGTCTGGGCGGTATCGCTGTCGCCTACTTTGGCCTCACCACACTGAGTGCCTCGGCGCTCGCCGTGGCCTTGGCGTTCATCATCAAGCCAGGATCCGGTGCGCAGACCCTTCAGTCCAGCGACCTGGGGCTGGAGGACTCGGGGCCTCCTCCTGTCCCCAAAGAAACGGTGGACTCTTTCCTCGACCTGGCCAGAAACCTGTTTCCCTCCAATCTTGTGGTTGCAGCTTTCCGTACGTATGCAACCGATTATAAAGTCGTGACCCAGAACAGCAGCTCTGGAAATGTAACCCATGAAAAGATCCCCATAGGCACTGAGATAGAAGGGATGAACATTTTAGGATTGGTCCTGTTTGCTCTGGTGTTAGGAGTGGCCTTAAAGAAACTAGGCTCCGAAGGAGAGGACCTCATCCGTTTCTTCAATTCCCTCAACGAGGCGACGATGGTGCTGGTGTCCTGGATTATGTGGTACGCGTCTTCAGGCTAGGTGGAGGCTCAGTG", false}
	fragment2 := Part{"GAAGTGCCATTCCGCCTGACCTGAAGACCAGTACGTACCTGTGGGCATCATGTTCCTTGTTGGAAGCAAGATCGTGGAAATGAAAGACATCATCGTGCTGGTGACCAGCCTGGGGAAATACATCTTCGCATCTATATTGGGCCACGTCATTCATGGTGGTATCGTCCTGCCGCTGATTTATTTTGTTTTCACACGAAAAAACCCATTCAGATTCCTCCTGGGCCTCCTCGCCCCATTTGCGACAGCATTTGCTACGTGCTCCAGCTCAGCGACCCTTCCCTCTATGATGAAGTGCATTGAAGAGAACAATGGTGTGGACAAGAGGATCTCCAGGTTTATTCTCCCCATCGGGGCCACCGTGAACATGGACGGAGCAGCCATCTTCCAGTGTGTGGCCGCGGTGTTCATTGCGCAACTCAACAACGTAGAGCTCAACGCAGGACAGATTTTCACCATTCTAGTGACTGCCACAGCGTCCAGTGTTGGAGCAGCAGGCGTGCCAGCTGGAGGGGTCCTCACCATTGCCATTATCCTGGAGGCCATTGGGCTGCCTACTCATGATCTGCCTCTGATCCTGGCTGTGGACTGGATTGTGGACCGGACCACCACGGTGGTGAATGTGGAAGGGGATGCCCTGGGTGCAGGCATTCTCCACCACCTGAATCAGAAGGCAACAAAGAAAGGCGAGCAGGAACTTGCTGAGGTGAAAGTGGAAGCCATCCCCAACTGCAAGTCTGAGGAGGAAACCTCGCCCCTGGTGACACACCAGAACCCCGCTGGCCCCGTGGCCAGTGCCCCAGAACTGGAATCCAAGGAGTCGGTTCTGTGAAGAGCTTAGAGACCGACGACTGCCTAAGGACATTCGCTGCGTCTTCAGGCTAGGTGGAGGCTCAGTG", false}

	Clones, _ := GoldenGate([]Part{fragment1, fragment2, popen}, "BbsI")

	fmt.Println(seqhash.RotateSequence(Clones[0].Sequence))
	// Output: AAAAAAAGGATCTCAAGAAGGCCTACTATTAGCAACAACGATCCTTTGATCTTTTCTACGGGGTCTGACGCTCAGTGGAACGAAAACTCACGTTAAGGGATTTTGGTCATGAGATTATCAAAAAGGATCTTCACCTAGATCCTTTTAAATTAAAAATGAAGTTTTAAATCAATCTAAAGTATATATGAGTAAACTTGGTCTGACAGTTACCAATGCTTAATCAGTGAGGCACCTATCTCAGCGATCTGTCTATTTCGTTCATCCATAGTTGCCTGACTCCCCGTCGTGTAGATAACTACGATACGGGAGGGCTTACCATCTGGCCCCAGTGCTGCAATGATACCGCGAGAACCACGCTCACCGGCTCCAGATTTATCAGCAATAAACCAGCCAGCCGGAAGGGCCGAGCGCAGAAGTGGTCCTGCAACTTTATCCGCCTCCATCCAGTCTATTAATTGTTGCCGGGAAGCTAGAGTAAGTAGTTCGCCAGTTAATAGTTTGCGCAACGTTGTTGCCATTGCTACAGGCATCGTGGTGTCACGCTCGTCGTTTGGTATGGCTTCATTCAGCTCCGGTTCCCAACGATCAAGGCGAGTTACATGATCCCCCATGTTGTGCAAAAAAGCGGTTAGCTCCTTCGGTCCTCCGATCGTTGTCAGAAGTAAGTTGGCCGCAGTGTTATCACTCATGGTTATGGCAGCACTGCATAATTCTCTTACTGTCATGCCATCCGTAAGATGCTTTTCTGTGACTGGTGAGTACTCAACCAAGTCATTCTGAGAATAGTGTATGCGGCGACCGAGTTGCTCTTGCCCGGCGTCAATACGGGATAATACCGCGCCACATAGCAGAACTTTAAAAGTGCTCATCATTGGAAAACGTTCTTCGGGGCGAAAACTCTCAAGGATCTTACCGCTGTTGAGATCCAGTTCGATGTAACCCACTCGTGCACCCAACTGATCTTCAGCATCTTTTACTTTCACCAGCGTTTCTGGGTGAGCAAAAACAGGAAGGCAAAATGCCGCAAAAAAGGGAATAAGGGCGACACGGAAATGTTGAATACTCATACTCTTCCTTTTTCAATATTATTGAAGCATTTATCAGGGTTATTGTCTCATGAGCGGATACATATTTGAATGTATTTAGAAAAATAAACAAATAGGGGTTCCGCGCACCTGCACCAGTCAGTAAAACGACGGCCAGTAGTCAAAAGCCTCCGACCGGAGGCTTTTGACTTGGTTCAGGTGGAGTGGGAGAAACACGTGGCAAACATTCCGGTCTCAAATGGAAAAGAGCAACGAAACCAACGGCTACCTTGACAGCGCTCAAGCCGGCCCTGCAGCTGGCCCGGGCGCTCCGGGTACCGCCGCGGGTCGTGCACGTCGTTGCGCGGGCTTCCTGCGGCGCCAAGCGCTGGTGCTGCTCACGGTGTCTGGTGTTCTGGCAGGCGCCGGTTTGGGCGCGGCACTGCGTGGGCTCAGCCTGAGCCGCACCCAGGTCACCTACCTGGCCTTCCCCGGCGAGATGCTGCTCCGCATGCTGCGCATGATCATCCTGCCGCTGGTGGTCTGCAGCCTGGTGTCGGGCGCCGCCTCCCTCGATGCCAGCTGCCTCGGGCGTCTGGGCGGTATCGCTGTCGCCTACTTTGGCCTCACCACACTGAGTGCCTCGGCGCTCGCCGTGGCCTTGGCGTTCATCATCAAGCCAGGATCCGGTGCGCAGACCCTTCAGTCCAGCGACCTGGGGCTGGAGGACTCGGGGCCTCCTCCTGTCCCCAAAGAAACGGTGGACTCTTTCCTCGACCTGGCCAGAAACCTGTTTCCCTCCAATCTTGTGGTTGCAGCTTTCCGTACGTATGCAACCGATTATAAAGTCGTGACCCAGAACAGCAGCTCTGGAAATGTAACCCATGAAAAGATCCCCATAGGCACTGAGATAGAAGGGATGAACATTTTAGGATTGGTCCTGTTTGCTCTGGTGTTAGGAGTGGCCTTAAAGAAACTAGGCTCCGAAGGAGAGGACCTCATCCGTTTCTTCAATTCCCTCAACGAGGCGACGATGGTGCTGGTGTCCTGGATTATGTGGTACGTACCTGTGGGCATCATGTTCCTTGTTGGAAGCAAGATCGTGGAAATGAAAGACATCATCGTGCTGGTGACCAGCCTGGGGAAATACATCTTCGCATCTATATTGGGCCACGTCATTCATGGTGGTATCGTCCTGCCGCTGATTTATTTTGTTTTCACACGAAAAAACCCATTCAGATTCCTCCTGGGCCTCCTCGCCCCATTTGCGACAGCATTTGCTACGTGCTCCAGCTCAGCGACCCTTCCCTCTATGATGAAGTGCATTGAAGAGAACAATGGTGTGGACAAGAGGATCTCCAGGTTTATTCTCCCCATCGGGGCCACCGTGAACATGGACGGAGCAGCCATCTTCCAGTGTGTGGCCGCGGTGTTCATTGCGCAACTCAACAACGTAGAGCTCAACGCAGGACAGATTTTCACCATTCTAGTGACTGCCACAGCGTCCAGTGTTGGAGCAGCAGGCGTGCCAGCTGGAGGGGTCCTCACCATTGCCATTATCCTGGAGGCCATTGGGCTGCCTACTCATGATCTGCCTCTGATCCTGGCTGTGGACTGGATTGTGGACCGGACCACCACGGTGGTGAATGTGGAAGGGGATGCCCTGGGTGCAGGCATTCTCCACCACCTGAATCAGAAGGCAACAAAGAAAGGCGAGCAGGAACTTGCTGAGGTGAAAGTGGAAGCCATCCCCAACTGCAAGTCTGAGGAGGAAACCTCGCCCCTGGTGACACACCAGAACCCCGCTGGCCCCGTGGCCAGTGCCCCAGAACTGGAATCCAAGGAGTCGGTTCTGTGAAGAGCTTAGAGACCGACGACTGCCTAAGGACATTCGCTGAGGTGTCAATCGTCGGAGCCGCTGAGCAATAACTAGCATAACCCCTTGGGGCCTCTAAACGGGTCTTGAGGGGTTTTTTGCATGGTCATAGCTGTTTCCTGAGAGCTTGGCAGGTGATGACACACATTAACAAATTTCGTGAGGAGTCTCCAGAAGAATGCCATTAATTTCCATAGGCTCCGCCCCCCTGACGAGCATCACAAAAATCGACGCTCAAGTCAGAGGTGGCGAAACCCGACAGGACTATAAAGATACCAGGCGTTTCCCCCTGGAAGCTCCCTCGTGCGCTCTCCTGTTCCGACCCTGCCGCTTACCGGATACCTGTCCGCCTTTCTCCCTTCGGGAAGCGTGGCGCTTTCTCATAGCTCACGCTGTAGGTATCTCAGTTCGGTGTAGGTCGTTCGCTCCAAGCTGGGCTGTGTGCACGAACCCCCCGTTCAGCCCGACCGCTGCGCCTTATCCGGTAACTATCGTCTTGAGTCCAACCCGGTAAGACACGACTTATCGCCACTGGCAGCAGCCACTGGTAACAGGATTAGCAGAGCGAGGTATGTAGGCGGTGCTACAGAGTTCTTGAAGTGGTGGCCTAACTACGGCTACACTAGAAGAACAGTATTTGGTATCTGCGCTCTGCTGAAGCCAGTTACCTTCGGAAAAAGAGTTGGTAGCTCTTGATCCGGCAAACAAACCACCGCTGGTAGCGGTGGTTTTTTTGTTTGCAAGCAGCAGATTACGCGCAG
}
