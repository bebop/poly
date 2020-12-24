package poly

import (
	"bytes"
	"fmt"
	"testing"

	"github.com/sergi/go-diff/diffmatchpatch"
	"lukechampine.com/blake3"
)

func ExampleHash() {
	puc19 := ReadGbk("data/puc19.gbk")
	fmt.Println(puc19.Hash(blake3.New(32, nil))) // passing new hash.Hash struct to Hasher

	// output: 4b0616d1b3fc632e42d78521deb38b44fba95cca9fde159e01cd567fa996ceb9
}

func TestHashRegression(t *testing.T) {
	puc19GbkBlake3Hash := "4b0616d1b3fc632e42d78521deb38b44fba95cca9fde159e01cd567fa996ceb9"
	puc19 := ReadGbk("data/puc19.gbk")
	if got := puc19.Hash(blake3.New(32, nil)); got != puc19GbkBlake3Hash {
		t.Errorf("TestHashRegression has failed. Blake3 sequence hash has returned %q, want %q", got, puc19GbkBlake3Hash)
	}

	// testing feature hash method
	ampRHash := "e2bc8192c23919ce4e2b2b03f7af644ab8b714865a429408c57d30fa57557a85"
	for _, feature := range puc19.Features {
		if feature.Attributes["label"] == "AmpR" {
			hash := feature.Hash(blake3.New(32, nil))
			if hash != ampRHash {
				t.Errorf("TestHashRegression has failed. Blake3 feature sequence hash has returned %q, want %q", hash, ampRHash)
			}
		}
	}
}

func ExampleRotateSequence() {
	sequence := ReadGbk("data/puc19.gbk")
	sequenceLength := len(sequence.Sequence)
	testSequence := sequence.Sequence[sequenceLength/2:] + sequence.Sequence[0:sequenceLength/2]

	fmt.Println(RotateSequence(sequence.Sequence) == RotateSequence(testSequence))
	// output: true
}

func TestLeastRotation(t *testing.T) {
	sequence := ReadGbk("data/puc19.gbk")
	var sequenceBuffer bytes.Buffer

	sequenceBuffer.WriteString(sequence.Sequence)
	bufferLength := sequenceBuffer.Len()

	var rotatedSequence string
	for elementIndex := 0; elementIndex < bufferLength; elementIndex++ {
		bufferElement, _, _ := sequenceBuffer.ReadRune()
		sequenceBuffer.WriteRune(bufferElement)
		if elementIndex == 0 {
			rotatedSequence = RotateSequence(sequenceBuffer.String())
		} else {
			newRotatedSequence := RotateSequence(sequenceBuffer.String())
			if rotatedSequence != newRotatedSequence {
				dmp := diffmatchpatch.New()
				diffs := dmp.DiffMain(rotatedSequence, newRotatedSequence, false)
				t.Errorf("TestLeastRotation() has failed. rotationSequence mutated.")
				fmt.Println(dmp.DiffPrettyText(diffs))
			}
		}
	}

}

/******************************************************************************

Seqhash related tests begin here.

******************************************************************************/

func ExampleSeqhash() {
	sequence := ReadGbk("data/puc19.gbk")

	seqhash, _ := Seqhash(sequence.Sequence, "DNA", true, true)
	fmt.Println(seqhash)
	// output: v1_DCD_4b0616d1b3fc632e42d78521deb38b44fba95cca9fde159e01cd567fa996ceb9
}

func ExampleSequenceSeqhash() {
	sequence := ReadGbk("data/puc19.gbk")

	// SequenceSeqhash assumes doubleStranded sequence and defaults to linear
	// if sequence.Meta.Locus.Circular is not set
	seqhash, _ := sequence.SequenceSeqhash()
	fmt.Println(seqhash)
	// output: v1_DCD_4b0616d1b3fc632e42d78521deb38b44fba95cca9fde159e01cd567fa996ceb9
}

func TestSeqhash(t *testing.T) {
	// Test TNA as sequenceType
	_, err := Seqhash("ATGGGCTAA", "TNA", true, true)
	if err == nil {
		t.Errorf("TestSeqhash() has failed. TNA is not a valid sequenceType.")
	}
	// Test X in DNA or RNA
	_, err = Seqhash("XTGGCCTAA", "DNA", true, true)
	if err == nil {
		t.Errorf("TestSeqhashSequenceString() has failed. X is not a valid DNA or RNA sequence character.")
	}
	// Test X in PROTEIN
	_, err = Seqhash("MGCB*", "PROTEIN", false, false)
	if err == nil {
		t.Errorf("TestSeqhashSequenceProteinString() has failed. B is not a valid PROTEIN sequence character.")
		fmt.Println(err)
	}
	// Test double stranded Protein
	_, err = Seqhash("MGCS*", "PROTEIN", false, true)
	if err == nil {
		t.Errorf("TestSeqhashProteinDoubleStranded() has failed. Proteins cannot be double stranded.")
	}

	// Test circular double stranded hashing
	seqhash, _ := Seqhash("TTAGCCCAT", "DNA", true, true)
	if seqhash != "v1_DCD_a376845b679740014f3eb501429b45e592ecc32a6ba8ba922cbe99217f6e9287" {
		t.Errorf("Circular double stranded hashing failed. Expected v1_DCD_a376845b679740014f3eb501429b45e592ecc32a6ba8ba922cbe99217f6e9287, got: " + seqhash)
	}
	// Test circular single stranded hashing
	seqhash, _ = Seqhash("TTAGCCCAT", "DNA", true, false)
	if seqhash != "v1_DCS_ef79b6e62394e22a176942dfc6a5e62eeef7b5281ffcb2686ecde208ec836ba4" {
		t.Errorf("Circular single stranded hashing failed. Expected v1_DCS_ef79b6e62394e22a176942dfc6a5e62eeef7b5281ffcb2686ecde208ec836ba4, got: " + seqhash)
	}
	// Test linear double stranded hashing
	seqhash, _ = Seqhash("TTAGCCCAT", "DNA", false, true)
	if seqhash != "v1_DLD_c2c9fc44df72035082a152e94b04492182331bc3be2f62729d203e072211bdbf" {
		t.Errorf("Linear double stranded hashing failed. Expected v1_DLD_c2c9fc44df72035082a152e94b04492182331bc3be2f62729d203e072211bdbf, got: " + seqhash)
	}
	// Test linear single stranded hashing
	seqhash, _ = Seqhash("TTAGCCCAT", "DNA", false, false)
	if seqhash != "v1_DLS_063ea37d1154351639f9a48546bdae62fd8a3c18f3d3d3061060c9a55352d967" {
		t.Errorf("Linear single stranded hashing failed. Expected v1_DLS_063ea37d1154351639f9a48546bdae62fd8a3c18f3d3d3061060c9a55352d967, got: " + seqhash)
	}

	// Test RNA Seqhash
	seqhash, _ = Seqhash("TTAGCCCAT", "RNA", false, false)
	if seqhash != "v1_RLS_063ea37d1154351639f9a48546bdae62fd8a3c18f3d3d3061060c9a55352d967" {
		t.Errorf("Linear single stranded hashing failed. Expected v1_RLS_063ea37d1154351639f9a48546bdae62fd8a3c18f3d3d3061060c9a55352d967, got: " + seqhash)
	}
	// Test Protein Seqhash
	seqhash, _ = Seqhash("MGC*", "PROTEIN", false, false)
	if seqhash != "v1_PLS_922ec11f5227ce77a42f07f565a7a1a479772b5cf3f1f6e93afc5ecbc0fd5955" {
		t.Errorf("Linear single stranded hashing failed. Expected v1_PLS_922ec11f5227ce77a42f07f565a7a1a479772b5cf3f1f6e93afc5ecbc0fd5955, got: " + seqhash)
	}

}

func TestSeqhashMethod(t *testing.T) {
	// Test no MoleculeType
	sequence := ReadGbk("data/puc19.gbk")
	sequence.Meta.Locus.MoleculeType = ""
	_, err := sequence.SequenceSeqhash()
	if err == nil {
		t.Errorf("Seqhash method should fail with no MoleculeType")
	}
	// Test RNA sequence feature
	sequence.Meta.Locus.MoleculeType = "RNA"
	seqhash, _ := sequence.SequenceSeqhash()
	if seqhash != "v1_RCD_4b0616d1b3fc632e42d78521deb38b44fba95cca9fde159e01cd567fa996ceb9" {
		t.Errorf("RNA hashing failed. Expected v1_RCD_4b0616d1b3fc632e42d78521deb38b44fba95cca9fde159e01cd567fa996ceb9, got: " + seqhash)
	}
	// Test bad MoleculeType
	sequence.Meta.Locus.MoleculeType = "TNA"
	_, err = sequence.SequenceSeqhash()
	if err == nil {
		t.Errorf("Seqhash method should fail with TNA (or other) MoleculeType")
	}
	// Test bad sequence
	sequence.Meta.Locus.MoleculeType = "DNA"
	sequence.Sequence = "gagatacctacagcgtgagctatgagaaagcgccacgcttcccgaagggagaaaggcggacaggtatccggtaagcggcagggtcggaacaggagagcgcacgagggagcttccagggggaaacgcctggtatctttatagtcctgtcgggtttcgccacctctgacttgagcgtcgatttttgtgatgctcgtcaggggggcggagcctatggaaaaacgccagcaacgcggcctttttacggttcctggccttttgctggccttttgctcacatgttctttcctgcgttatcccctgattctgtggataaccgtattaccgcctttgagtgagctgataccgctcgccgcagccgaacgaccgagcgcagcgagtcagtgagcgaggaagcggaagagcgcccaatacgcaaaccgcctctccccgcgcgttggccgattcattaatgcagctggcacgacaggtttcccgactggaaagcgggcagtgagcgcaacgcaattaatgtgagttagctcactcattaggcaccccaggctttacactttatgcttccggctcgtatgttgtgtggaattgtgagcggataacaatttcacacaggaaacagctatgaccatgattacgccaagcttgcatgcctgcaggtcgactctagaggatccccgggtaccgagctcgaattcactggccgtcgttttacaacgtcgtgactgggaaaaccctggcgttacccaacttaatcgccttgcagcacatccccctttcgccagctggcgtaatagcgaagaggcccgcaccgatcgcccttcccaacagttgcgcagcctgaatggcgaatggcgcctgatgcggtattttctccttacgcatctgtgcggtatttcacaccgcatatggtgcactctcagtacaatctgctctgatgccgcatagttaagccagccccgacacccgccaacacccgctgacgcgccctgacgggcttgtctgctcccggcatccgcttacagacaagctgtgaccgtctccgggagctgcatgtgtcagaggttttcaccgtcatcaccgaaacgcgcgagacgaaagggcctcgtgatacgcctatttttataggttaatgtcatgataataatggtttcttagacgtcaggtggcacttttcggggaaatgtgcgcggaacccctatttgtttatttttctaaatacattcaaatatgtatccgctcatgagacaataaccctgataaatgcttcaataatattgaaaaaggaagagtatgagtattcaacatttccgtgtcgcccttattcccttttttgcggcattttgccttcctgtttttgctcacccagaaacgctggtgaaagtaaaagatgctgaagatcagttgggtgcacgagtgggttacatcgaactggatctcaacagcggtaagatccttgagagttttcgccccgaagaacgttttccaatgatgagcacttttaaagttctgctatgtggcgcggtattatcccgtattgacgccgggcaagagcaactcggtcgccgcatacactattctcagaatgacttggttgagtactcaccagtcacagaaaagcatcttacggatggcatgacagtaagagaattatgcagtgctgccataaccatgagtgataacactgcggccaacttacttctgacaacgatcggaggaccgaaggagctaaccgcttttttgcacaacatgggggatcatgtaactcgccttgatcgttgggaaccggagctgaatgaagccataccaaacgacgagcgtgacaccacgatgcctgtagcaatggcaacaacgttgcgcaaactattaactggcgaactacttactctagcttcccggcaacaattaatagactggatggaggcggataaagttgcaggaccacttctgcgctcggcccttccggctggctggtttattgctgataaatctggagccggtgagcgtgggtctcgcggtatcattgcagcactggggccagatggtaagccctcccgtatcgtagttatctacacgacggggagtcaggcaactatggatgaacgaaatagacagatcgctgagataggtgcctcactgattaagcattggtaactgtcagaccaagtttactcatatatactttagattgatttaaaacttcatttttaatttaaaaggatctaggtgaagatcctttttgataatctcatgaccaaaatcccttaacgtgagttttcgttccactgagcgtcagaccccgtagaaaagatcaaaggatcttcttgagatcctttttttctgcgcgtaatctgctgcttgcaaacaaaaaaaccaccgctaccagcggtggtttgtttgccggatcaagagctaccaactctttttccgaaggtaactggcttcagcagagcgcagataccaaatactgttcttctagtgtagccgtagttaggccaccacttcaagaactctgtagcaccgcctacatacctcgctctgctaatcctgttaccagtggctgctgccagtggcgataagtcgtgtcttaccgggttggactcaagacgatagttaccggataaggcgcagcggtcgggctgaacggggggttcgtgcacacagcccagcttggagcgaacgacctacaccgaactx"
	_, err = sequence.SequenceSeqhash()
	if err == nil {
		t.Errorf("Seqhash method should fail with X nucleotide")
	}
}
