package seqhash

import (
	"bytes"
	"fmt"
	"os"
	"testing"

	"github.com/TimothyStiles/poly/bio"
	"github.com/sergi/go-diff/diffmatchpatch"
)

func TestHash(t *testing.T) {
	// Test TNA as sequenceType
	_, err := Hash("ATGGGCTAA", "TNA", true, true)
	if err == nil {
		t.Errorf("TestHash() has failed. TNA is not a valid sequenceType.")
	}
	// Test X in DNA or RNA
	_, err = Hash("XTGGCCTAA", "DNA", true, true)
	if err == nil {
		t.Errorf("TestSeqhashSequenceString() has failed. X is not a valid DNA or RNA sequence character.")
	}
	// Test X in PROTEIN
	_, err = Hash("MGCJ*", "PROTEIN", false, false)
	if err == nil {
		t.Errorf("TestSeqhashSequenceProteinString() has failed. J is not a valid PROTEIN sequence character.")
		fmt.Println(err)
	}
	// Test double stranded Protein
	_, err = Hash("MGCS*", "PROTEIN", false, true)
	if err == nil {
		t.Errorf("TestSeqhashProteinDoubleStranded() has failed. Proteins cannot be double stranded.")
	}

	// Test circular double stranded hashing
	seqhash, _ := Hash("TTAGCCCAT", "DNA", true, true)
	if seqhash != "v1_DCD_a376845b679740014f3eb501429b45e592ecc32a6ba8ba922cbe99217f6e9287" {
		t.Errorf("Circular double stranded hashing failed. Expected v1_DCD_a376845b679740014f3eb501429b45e592ecc32a6ba8ba922cbe99217f6e9287, got: " + seqhash)
	}
	// Test circular single stranded hashing
	seqhash, _ = Hash("TTAGCCCAT", "DNA", true, false)
	if seqhash != "v1_DCS_ef79b6e62394e22a176942dfc6a5e62eeef7b5281ffcb2686ecde208ec836ba4" {
		t.Errorf("Circular single stranded hashing failed. Expected v1_DCS_ef79b6e62394e22a176942dfc6a5e62eeef7b5281ffcb2686ecde208ec836ba4, got: " + seqhash)
	}
	// Test linear double stranded hashing
	seqhash, _ = Hash("TTAGCCCAT", "DNA", false, true)
	if seqhash != "v1_DLD_c2c9fc44df72035082a152e94b04492182331bc3be2f62729d203e072211bdbf" {
		t.Errorf("Linear double stranded hashing failed. Expected v1_DLD_c2c9fc44df72035082a152e94b04492182331bc3be2f62729d203e072211bdbf, got: " + seqhash)
	}
	// Test linear single stranded hashing
	seqhash, _ = Hash("TTAGCCCAT", "DNA", false, false)
	if seqhash != "v1_DLS_063ea37d1154351639f9a48546bdae62fd8a3c18f3d3d3061060c9a55352d967" {
		t.Errorf("Linear single stranded hashing failed. Expected v1_DLS_063ea37d1154351639f9a48546bdae62fd8a3c18f3d3d3061060c9a55352d967, got: " + seqhash)
	}

	// Test RNA Seqhash
	seqhash, _ = Hash("TTAGCCCAT", "RNA", false, false)
	if seqhash != "v1_RLS_063ea37d1154351639f9a48546bdae62fd8a3c18f3d3d3061060c9a55352d967" {
		t.Errorf("Linear single stranded hashing failed. Expected v1_RLS_063ea37d1154351639f9a48546bdae62fd8a3c18f3d3d3061060c9a55352d967, got: " + seqhash)
	}
	// Test Protein Seqhash
	seqhash, _ = Hash("MGC*", "PROTEIN", false, false)
	if seqhash != "v1_PLS_922ec11f5227ce77a42f07f565a7a1a479772b5cf3f1f6e93afc5ecbc0fd5955" {
		t.Errorf("Linear single stranded hashing failed. Expected v1_PLS_922ec11f5227ce77a42f07f565a7a1a479772b5cf3f1f6e93afc5ecbc0fd5955, got: " + seqhash)
	}
}

func TestLeastRotation(t *testing.T) {
	file, _ := os.Open("../data/puc19.gbk")
	defer file.Close()
	parser, _ := bio.NewGenbankParser(file)
	sequence, _ := parser.Next()
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
