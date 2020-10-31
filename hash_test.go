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
