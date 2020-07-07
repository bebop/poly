package main

import (
	"bytes"
	"fmt"
	"testing"

	"github.com/sergi/go-diff/diffmatchpatch"
)

func TestBlake3HashRegression(t *testing.T) {
	puc19GbkBlake3Hash := "4b0616d1b3fc632e42d78521deb38b44fba95cca9fde159e01cd567fa996ceb9"
	puc19 := ReadGbk("data/puc19.gbk")
	if got := puc19.blake3Hash(); got != puc19GbkBlake3Hash {
		t.Errorf("TestBlake3HashRegression has failed. Blake3 has returned %q, want %q", got, puc19GbkBlake3Hash)
	}
}

func TestLeastRotation(t *testing.T) {
	annotatedSequence := ReadGbk("data/puc19.gbk")
	var sequenceBuffer bytes.Buffer

	sequenceBuffer.WriteString(annotatedSequence.Sequence.Sequence)
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
