package main

import (
	"bytes"
	"fmt"
	"testing"

	"github.com/sergi/go-diff/diffmatchpatch"
)

func TestBlake3HashRegression(t *testing.T) {
	puc19GbkBlake3Hash := "4031e1971acc8ff1bf0aa4ed623bc58beefc15e043075866a0854d592d80b28b"
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
