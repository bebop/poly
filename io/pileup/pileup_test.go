package pileup

import (
	"errors"
	"io"
	"os"
	"testing"
)

func TestParse(t *testing.T) {
	file, err := os.Open("data/test3.pileup")
	if err != nil {
		t.Errorf("Failed to open test.pileup: %s", err)
	}
	parser := NewParser(file)
	var pileupReads []Pileup
	for {
		pileupRead, err := parser.ParseNext()
		if err != nil {
			if !errors.Is(err, io.EOF) {
				t.Errorf("Got unknown error: %s", err)
			}
			break
		}
		pileupReads = append(pileupReads, pileupRead)
	}
	//if pileupReads[2].ReadCount != 1191 {
	//	t.Errorf("Expected 1191 read counts. Got: %d", pileupReads[2].ReadCount)
	//}
}
