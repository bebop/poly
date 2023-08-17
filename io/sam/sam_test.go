package sam

import (
	"errors"
	"io"
	"os"
	"testing"
)

func TestParse(t *testing.T) {
	file, err := os.Open("data/aln.sam")
	if err != nil {
		t.Errorf("Failed to open aln.sam: %s", err)
	}
	parser, header, err := NewParser(file, DefaultMaxLineSize)
	if len(header.HD) != 3 {
		t.Errorf("HD should have 3 TAG:DATA pairs")
	}
	for {
		_, err := parser.ParseNext()
		if err != nil {
			if !errors.Is(err, io.EOF) {
				t.Errorf("Got unknown error: %s", err)
			}
			break
		}
	}
}
