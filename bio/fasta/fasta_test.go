package fasta

import (
	"bytes"
	"errors"
	"io"
	"strings"
	"testing"
)

func TestHeader(t *testing.T) {
	_, err := NewParser(nil, 256).Header()
	if err != nil {
		t.Errorf("Unexpected error while parsing header: %v", err)
	}
}

func TestParser(t *testing.T) {
	for testIndex, test := range []struct {
		content  string
		expected []Record
	}{
		{
			content:  ">humen\nGATTACA\nCATGAT", // EOF-ended Fasta is valid
			expected: []Record{{Identifier: "humen", Sequence: "GATTACACATGAT"}},
		},
		{
			content:  ">humen\nGATTACA\nCATGAT\n",
			expected: []Record{{Identifier: "humen", Sequence: "GATTACACATGAT"}},
		},
		{
			content: ">doggy or something\nGATTACA\n\nCATGAT\n\n;a fun comment\n" +
				">homunculus\nAAAA\n",
			expected: []Record{
				{Identifier: "doggy or something", Sequence: "GATTACACATGAT"},
				{Identifier: "homunculus", Sequence: "AAAA"},
			},
		},
	} {
		var fastas []Record
		parser := NewParser(strings.NewReader(test.content), 256)
		for {
			fa, err := parser.Next()
			if err != nil {
				if !errors.Is(err, io.EOF) {
					t.Errorf("Got error: %s", err)
				}
				break
			}
			fastas = append(fastas, *fa)
		}
		if len(fastas) != len(test.expected) {
			t.Errorf("case index %d: got %d fastas, expected %d", testIndex, len(fastas), len(test.expected))
			continue
		}
		for index, gotFasta := range fastas {
			expected := test.expected[index]
			if expected != gotFasta {
				t.Errorf("got!=expected: %+v != %+v", gotFasta, expected)
			}
		}
	}
}

// TestReadEmptyFasta tests that an empty fasta file is parsed correctly.
func TestReadEmptyFasta(t *testing.T) {
	var fastas []Record
	var targetError error
	emptyFasta := "testing\natagtagtagtagtagatgatgatgatgagatg\n\n\n\n\n\n\n\n\n\n\n"
	parser := NewParser(strings.NewReader(emptyFasta), 256)
	for {
		fa, err := parser.Next()
		if err != nil {
			if errors.Is(err, io.EOF) {
				err = nil // EOF not treated as parsing error.
			}
			targetError = err
			break
		}
		fastas = append(fastas, *fa)
	}
	if targetError == nil {
		t.Errorf("expected error reading empty fasta stream")
	}
	if len(fastas) != 0 {
		t.Errorf("expected 1 fastas, got %d", len(fastas))
	}
}

func TestReadEmptySequence(t *testing.T) {
	var targetError error
	emptyFasta := ">testing\natagtagtagtagtagatgatgatgatgagatg\n>testing2\n\n\n\n\n\n\n\n\n\n"
	parser := NewParser(strings.NewReader(emptyFasta), 256)
	for {
		_, err := parser.Next()
		if err != nil {
			if errors.Is(err, io.EOF) {
				err = nil // EOF not treated as parsing error.
			}
			targetError = err
			break
		}
	}
	if targetError == nil {
		t.Errorf("expected error reading empty fasta sequence stream: %s", targetError)
	}
}

func TestBufferSmall(t *testing.T) {
	var targetError error
	emptyFasta := ">test\natagtagtagtagtagatgatgatgatgagatg\n>test\n\n\n\n\n\n\n\n\n\n"
	parser := NewParser(strings.NewReader(emptyFasta), 8)
	for {
		_, err := parser.Next()
		if err != nil {
			if errors.Is(err, io.EOF) {
				err = nil // EOF not treated as parsing error.
			}
			targetError = err
			break
		}
	}
	if targetError == nil {
		t.Errorf("expected error with too small of a buffer")
	}
}

// The following functions help test for writing with a limit in order to get
// that sweet sweet test coverage.
// LimitedWriter wraps another io.Writer and returns an error if more than maxSize bytes are written.
type LimitedWriter struct {
	w       io.Writer
	written int64
	maxSize int64
}

func NewLimitedWriter(w io.Writer, maxSize int64) *LimitedWriter {
	return &LimitedWriter{
		w:       w,
		written: 0,
		maxSize: maxSize,
	}
}

func (lw *LimitedWriter) Write(p []byte) (int, error) {
	if int64(len(p)) > (lw.maxSize - lw.written) {
		return 0, errors.New("write exceeds maximum size")
	}
	n, err := lw.w.Write(p)
	lw.written += int64(n)
	return n, err
}

func TestWrite(t *testing.T) {
	// 81 polyA
	s := ">test\nAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA\n"
	parser := NewParser(strings.NewReader(s), 1024)
	fa, err := parser.Next()
	if err != nil {
		t.Errorf("Failed to read polyA: %s", err)
	}
	byteSizes := []int{0, 1, 5, 6, 7, 85, 86, 87, 89, 128}
	for _, i := range byteSizes {
		var buf bytes.Buffer
		writer := NewLimitedWriter(&buf, int64(i))
		_, _ = fa.WriteTo(writer)
	}
}
