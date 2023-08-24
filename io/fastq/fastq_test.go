package fastq

import (
	"errors"
	"fmt"
	"os"
	"strings"
	"testing"
)

func TestParseNLow(t *testing.T) {
	file, err := os.Open("data/nanosavseq.fastq")
	if err != nil {
		t.Errorf("Failed to read nanosavseq.fastq. Got error: %s", err)
	}
	const maxLineSize = 2 * 32 * 1024
	parser := NewParser(file, maxLineSize)
	_, err = parser.ParseN(0)
	if err != nil {
		t.Errorf("Failed to parse 0 fastqs. Got error: %s", err)
	}
}

func TestParseSmallLine(t *testing.T) {
	file, _ := os.Open("data/nanosavseq.fastq")
	parser := NewParser(file, 10)
	_, err := parser.ParseAll()
	if err == nil {
		t.Errorf("Should have encountered a maxLine error")
	}
	parser.Reset(file)
}

func TestRead(t *testing.T) {
	_, err := Read("data/doesntexist.fastq")
	if err == nil {
		t.Errorf("Should have failed to read non-existent file")
	}
	_, err = ReadGz("data/doesntexist.fastq.gz")
	if err == nil {
		t.Errorf("Should have failed to read non-existent gz file")
	}
	_, err = ReadGz("data/nanosavseq.fastq")
	if err == nil {
		t.Errorf("Should have failed to read a file that is not gz'ed")
	}
}

func testException(t *testing.T, filePath string, errorString string) {
	file, err := os.Open(filePath)
	if err != nil {
		t.Errorf("Failed to read %s. Got error: %s", filePath, err)
	}
	const maxLineSize = 2 * 32 * 1024
	parser := NewParser(file, maxLineSize)
	_, err = parser.ParseAll()
	if err == nil {
		t.Errorf("%s parser should have gotten error: %s", filePath, errorString)
	}
}

func TestParseExceptions(t *testing.T) {
	testException(t, "data/nanosavseq_noseq.fastq", "no seq")
	testException(t, "data/nanosavseq_noquality.fastq", "no quality")
	testException(t, "data/nanosavseq_noidentifier.fastq", "no identifier")
	testException(t, "data/nanosavseq_emptyseq.fastq", "empty seq")
	testException(t, "data/nanosavseq_noplus.fastq", "no plus EOF")
	testException(t, "data/nanosavseq_noquality2.fastq", "no quality EOF")
}

func TestGzpOpen(t *testing.T) {
	fastqs, err := ReadGz("data/pOpen_v3.fastq.gz")
	if err != nil {
		t.Errorf("Failed to open pOpen_v3: %s", err)
	}
	for _, read := range fastqs {
		id := read.Identifier
		sequence := read.Sequence
		quality := read.Quality
		for _, char := range sequence {
			if !strings.Contains("ATGCYRSWKMBDHVNZ", string(char)) {
				fmt.Println(id, sequence, quality)
				t.Fatalf("%s", errors.New("Only letters ATUGCYRSWKMBDHVNZ are allowed for DNA/RNA. Got letter: "+string(char)))
			}
		}
	}
}
