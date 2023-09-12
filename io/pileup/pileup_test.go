package pileup

import (
	"bytes"
	"errors"
	"fmt"
	"io"
	"os"
	"strings"
	"testing"
)

func TestParse(t *testing.T) {
	file, err := os.Open("data/test.pileup")
	if err != nil {
		t.Errorf("Failed to open test.pileup: %s", err)
	}
	defer file.Close()
	const maxLineSize = 2 * 32 * 1024
	parser := NewParser(file, maxLineSize)
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
	if pileupReads[2].ReadCount != 1401 {
		t.Errorf("Expected 1401 read counts. Got: %d", pileupReads[2].ReadCount)
	}

	// Test not enough fields
	file2, err := os.Open("data/test_not_enough_fields.pileup")
	if err != nil {
		t.Errorf("Failed to open test_not_enough_fields.pileup: %s", err)
	}
	defer file2.Close()
	parser = NewParser(file2, maxLineSize)
	for {
		_, err = parser.ParseNext()
		if err != nil {
			if !strings.Contains(fmt.Sprint(err), "values, expected 6.") {
				t.Errorf("Got unknown error: %s", err)
			}
			break
		}
	}

	// Test Position Int
	file3, err := os.Open("data/test_position_non_int.pileup")
	if err != nil {
		t.Errorf("Failed to open test_position_non_int.pileup: %s", err)
	}
	defer file3.Close()
	parser = NewParser(file3, maxLineSize)
	for {
		_, err = parser.ParseNext()
		if err != nil {
			if !strings.Contains(fmt.Sprint(err), "Error on line") {
				t.Errorf("Got unknown error: %s", err)
			}
			break
		}
	}

	// Test ReadCount Int
	file4, err := os.Open("data/test_readcount_non_int.pileup")
	if err != nil {
		t.Errorf("Failed to open test_readcount_non_int.pileup: %s", err)
	}
	defer file4.Close()
	parser = NewParser(file4, maxLineSize)
	for {
		_, err = parser.ParseNext()
		if err != nil {
			if !strings.Contains(fmt.Sprint(err), "Error on line") {
				t.Errorf("Got unknown error: %s", err)
			}
			break
		}
	}

	// Test +,- Rune
	file5, err := os.Open("data/test_read_off.pileup")
	if err != nil {
		t.Errorf("Failed to open test_read_off.pileup: %s", err)
	}
	defer file5.Close()
	parser = NewParser(file5, maxLineSize)
	for {
		_, err = parser.ParseNext()
		if err != nil {
			if !strings.Contains(fmt.Sprint(err), "Rune within +,- not found on line") {
				t.Errorf("Got unknown error: %s", err)
			}
			break
		}
	}

	// Test Unknown Rune
	file6, err := os.Open("data/test_unknown_rune.pileup")
	if err != nil {
		t.Errorf("Failed to open test_unknown_rune.pileup: %s", err)
	}
	defer file6.Close()
	parser = NewParser(file6, maxLineSize)
	for {
		_, err = parser.ParseNext()
		if err != nil {
			if !strings.Contains(fmt.Sprint(err), "Rune not found on line") {
				t.Errorf("Got unknown error: %s", err)
			}
			break
		}
	}

	// Test ParseN
	file7, err := os.Open("data/test.pileup")
	if err != nil {
		t.Errorf("Failed to open test.pileup: %s", err)
	}
	defer file7.Close()
	parser = NewParser(file7, maxLineSize)
	_, err = parser.ParseN(2)
	if err != nil {
		t.Errorf("Got unknown error: %s", err)
	}
}

func TestBuild(t *testing.T) {
	file, err := os.Open("data/test.pileup")
	if err != nil {
		t.Errorf("Failed to open test.pileup: %s", err)
	}
	const maxLineSize = 2 * 32 * 1024
	parser := NewParser(file, maxLineSize)
	pileups, err := parser.ParseAll() // Parse all
	if err != nil {
		t.Errorf("Failed to parse pileup: %s", err)
	}
	parser.Reset(file)
	var b bytes.Buffer
	if err := WritePileups(pileups, &b); err != nil {
		t.Errorf("Failed to write pileups: %v", err)
	}
	newBody := b.Bytes()
	file.Close()

	// Open file again, read all
	file, err = os.Open("data/test.pileup")
	if err != nil {
		t.Errorf("Failed to open test.pileup: %s", err)
	}
	body, err := io.ReadAll(file) // Read them bytes
	if err != nil {
		t.Errorf("Failed to readall: %s", err)
	}
	if string(body) != string(newBody) {
		t.Errorf("input bytes do not equal output bytes of build")
	}
}

func TestRead(t *testing.T) {
	_, err := Read("data/DOESNT_EXIST.pileup")
	if err == nil {
		t.Errorf("Should have failed on empty file")
	}
	_, err = Read("data/test.pileup")
	if err != nil {
		t.Errorf("Should have succeeded in reading pileup file")
	}
}

func TestWrite(t *testing.T) {
	pileups, _ := Read("data/test.pileup")
	err := Write(pileups, "tmp.pileup")
	if err != nil {
		t.Errorf("Failed to write temporary pileup")
	}
	err = os.Remove("tmp.pileup")
	if err != nil {
		t.Errorf("Failed to delete temporary pileup")
	}
}
