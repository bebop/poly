package uniprot

import (
	"compress/gzip"
	"fmt"
	"os"
	"testing"
)

func ExampleReadUniprot() {
	entries, _, _ := ReadUniprot("data/uniprot_sprot_mini.xml.gz")

	var entry Entry
	for singleEntry := range entries {
		entry = singleEntry
	}
	fmt.Println(entry.Accession[0])
	// Output: O55723
}

func ExampleParseUniprot() {
	xmlFile, _ := os.Open("data/uniprot_sprot_mini.xml.gz")
	unzippedBytes, _ := gzip.NewReader(xmlFile)

	entries := make(chan Entry, 100) // if you don't have a buffered channel, nothing will be read in loops on the channel.
	decoderErrors := make(chan error, 100)
	go ParseUniprot(unzippedBytes, entries, decoderErrors)

	var entry Entry
	for singleEntry := range entries {
		entry = singleEntry
	}
	fmt.Println(entry.Accession[0])
	// Output: O55723
}

func TestReadUniprot(t *testing.T) {
	_, _, err := ReadUniprot("data/test")
	if err == nil {
		t.Errorf("Failed to fail on non-gzipped file")
	}

	_, _, err = ReadUniprot("data/FAKE")
	if err == nil {
		t.Errorf("Failed to fail on empty file")
	}

	_, errors, err := ReadUniprot("data/uniprot_sprot_mini.xml.gz")
	if err != nil {
		t.Errorf("Failed on real file with error: %v", err)
	}

	for err := range errors {
		if err != nil {
			t.Errorf("Failed during parsing with error: %v", err)
		}
	}
}
