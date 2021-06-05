package uniprot

import (
	"compress/gzip"
	"fmt"
	"os"
	"testing"
)

func ExampleRead() {
	entries, _, _ := Read("data/uniprot_sprot_mini.xml.gz")

	var entry Entry
	for singleEntry := range entries {
		entry = singleEntry
	}
	fmt.Println(entry.Accession[0])
	// Output: O55723
}

func ExampleParse() {
	xmlFile, _ := os.Open("data/uniprot_sprot_mini.xml.gz")
	unzippedBytes, _ := gzip.NewReader(xmlFile)

	entries := make(chan Entry, 100) // if you don't have a buffered channel, nothing will be read in loops on the channel.
	decoderErrors := make(chan error, 100)
	go Parse(unzippedBytes, entries, decoderErrors)

	var entry Entry
	for singleEntry := range entries {
		entry = singleEntry
	}
	fmt.Println(entry.Accession[0])
	// Output: O55723
}

func TestRead(t *testing.T) {
	_, _, err := Read("data/test")
	if err == nil {
		t.Errorf("Failed to fail on non-gzipped file")
	}

	_, _, err = Read("data/FAKE")
	if err == nil {
		t.Errorf("Failed to fail on empty file")
	}
}
