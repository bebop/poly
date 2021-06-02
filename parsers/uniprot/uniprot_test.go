package uniprot

import (
	"fmt"
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

func TestReadUniprot(t *testing.T) {
	_, _, err := ReadUniprot("data/test")
	if err == nil {
		t.Errorf("Failed to fail on non-gzipped file")
	}

	_, _, err = ReadUniprot("data/FAKE")
	if err == nil {
		t.Errorf("Failed to fail on empty file")
	}
}
