package uniprot_test

import (
	"fmt"

	"github.com/TimothyStiles/poly/bio/uniprot"
)

// This example shows how to open a uniprot data dump file and read the results
// into a list. Directly using the channel without converting to an array
// should be used for the Trembl data dump
func Example_basic() {
	entries, _, _ := uniprot.Read("data/uniprot_sprot_mini.xml.gz")

	var entry uniprot.Entry
	for singleEntry := range entries {
		entry = singleEntry
	}
	fmt.Println(entry.Accession[0])
	// Output: O55723
}
