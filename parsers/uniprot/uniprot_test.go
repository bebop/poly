package uniprot

import (
	"fmt"
)

func ExampleReadUniprot() {
	entries := make(chan Entry, 100000)
	go ReadUniprot("data/uniprot_sprot_mini.xml.gz", entries)
	var entry Entry
	for e := range entries {
		entry = e
	}

	fmt.Println(entry.Accession[0])
	// Output: O55723
}
