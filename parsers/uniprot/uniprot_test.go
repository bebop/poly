package uniprot

import (
	"fmt"
)

func ExampleReadUniprot() {
	entries := make(chan Entry, 100000)
	go ReadUniprot("data/uniprot_sprot_mini.xml.gz", entries)

	entry := <-entries
	fmt.Println(entry.Accession[0])
	// Output: P0C9F0
}
