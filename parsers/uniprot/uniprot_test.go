package uniprot

import (
	"fmt"
)

func ExampleReadUniprot() {
	entries, _, _ := ReadUniprot("data/uniprot_sprot_mini.xml.gz")

	entry := <-entries
	fmt.Println(entry.Accession[0])
	// Output: P0C9F0
}
