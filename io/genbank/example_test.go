package genbank_test

import (
	"fmt"

	"github.com/TimothyStiles/poly/io/genbank"
)

// This example shows how to open a genbank file and search for a gene given
// its name. After finding it, notes about the particular gene are read.
func Example_basic() {
	sequence, _ := genbank.Read("../../data/puc19.gbk")
	for _, feature := range sequence.Features {
		if feature.Attributes["gene"] == "bla" {
			fmt.Println(feature.Attributes["note"])
		}
	}
	// Output: confers resistance to ampicillin, carbenicillin, andrelated antibiotics
}
