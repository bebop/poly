package gff_test

import (
	"fmt"

	"github.com/TimothyStiles/poly/io/gff"
)

// This example shows how to open a gff file and search for a gene given its
// locus tag. We then display the EC number of that particular gene.
func Example_basic() {
	sequence, _ := gff.Read("../../data/ecoli-mg1655-short.gff")
	for _, feature := range sequence.Features {
		if feature.Attributes["locus_tag"] == "b0003" {
			fmt.Println(feature.Attributes["EC_number"])
		}
	}
	// Output: 2.7.1.39
}
