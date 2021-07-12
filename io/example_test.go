package io

import (
	"github.com/TimothyStiles/poly/io/fasta"
	"github.com/TimothyStiles/poly/io/genbank"
	"github.com/TimothyStiles/poly/io/gff"
	"github.com/TimothyStiles/poly/io/polyjson"
)

// This is where the integration tests that make effed up cyclic dependencies go.

func Example() {

	// Poly can take in basic gff, gbk, fasta, and JSON.
	// We call the json package "pson" (poly JSON) to prevent namespace collision with Go's standard json package.

	gffInput := gff.Read("../data/ecoli-mg1655-short.gff")
	gbkInput := genbank.Read("../data/puc19.gbk")
	fastaInput := fasta.Read("fasta/data/base.fasta")
	jsonInput := polyjson.Read("../data/puc19static.json")

	// Poly can also output these file formats though I wouldn't try doing gbk<->gff or anything like that unless it's JSON.

	gff.Write(gffInput, "test.gff")
	genbank.Write(gbkInput, "test.gbk")
	fasta.Write(fastaInput, "test.fasta")
	polyjson.Write(jsonInput, "test.json")

	// Poly can also output digested gff and gbk as JSON.

	polyjson.Write(gffInput, "gff_test.json")
	polyjson.Write(gbkInput, "gbk_test.json")

}
