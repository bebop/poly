package io_test

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

	gffInput, _ := gff.Read("../data/ecoli-mg1655-short.gff")
	gbkInput, _ := genbank.Read("../data/puc19.gbk")
	fastaInput, _ := fasta.Read("fasta/data/base.fasta")
	jsonInput, _ := polyjson.Read("../data/cat.json")

	// Poly can also output these file formats. Every file format has a corresponding Write function.
	_ = gff.Write(gffInput, "test.gff")
	_ = genbank.Write(gbkInput, "test.gbk")
	_ = fasta.Write(fastaInput, "test.fasta")
	_ = polyjson.Write(jsonInput, "test.json")

	// Extra tips:

	// 1. All of these file formats can be read and written in JSON format using their native schemas.
	// 2. If you want to convert from one format to another (e.g. genbank to polyjson), you can easily do so with a for-loop and some field mapping.
	// 3. Every file format is unique but they all share a common interface so you can use them with almost every native function in Poly.
}
