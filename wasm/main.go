/*
Package main compiles poly to have convenient wasm bindings.
*/
package main

import (
	sq "github.com/TimothyStiles/poly/seqhash"
)

//go:wasm-module poly
//export seqhash
func seqhash(sequence string, sequenceType string, circular bool, doubleStranded bool) (string, error) {
	sequenceTypeEnum := sq.SequenceType(sequenceType)
	return sq.Hash(sequence, sequenceTypeEnum, circular, doubleStranded)
}

// main is required for the `wasi` target, even if it isn't used.
func main() {}
