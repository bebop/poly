package seqhash_test

import (
	"fmt"

	"github.com/TimothyStiles/poly/bio/genbank"
	"github.com/TimothyStiles/poly/seqhash"
)

// This example shows how to seqhash a sequence.
func Example_basic() {
	sequence := "ATGC"
	sequenceType := seqhash.DNA
	circular := false
	doubleStranded := true

	sequenceSeqhash, _ := seqhash.Hash(sequence, sequenceType, circular, doubleStranded)
	fmt.Println(sequenceSeqhash)
	// Output: v1_DLD_f4028f93e08c5c23cbb8daa189b0a9802b378f1a1c919dcbcf1608a615f46350
}

func ExampleHash() {
	sequence := "ATGC"
	sequenceType := seqhash.DNA
	circular := false
	doubleStranded := true

	sequenceSeqhash, _ := seqhash.Hash(sequence, sequenceType, circular, doubleStranded)
	fmt.Println(sequenceSeqhash)
	// Output: v1_DLD_f4028f93e08c5c23cbb8daa189b0a9802b378f1a1c919dcbcf1608a615f46350
}

func ExampleRotateSequence() {
	sequence, _ := genbank.Read("../data/puc19.gbk")
	sequenceLength := len(sequence.Sequence)
	testSequence := sequence.Sequence[sequenceLength/2:] + sequence.Sequence[0:sequenceLength/2]

	fmt.Println(seqhash.RotateSequence(sequence.Sequence) == seqhash.RotateSequence(testSequence))
	// output: true
}
