package seqhash_test

import (
	"fmt"
	"github.com/TimothyStiles/poly/seqhash"
)

// This example shows how to seqhash a sequence.
func Example_basic() {
	sequence := "ATGC"
	sequenceType := "DNA"
	circular := false
	doubleStranded := true

	sequenceSeqhash, _ := seqhash.Hash(sequence, sequenceType, circular, doubleStranded)
	fmt.Println(sequenceSeqhash)
	// Output: v1_DLD_f4028f93e08c5c23cbb8daa189b0a9802b378f1a1c919dcbcf1608a615f46350
}
