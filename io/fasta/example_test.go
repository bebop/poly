package fasta_test

import (
	"fmt"
	"github.com/TimothyStiles/poly/io/fasta"
)

// This example shows how to open a file with the fasta parser. The sequences
// within that file can then be analyzed further with different software.
func Example_basic() {
	fastas := fasta.Read("data/base.fasta")
	fmt.Println(fastas[1].Sequence)
	// Output: ADQLTEEQIAEFKEAFSLFDKDGDGTITTKELGTVMRSLGQNPTEAELQDMINEVDADGNGTIDFPEFLTMMARKMKDTDSEEEIREAFRVFDKDGNGYISAAELRHVMTNLGEKLTDEEVDEMIREADIDGDGQVNYEEFVQMMTAK*
}
