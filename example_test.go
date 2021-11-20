package poly_test

import (
	"fmt"

	"github.com/TimothyStiles/poly/io/fasta"
)

func Example_getting_started() {
	fastas := fasta.Read("data/base.fasta")
	fmt.Println(fastas[1].Sequence)
	// Output: ADQLTEEQIAEFKEAFSLFDKDGDGTITTKELGTVMRSLGQNPTEAELQDMINEVDADGNGTIDFPEFLTMMARKMKDTDSEEEIREAFRVFDKDGNGYISAAELRHVMTNLGEKLTDEEVDEMIREADIDGDGQVNYEEFVQMMTAK*
}
