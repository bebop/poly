package rebase_test

import (
	"fmt"
	"github.com/TimothyStiles/poly/io/rebase"
)

// This example reads rebase into an enzymeMap and returns the AarI recognition
// sequence.
func Example_basic() {
	enzymeMap, _ := rebase.Read("data/rebase_test.txt")
	fmt.Println(enzymeMap["AarI"].RecognitionSequence)
	// Output: CACCTGC(4/8)
}
