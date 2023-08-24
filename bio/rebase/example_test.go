package rebase_test

import (
	"fmt"

	"github.com/TimothyStiles/poly/bio/rebase"
)

// This example reads rebase into an enzymeMap and returns the AarI recognition
// sequence.
func Example_basic() {
	enzymeMap, _ := rebase.Read("data/rebase_test.txt")
	fmt.Println(enzymeMap["AarI"].RecognitionSequence)
	// Output: CACCTGC(4/8)
}

func ExampleRead() {
	enzymeMap, _ := rebase.Read("data/rebase_test.txt")
	fmt.Println(enzymeMap["AarI"].MicroOrganism)
	// Output: Arthrobacter aurescens SS2-322
}

func ExampleExport() {
	enzymeMap, _ := rebase.Read("data/rebase_test.txt")
	enzymeJSON, _ := rebase.Export(enzymeMap)
	fmt.Println(string(enzymeJSON)[:100])
	// Output: {"AaaI":{"name":"AaaI","isoschizomers":["XmaIII","BseX3I","BsoDI","BstZI","EagI","EclXI","Eco52I","S
}
