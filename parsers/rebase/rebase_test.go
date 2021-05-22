package rebase

import (
	"fmt"
)

func ExampleReadRebase() {
	enzymeMap, _ := ReadRebase("data/rebase_test.txt")
	fmt.Println(enzymeMap["AarI"].MicroOrganism)
	// Output: Arthrobacter aurescens SS2-322
}

func ExampleExportRebase() {
	enzymeMap, _ := ReadRebase("data/rebase_test.txt")
	enzymeJson, _ := ExportRebase(enzymeMap)
	fmt.Println(string(enzymeJson)[:100])
	// Output: {"AaaI":{"name":"AaaI","isoschizomers":["XmaIII","BseX3I","BsoDI","BstZI","EagI","EclXI","Eco52I","S
}
