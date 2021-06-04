package rebase

import (
	"fmt"
	"testing"
)

func ExampleReadRebase() {
	enzymeMap, _ := ReadRebase("data/rebase_test.txt")
	fmt.Println(enzymeMap["AarI"].MicroOrganism)
	// Output: Arthrobacter aurescens SS2-322
}

func ExampleExportRebase() {
	enzymeMap, _ := ReadRebase("data/rebase_test.txt")
	enzymeJson := ExportRebase(enzymeMap)
	fmt.Println(string(enzymeJson)[:100])
	// Output: {"AaaI":{"name":"AaaI","isoschizomers":["XmaIII","BseX3I","BsoDI","BstZI","EagI","EclXI","Eco52I","S
}

func TestReadRebase(t *testing.T) {
	_, err := ReadRebase("data/FAKE.txt")
	if err == nil {
		t.Errorf("Failed to error on fake file")
	}
}
