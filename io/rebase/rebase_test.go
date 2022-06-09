package rebase

import (
	"fmt"
	"testing"
)

func ExampleRead() {
	enzymeMap, _ := Read("data/rebase_test.txt")
	fmt.Println(enzymeMap["AarI"].MicroOrganism)
	// Output: Arthrobacter aurescens SS2-322
}

func ExampleExport() {
	enzymeMap, _ := Read("data/rebase_test.txt")
	enzymeJSON := Export(enzymeMap)
	fmt.Println(string(enzymeJSON)[:100])
	// Output: {"AaaI":{"name":"AaaI","isoschizomers":["XmaIII","BseX3I","BsoDI","BstZI","EagI","EclXI","Eco52I","S
}

func TestRead(t *testing.T) {
	_, err := Read("data/FAKE.txt")
	if err == nil {
		t.Errorf("Failed to error on fake file")
	}
}
