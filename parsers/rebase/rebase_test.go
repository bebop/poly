package rebase

import (
	"fmt"
)

func ExampleReadRebase() {
	enzymeMap, _ := ReadRebase("data/link_withrefm.txt")
	fmt.Println(enzymeMap["BsaI"].MicroOrganism)
	// Output: Bacillus stearothermophilus 6-55
}

func ExampleExportRebase() {
	enzymeMap, _ := ReadRebase("data/link_withrefm.txt")
	enzymeJson, _ := ExportRebase(enzymeMap)
	fmt.Println(string(enzymeJson)[:100])
	// Output: {"AaaI":{"name":"AaaI","isoschizomers":["XmaIII","BseX3I","BsoDI","BstZI","EagI","EclXI","Eco52I","S
}
