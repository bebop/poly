package rbs_calculator

import "fmt"

func ExampleRibosomeBindingSite() {
	bindingSite, _ := RibosomeBindingSite("TATCGGGGCTTCCGTCGGCCATAAGGAGGTAAAAAATGGCGAGCTCTGAAGACGTTATCAAAGAGTTCAT")
	// if err != nil {
	// 	fmt.Printf("Failed to initialize lookup table: %s", err)
	// 	return
	// }

	fmt.Printf("%v\n", bindingSite)
	// Output: -11.9
}
