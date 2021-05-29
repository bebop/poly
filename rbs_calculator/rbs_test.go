package rbs_calculator

import "fmt"

func ExampleRibosomeBindingSite() {
	bindingSite, _ := RibosomeBindingSite("TATCGGGGCTTCCGTCGGCCATAAGGAGGTAAAAAATGGCGAGCTCTGAAGACGTTATCAAAGAGTTCAT")
	// if err != nil {
	// 	fmt.Printf("Failed to initialize lookup table: %s", err)
	// 	return
	// }

	bindingSiteSequence := bindingSite.mrna[bindingSite.fivePrimeIdx:bindingSite.threePrimeIdx]
	fmt.Printf("binding site sequence: %v (%v,%v)\n", bindingSiteSequence, bindingSite.fivePrimeIdx, bindingSite.threePrimeIdx)
	fmt.Printf("mfe: %v\n", bindingSite.minimumFreeEnergy)
	// Output:
	// binding site sequence: GGUAAAAAAUG (27,38)
	// mfe: -18.2
}
