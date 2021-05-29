package rbs_calculator

import (
	"fmt"
)

func ExampleRibosomeBindingSite() {
	bindingSite, err := RibosomeBindingSite(EColiRNA, "TATCGGGGCTTCCGTCGGCCATAAGGAGGTAAAAAATGGCGAGCTCTGAAGACGTTATCAAAGAGTTCAT")
	if err != nil {
		fmt.Printf("Failed to initialize lookup table: %s", err)
		return
	}

	bindingSiteSequence := bindingSite.MRNA[bindingSite.FivePrimeIdx:bindingSite.ThreePrimeIdx]
	fmt.Printf("binding site sequence: %v (%v,%v)\n", bindingSiteSequence, bindingSite.FivePrimeIdx, bindingSite.ThreePrimeIdx)
	fmt.Printf("mfe: %v\n", bindingSite.MinimumFreeEnergy)
	fmt.Printf("translationInitiationRate: %v\n", bindingSite.TranslationInitiationRate)
	// Output:
	// binding site sequence: GGUAAAAAAUG (27,38)
	// mfe: -18.2
	// translationInitiationRate: 1022.4939796226361
}
