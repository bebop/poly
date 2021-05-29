package rbs_calculator

import (
	"fmt"
	"math"
	"strings"

	"github.com/TimothyStiles/poly"
	"github.com/TimothyStiles/poly/mfe"
)

type BindingSite struct {
	TranslationInitiationRate, MinimumFreeEnergy float64
	FivePrimeIdx, ThreePrimeIdx                  int
	RRNA, MRNA                                   string
}

var EColiRNA = "ACCTCCTTA"

// This file goes agains Go's style conventions and uses the snake case
// convention for variable names (functions, types, and struct fields still
// follow Go's style conventions).
func RibosomeBindingSite(rRNA, mRNA string) (*BindingSite, error) {
	len_mRNA := len(mRNA)
	if len_mRNA < 11 {
		return nil, fmt.Errorf("length of mRNA (%v) must be greater than 11", len_mRNA)
	}

	var binding_site_three_prime int
	var least_dG_rRNA_mRNA float64 = 1000000.0
	lookup_table, _ := LookupTable()
	for i := 11; i <= len_mRNA; i++ {
		potential_binding_site := mRNA[i-11 : i]
		dG_rRNA_mRNA := lookup_table[rRNA][potential_binding_site]
		if dG_rRNA_mRNA < least_dG_rRNA_mRNA {
			least_dG_rRNA_mRNA = dG_rRNA_mRNA
			binding_site_three_prime = i
		}
	}

	// Special case is cannot RBS that drags into the CDS
	// The third variable they use on

	mRNA = strings.ReplaceAll(mRNA, "T", "U")
	mRNA_structure, _ := poly.LinearFold(mRNA, poly.DefaultBeamSize)
	// log.Printf("sequence: %v\n", mRNA)
	// log.Printf("structure: %v\n", mRNA_structure)

	dG_mRNA, _, _ := mfe.MinimumFreeEnergy(mRNA, mRNA_structure, mfe.DefaultTemperature)
	translationInitiationRate := translationInitiationRate(dG_mRNA)

	total_mfe := least_dG_rRNA_mRNA + dG_mRNA
	return &BindingSite{
		TranslationInitiationRate: translationInitiationRate,
		MinimumFreeEnergy:         total_mfe,
		ThreePrimeIdx:             binding_site_three_prime,
		FivePrimeIdx:              binding_site_three_prime - 11,
		RRNA:                      rRNA,
		MRNA:                      mRNA,
	}, nil
}

func translationInitiationRate(freeEnergy float64) float64 {
	probabilisticFreeEnergyAtEquibiliriumConditionsConverstionFactor := 0.45
	freeEnergyAtEquibiliriumConditions := -probabilisticFreeEnergyAtEquibiliriumConditionsConverstionFactor * freeEnergy
	return math.Exp(freeEnergyAtEquibiliriumConditions)
}
