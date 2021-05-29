package rbs_calculator

import (
	"fmt"

	"github.com/TimothyStiles/poly"
	"github.com/TimothyStiles/poly/mfe"
)

type BindingSite struct {
	translationInitiationRate, minimumFreeEnergy float64
	fivePrimeIdx, threePrimeIdx                  int
	rrna, mrna                                   string
}

// This file goes agains Go's style conventions and uses the snake case
// convention for variable names (functions, types, and struct fields still
// follow Go's style conventions).
func RibosomeBindingSite(mRNA string) (*BindingSite, error) {
	len_mRNA := len(mRNA)
	if len_mRNA < 11 {
		return nil, fmt.Errorf("length of mRNA (%v) must be greater than 11", len_mRNA)
	}

	rRNA := "ACCTCCTTA"

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

	mRNA_structure, _ := poly.LinearFold(mRNA)
	dG_mRNA, _, _ := mfe.MinimumFreeEnergy(mRNA, mRNA_structure, mfe.DefaultTemperature)

	total_mfe := least_dG_rRNA_mRNA + dG_mRNA
	return &BindingSite{
		minimumFreeEnergy: total_mfe,
		threePrimeIdx:     binding_site_three_prime,
		fivePrimeIdx:      binding_site_three_prime - 11,
		rrna:              rRNA,
		mrna:              mRNA,
	}, nil
}
