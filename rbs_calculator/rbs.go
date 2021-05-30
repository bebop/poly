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
	MRNAFreeEnergy, BindingSiteFreeEnergy        float64
}

var EColiRNA = "ACCTCCTTA"
var defaultLenBindingSite = 11

// This file goes agains Go's style conventions and uses the snake case
// convention for variable names (functions, types, and struct fields still
// follow Go's style conventions).
func RibosomeBindingSite(rRNA, mRNA string) (*BindingSite, error) {
	mRNA = strings.ToUpper(mRNA)
	lenMRNA := len(mRNA)
	if lenMRNA < defaultLenBindingSite {
		return nil, fmt.Errorf("length of mRNA (%v) must be greater than 11", lenMRNA)
	}

	var bindingSiteThreePrime int = defaultLenBindingSite
	var leastBindingSiteFreeEnergy float64 = 1000000.0
	lookupTable, _ := LookupTable()
	for i := defaultLenBindingSite; i <= lenMRNA; i++ {
		potentialBindingSite := mRNA[i-defaultLenBindingSite : i]
		fmt.Printf("%v, %v\n", rRNA, potentialBindingSite)
		bindingSiteFreeEnergy, ok := lookupTable[rRNA][potentialBindingSite]
		if potentialBindingSite == "GGTAAAAAATG" {
			panic(fmt.Sprintf("%v, %v", bindingSiteFreeEnergy, ok))
		}
		if ok {
			if bindingSiteFreeEnergy < leastBindingSiteFreeEnergy {
				leastBindingSiteFreeEnergy = bindingSiteFreeEnergy
				bindingSiteThreePrime = i
			}
		}
		// fmt.Println(i)
	}

	// Special case is cannot RBS that drags into the CDS
	// The third variable they use on
	mRNA = strings.ReplaceAll(mRNA, "T", "U")
	mRNAStructure, _ := poly.LinearFold(mRNA, poly.DefaultBeamSize)
	mRNAFreeEnergy, _, _ := mfe.MinimumFreeEnergy(mRNA, mRNAStructure, mfe.DefaultTemperature)

	totalMFE := leastBindingSiteFreeEnergy + mRNAFreeEnergy
	return &BindingSite{
		TranslationInitiationRate: translationInitiationRate(mRNAFreeEnergy),
		MinimumFreeEnergy:         totalMFE,
		ThreePrimeIdx:             bindingSiteThreePrime,
		FivePrimeIdx:              bindingSiteThreePrime - 11,
		RRNA:                      rRNA,
		MRNA:                      mRNA,
		MRNAFreeEnergy:            mRNAFreeEnergy,
		BindingSiteFreeEnergy:     leastBindingSiteFreeEnergy,
	}, nil
}

func translationInitiationRate(freeEnergy float64) float64 {
	probabilisticFreeEnergyAtEquibiliriumConditionsConverstionFactor := 0.45
	freeEnergyAtEquibiliriumConditions := -probabilisticFreeEnergyAtEquibiliriumConditionsConverstionFactor * freeEnergy
	return math.Exp(freeEnergyAtEquibiliriumConditions)
}
