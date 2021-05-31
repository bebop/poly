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
	RRNA, MRNA, MRNAStructure                    string
	MRNAFreeEnergy, BindingSiteFreeEnergy        float64
	FoundBindingSite                             bool
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
	mRNA = strings.ReplaceAll(mRNA, "U", "T")

	var bindingSiteThreePrime int = defaultLenBindingSite
	var leastBindingSiteFreeEnergy float64 = 0
	firstBindingSite := true
	lookupTable, _ := LookupTable()
	for i := defaultLenBindingSite; i <= lenMRNA; i++ {
		potentialBindingSite := mRNA[i-defaultLenBindingSite : i]
		// fmt.Println(i, potentialBindingSite)
		bindingSiteFreeEnergy, ok := lookupTable[rRNA][potentialBindingSite]
		if ok {
			if bindingSiteFreeEnergy < leastBindingSiteFreeEnergy || firstBindingSite {
				leastBindingSiteFreeEnergy = bindingSiteFreeEnergy
				bindingSiteThreePrime = i
				firstBindingSite = false
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
		MRNAStructure:             mRNAStructure,
		BindingSiteFreeEnergy:     leastBindingSiteFreeEnergy,
		FoundBindingSite:          !firstBindingSite,
	}, nil
}

func translationInitiationRate(freeEnergy float64) float64 {
	probabilisticFreeEnergyAtEquibiliriumConditionsConverstionFactor := 0.45
	freeEnergyAtEquibiliriumConditions := -probabilisticFreeEnergyAtEquibiliriumConditionsConverstionFactor * freeEnergy
	return math.Exp(freeEnergyAtEquibiliriumConditions)
}
