package rbs_calculator

import (
	"fmt"
	"math"
	"strings"

	"github.com/TimothyStiles/poly"
	"github.com/TimothyStiles/poly/mfe"
)

/*****************************************************
June, 5, 2021

The Ribosome Binding Site Calculator

Process to develop the model

Current research provides intuition as to the process behind a ribosome binding
to a mRNA, but there is no concensus on an actual model of this interaction.

Thus, to develop a model, we have created a Go Jupyter notebook that allows
one to examine different primary and secondary structure features to come up
with variables and their relationship to the end result of translation initation
rate.

To create the model, we follow a Data -> Model Theory, Variables & Assumptions
-> Visualization -> Coefficient Calculation & Model Equation workflow.


Data:
We need a high quality dataset that includes the mRNA sequence, its structure,
the organism the sequence was inserted into and the translation initiation rate.
The mRNA sequence is split into its 5' untranslated region (5' UTR) and its
protein coding sequence (CDS).

Unfortunately, current datasets don't include the mRNA's structure so we have to
estimate its structure using folding algorithms like LinearFold.

Given that the mRNA's structure plays an important role in the process of a
ribosome binding to a sequence, using a different folding algorithm (such as
ViennaRNA's folding algorithm) will result in a completely different model due
to a different secondary structure.

To combat this, we treat the folded structure as non-ground truth data and
add the estimated folded structure to the dataset above using xxxx.go. This
creates a csv xxxx.go in xxx folder which is then used as the dataset for the
subsequent steps.

Note that since our ribosome binding site model depends on using LinearFold as
the folding algorithm, we only allow a user to input the mRNA sequence, but
not its structure (which is internally computed using LinearFold).

If the folding alogrihm is changed, the xxx.csv dataset must be re-genrated and
all the following steps in the process need to be re-run to generate an updated
model.


Model Theory & Variables:

This model assumeThe variables derived from the theory is put into square
brackets below the relevant paragraph.

We think that ribosome binding occurs in two steps:
1. A ribosome "loosely-binds" to a mRNA.
2. The ribosome then fully binds to a mRNA when a start codon is encountered and
starts translating a protein coding sequence until a terminator codon is
encountered.

A ribosome first "loosely-binds" when the protein S1 on the ribosome encounters
an pyrimidine-rich (AU, or CA) single-stranded region of the mRNA which located
approximately 11nt upstream[123] of the Shine Dalgarno sequence (explained
below). The ribosome then "waits" at the standby site until the mRNA structure
briefly unfolds. When a start codon is encountered, the ribosome then "fully"
binds to the mRNA.
[length of As, Cs and Us in single-stranded region]

The structure of the ribosome is such that its fMet-tRNA (3'-UAC-5') and its 3'
end of the 16S rRNA (3'-ATTCCTCCA-5'), called the anti-Shine Dalgarno sequence,
are approximately 4-9 nucleotides apart. Many studies suggest the presence of
the Shine Dalgarno sequence 4-9 nucleotides upstream (in the 5' direction) of a
start codon in the mRNA is responsible for the full binding of the ribosome to
the mRNA. However, numerous prokaryotic mRNAs don't possess a SD sequence at
all. Hence we don't penalise for the presence or lack of SD sequence. But many
papers do suggest that the presence of long SD sequences actually inhibit
translation due to higher binding affinity rendering the ribosome unable to move
in along the mRNA seqeunce. Thus, we penalise for long SD sequence.
[length of Shine Dalgarno sequence]

Note that the ribosome protein S1 is only present for gram-negative bacteria.
Thus, the binding mechanism for gram-positive species will use a different model
than the one used for this RBS calculator.

123 - https://academic.oup.com/nar/article/46/8/4188/4840235#116447296


Model Assumptions:

We assume that
* if a mRNA contains more than one coding sequences, only one of them is
translated into a protein. Thus, if there are two start codons in a sequence and
a ribosome binds to the mRNA sequence before the first start codon, we assume
that the coding sequence encoded by the second start codon is not translated.
* a ribosome binds randomly to a mRNA. Thus, the ribosome calculator returns a
slice of all possible binding sites.


Visualization:

To understand a relationship between a variable and the results, we have a Go
Jupyter notebook that creates graphs between a variable and the translated
protein levels. By trial and error, we can see the relationship between many
different types of properties of a mRNA sequence and the levels of translated
protein.

Once we can see a (linear, quadratic, cubic) relationship between a variable
and the outcome, we add the variable to our model and calculate its coefficients.


Coefficient Calculation & Model Equation:

Coefficient calculation is done in the same jupyter notebook as the
visualization.

Although the previous step gives us the relationships between individual
variables and the outcome, it could be the case that adding all the variables
together leads to worse prediction of results.

Thus, as a final step we run an analysis that tries out all the possible
combinations of variables in the model to see if removing some terms could
actually improve the predicted results. The least set of variables with the
highest correlation with actual results are then used in the model.


*****************************************************/

type BindingSite struct {
	TranslationInitiationRate, MinimumFreeEnergy float64
	FivePrimeIdx, ThreePrimeIdx                  int
	RRNA, MRNA, MRNAStructure                    string
	MRNAFreeEnergy, MRNARRNAFreeEnergy           float64
	FoundBindingSite                             bool
}

var EColiRNA = "ACCTCCTTA"
var SalisRNA = "ATTCCTCCA"

var defaultLenBindingSite = 11

// RBS calculations according to Keoni's ideas
func RibosomeBindingSiteV1(rRNA, mRNA string) (*BindingSite, error) {
	mRNA = strings.ToUpper(mRNA)
	lenMRNA := len(mRNA)
	if lenMRNA < defaultLenBindingSite {
		return nil, fmt.Errorf("length of mRNA (%v) must be greater than 11", lenMRNA)
	}
	mRNA = strings.ReplaceAll(mRNA, "U", "T")

	var bindingSiteThreePrime int = defaultLenBindingSite
	var mRNArRNAFreeEnergy float64 = 0
	firstBindingSite := true
	lookupTable, _ := LookupTable()
	for i := defaultLenBindingSite; i <= lenMRNA; i++ {
		potentialBindingSite := mRNA[i-defaultLenBindingSite : i]
		// fmt.Println(i, potentialBindingSite)
		bindingSiteFreeEnergy, ok := lookupTable[rRNA][potentialBindingSite]
		if ok {
			if bindingSiteFreeEnergy < mRNArRNAFreeEnergy || firstBindingSite {
				mRNArRNAFreeEnergy = bindingSiteFreeEnergy
				bindingSiteThreePrime = i
				firstBindingSite = false
			}
		}

	}

	// Special case is cannot RBS that drags into the CDS
	// The third variable they use on
	mRNA = strings.ReplaceAll(mRNA, "T", "U")
	// rRNA = strings.ReplaceAll(rRNA, "T", "U")

	// mRNArRNAStructure, _ := poly.LinearFold(mRNA+rRNA, poly.DefaultBeamSize)
	// mRNArRNAFreeEnergy, _, _ := mfe.MinimumFreeEnergy(mRNA+rRNA, mRNArRNAStructure, mfe.DefaultTemperature)

	mRNAStructure, _ := poly.LinearFold(mRNA, poly.DefaultBeamSize)
	mRNAFreeEnergy, _, _ := mfe.MinimumFreeEnergy(mRNA, mRNAStructure, mfe.DefaultTemperature)

	// totalMFE := leastBindingSiteFreeEnergy + mRNAFreeEnergy
	// ΔG_total=ΔG_standby+ΔG_mRNA−rRNA+ΔG_spacing+ΔG_start+ΔG_stacking−ΔG_mRNA
	totalMFE := mRNArRNAFreeEnergy - mRNAFreeEnergy
	return &BindingSite{
		TranslationInitiationRate: translationInitiationRate(totalMFE),
		MinimumFreeEnergy:         totalMFE,
		ThreePrimeIdx:             bindingSiteThreePrime,
		FivePrimeIdx:              bindingSiteThreePrime - 11,
		RRNA:                      rRNA,
		MRNA:                      mRNA,
		MRNAFreeEnergy:            mRNAFreeEnergy,
		MRNAStructure:             mRNAStructure,
		MRNARRNAFreeEnergy:        mRNArRNAFreeEnergy,
		FoundBindingSite:          !firstBindingSite,
	}, nil
}

// RBS calculations that calculates dG_mRNA_rRNA by co-folding the mRNA and
// rRNA sequences
func RibosomeBindingSiteV2(rRNA, mRNA string) (*BindingSite, error) {
	mRNA = strings.ToUpper(mRNA)
	// Special case is cannot RBS that drags into the CDS
	// The third variable they use on
	mRNA = strings.ReplaceAll(mRNA, "T", "U")
	rRNA = strings.ReplaceAll(rRNA, "T", "U")

	mRNArRNAStructure, _ := poly.LinearFold(mRNA+rRNA, poly.DefaultBeamSize)
	mRNArRNAFreeEnergy, _, _ := mfe.MinimumFreeEnergy(mRNA+rRNA, mRNArRNAStructure, mfe.DefaultTemperature)

	mRNAStructure, _ := poly.LinearFold(mRNA, poly.DefaultBeamSize)
	mRNAFreeEnergy, _, _ := mfe.MinimumFreeEnergy(mRNA, mRNAStructure, mfe.DefaultTemperature)

	// ΔG_total=ΔG_standby+ΔG_mRNA−rRNA+ΔG_spacing+ΔG_start+ΔG_stacking−ΔG_mRNA
	totalMFE := mRNArRNAFreeEnergy - mRNAFreeEnergy
	return &BindingSite{
		TranslationInitiationRate: translationInitiationRate(totalMFE),
		MinimumFreeEnergy:         totalMFE,
		RRNA:                      rRNA,
		MRNA:                      mRNA,
		MRNAFreeEnergy:            mRNAFreeEnergy,
		MRNAStructure:             mRNAStructure,
		MRNARRNAFreeEnergy:        mRNArRNAFreeEnergy,
		// irrelevant fields
		ThreePrimeIdx:    0,
		FivePrimeIdx:     0,
		FoundBindingSite: true,
	}, nil
}

func translationInitiationRate(freeEnergy float64) float64 {
	probabilisticFreeEnergyAtEquibiliriumConditionsConverstionFactor := 0.45
	freeEnergyAtEquibiliriumConditions := -probabilisticFreeEnergyAtEquibiliriumConditionsConverstionFactor * freeEnergy
	return math.Exp(freeEnergyAtEquibiliriumConditions)
}
