package poly

import (
	"index/suffixarray"
	"math"
	"sort"
	"strings"
)

// For reference: https://www.sigmaaldrich.com/technical-documents/articles/biology/oligos-melting-temp.html

// thermodynamics stores enthalpy (dH, kcal/mol) and entropy (dS, cal/mol-K) values for nucleotide pairs
type thermodynamics struct{ H, S float64 }

/******************************************************************************
This section contains various penalties applied when calculating primer melting
temperature using the SantaLucia algorithm.
******************************************************************************/

// penalties for nearest neighbor effects
var nearestNeighborsThermodynamics = map[string]thermodynamics{
	"AA": {-7.6, -21.3},
	"TT": {-7.6, -21.3},
	"AT": {-7.2, -20.4},
	"TA": {-7.2, -21.3},
	"CA": {-8.5, -22.7},
	"TG": {-8.5, -22.7},
	"GT": {-8.4, -22.4},
	"AC": {-8.4, -22.4},
	"CT": {-7.8, -21.0},
	"AG": {-7.8, -21.0},
	"GA": {-8.2, -22.2},
	"TC": {-8.2, -22.2},
	"CG": {-10.6, -27.2},
	"GC": {-9.8, -24.4},
	"GG": {-8.0, -19.9},
	"CC": {-8.0, -19.9},
}

var initialThermodynamicPenalty = thermodynamics{0.2, -5.7}   // penalty for initiating helix
var symmetryThermodynamicPenalty = thermodynamics{0, -1.4}    // penalty for self-complementarity
var terminalATThermodynamicPenalty = thermodynamics{2.2, 6.9} // penalty for 3' AT

/******************************************************************************
End of melting temp penalties section for SantaLucia melting temp algorithm.
******************************************************************************/

// SantaLucia calculates the melting point of a short DNA sequence (15-200 bp), using the Nearest Neighbors method [SantaLucia, J. (1998) PNAS, doi:10.1073/pnas.95.4.1460]
func SantaLucia(sequence string, primerConcentration, saltConcentration, magnesiumConcentration float64) (meltingTemp, dH, dS float64) {
	sequence = strings.ToUpper(sequence)

	const gasConstant = 1.9872 // gas constant (cal / mol - K)

	var symmetryFactor float64 // symmetry factor

	// apply initialization penalty
	dH += initialThermodynamicPenalty.H
	dS += initialThermodynamicPenalty.S
	// apply symmetry penalty if sequence is self-complementary
	if sequence == ReverseComplement(sequence) {
		dH += symmetryThermodynamicPenalty.H
		dS += symmetryThermodynamicPenalty.S
		symmetryFactor = 1
	} else {
		symmetryFactor = 4
	}
	// apply penalty if 3' nucleotides are A or T
	if sequence[len(sequence)-1] == 'A' || sequence[len(sequence)-1] == 'T' {
		dH += terminalATThermodynamicPenalty.H
		dS += terminalATThermodynamicPenalty.S
	}
	// apply salt penalty ; von Ahsen et al 1999
	saltEffect := saltConcentration + (magnesiumConcentration * 140)
	dS += (0.368 * float64(len(sequence)-1) * math.Log(saltEffect))
	// calculate penalty for nearest neighbor effects
	for i := 0; i+1 < len(sequence); i++ {
		dT := nearestNeighborsThermodynamics[sequence[i:i+2]]
		dH += dT.H
		dS += dT.S
	}

	meltingTemp = dH*1000/(dS+gasConstant*math.Log(primerConcentration/symmetryFactor)) - 273.15
	return meltingTemp, dH, dS
}

// MarmurDoty calculates the melting point of an extremely short DNA sequence (<15 bp) using a modified Marmur Doty formula [Marmur J & Doty P (1962). Determination of the base composition of deoxyribonucleic acid from its thermal denaturation temperature. J Mol Biol, 5, 109-118.]
func MarmurDoty(sequence string) float64 {
	sequence = strings.ToUpper(sequence)

	aCount := float64(strings.Count(sequence, "A"))
	tCount := float64(strings.Count(sequence, "T"))
	cCount := float64(strings.Count(sequence, "C"))
	gCount := float64(strings.Count(sequence, "G"))

	meltingTemp := 2*(aCount+tCount) + 4*(cCount+gCount) - 7.0
	return meltingTemp
}

// MeltingTemp calls SantaLucia with default inputs for primer and salt concentration.
func MeltingTemp(sequence string) float64 {
	primerConcentration := 500e-9 // 500 nM (nanomolar) primer concentration
	saltConcentration := 50e-3    // 50 mM (millimolar) sodium concentration
	magnesiumConcentration := 0.0 // 0 mM (millimolar) magnesium concentration

	meltingTemp, _, _ := SantaLucia(sequence, primerConcentration, saltConcentration, magnesiumConcentration)
	return meltingTemp
}

/******************************************************************************

Start of PCR functions

The PCR functions basically take a sequence and simulate PCRs. The PCR template
is first converted to a suffixarray to enable fast searching. In the case of
small templates, suffixarray generation does not take a significant amount of time.
In the case of large templates with many PCRs (such as PCRs from an oligo pool),
the suffixarray saves a significant amount of time when simulating. A large template
with a single PCR may take longer to simulate, but will still work fine.

******************************************************************************/

// MakePrimers makes primers for a given sequence at a targetTm. It uses the default values from MeltingTemp
func MakePrimers(sequence string, targetTm float64) (string, string) {
	minimalPrimerLength := 15
	forwardPrimer := sequence[0:minimalPrimerLength]
	for i := 0; MeltingTemp(forwardPrimer) < targetTm; i++ {
		forwardPrimer = sequence[0 : minimalPrimerLength+i]
	}
	reversePrimer := ReverseComplement(sequence[len(sequence)-minimalPrimerLength:])
	for i := 0; MeltingTemp(reversePrimer) < targetTm; i++ {
		reversePrimer = ReverseComplement(sequence[len(sequence)-(minimalPrimerLength+i):])
	}
	return forwardPrimer, reversePrimer
}

// Pcr simulates a PCR reaction in both the forward and reverse direction
func Pcr(sequence string, circular bool, targetTm float64, forward string, reverse string) []string {
	// Setup forward and reverse indexes for this PCR reaction
	indexForward := suffixarray.New([]byte(sequence))
	sequenceReverse := ReverseComplement(sequence)
	indexReverse := suffixarray.New([]byte(sequenceReverse))

	// Run the PCR reaction with PcrSa
	return PcrSa(sequence, indexForward, sequenceReverse, indexReverse, circular, targetTm, forward, reverse)
}

// PcrSa simulates a PCR reaction, given suffix arrays of the forward sequence and reverse sequence of a given DNA molecule
func PcrSa(sequenceForward string, indexForward *suffixarray.Index, sequenceReverse string, indexReverse *suffixarray.Index, circular bool, targetTm float64, forward string, reverse string) []string {
	return append(PcrSaUnidirectional(sequenceForward, indexForward, circular, targetTm, forward, reverse), PcrSaUnidirectional(sequenceReverse, indexReverse, circular, targetTm, forward, reverse)...)
}

// PcrSaUnidirectional simulates a PCR reaction with a prebuilt suffix array (https://golang.org/pkg/index/suffixarray/) in a single direction
func PcrSaUnidirectional(sequence string, sequenceIndex *suffixarray.Index, circular bool, targetTm float64, forward string, reverse string) []string {
	// For the forward primer, index from right to left until the targetTm is hit
	var minimalForwardLength int
	for i := 10; MeltingTemp(forward[len(forward)-i:]) < targetTm; i++ {
		minimalForwardLength = i + 1
		if forward[len(forward)-i:] == forward {
			return []string{}
		}
	}
	// Once the targetTm is hit, search for the sequence in the sequenceIndex
	forwardSequence := forward[len(forward)-minimalForwardLength:]
	forwardLocations := sequenceIndex.Lookup([]byte(forwardSequence), -1)

	// Do the same for the reverse primer
	var minimalReverseLength int
	for i := 10; MeltingTemp(reverse[:i]) < targetTm; i++ {
		minimalReverseLength = i + 1
		if reverse == reverse[:i] {
			return []string{}
		}
	}
	reverseSequence := ReverseComplement(reverse[:minimalReverseLength])
	reverseLocations := sequenceIndex.Lookup([]byte(reverseSequence), -1)

	// Now that there are forwardLocations and reverseLocations, we need to find forward primer / reverse primer pairs. We do this
	// by finding a reverse primer after a forward primer that is smaller than the next forward pair. Once at the end of the list,
	// we need to check circularity, and if we are working with a circular genome, check that there aren't primers overlapping at
	// the origin of the sequence
	var pcrFragment string
	var pcrFragments []string

	// First, we sort the forwardLocations and reverseLocations lists
	sort.Ints(forwardLocations)
	sort.Ints(reverseLocations)

	// Next, iterate through the forwardLocations list
	for i, forwardLocation := range forwardLocations {
		// First, make sure that this isn't the last element in forwardLocations
		if i+1 != len(forwardLocations) {
			// If this isn't the last element in forwardLocations, then we can select the first reverseLocation that is less than the next forwardLocation
			for _, reverseLocation := range reverseLocations {
				if (forwardLocation < reverseLocation) && (reverseLocation < forwardLocations[i+1]) {
					// If both are true, we have found the sequence we are aiming to PCR!
					pcrFragment = forward[:len(forward)-minimalForwardLength] + sequence[forwardLocation:reverseLocation] + ReverseComplement(reverse)
					// Append to pcrFragments
					pcrFragments = append(pcrFragments, pcrFragment)
					break
				}
			}
		} else {
			foundFragment := false
			for _, reverseLocation := range reverseLocations {
				if forwardLocation < reverseLocation {
					pcrFragment = forward[:len(forward)-minimalForwardLength] + sequence[forwardLocation:reverseLocation] + ReverseComplement(reverse)
					pcrFragments = append(pcrFragments, pcrFragment)
					foundFragment = true
				}
			}
			// If the sequence is circular and we haven't found a fragment yet, check the other side of the origin
			if (circular == true) && (foundFragment == false) {
				for _, reverseLocation := range reverseLocations {
					if forwardLocations[0] > reverseLocation {
						// If either one of these are true, create a new pcrFragment and append to pcrFragments
						pcrFragment = forward[:len(forward)-minimalForwardLength] + sequence[forwardLocation:] + sequence[:reverseLocation] + ReverseComplement(reverse)
						pcrFragments = append(pcrFragments, pcrFragment)
					}
				}
			}
		}
	}
	return pcrFragments
}
