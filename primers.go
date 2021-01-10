package poly

import (
	"errors"
	"math"
	"strconv"
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

// TODO make custom function for phusion according to https://tmcalculator.neb.com/#!/help
// MeltingTemp calls SantaLucia with default inputs for primer and salt concentration.
func MeltingTemp(sequence string) float64 {
	primerConcentration := 500e-9 // 500 nM (nanomolar) primer concentration
	saltConcentration := 50e-3    // 50 mM (millimolar) sodium concentration
	magnesiumConcentration := 0.0 // 0 mM (millimolar) magnesium concentration

	meltingTemp, _, _ := SantaLucia(sequence, primerConcentration, saltConcentration, magnesiumConcentration)
	return meltingTemp
}

// MakePrimers makes primers for a given sequence at a targetTm. It uses the default values from MeltingTemp
func MakePrimers(sequence string, targetTm float64) (string, string) {
	forwardPrimer := sequence[0:15]
	for i := 0; MeltingTemp(forwardPrimer) < targetTm; i++ {
		forwardPrimer = sequence[0 : 15+i]
	}
	reversePrimer := ReverseComplement(sequence[len(sequence)-15:])
	for i := 0; MeltingTemp(reversePrimer) < targetTm; i++ {
		reversePrimer = ReverseComplement(sequence[len(sequence)-(15+i):])
	}
	return forwardPrimer, reversePrimer
}

// AminoAcidMutation holds information for a given amino acid mutation. Indexing of positions starts at 1, as per convention
type AminoAcidMutation struct {
	StartPosition     int
	EndPosition       int
	InitialAminoAcids string
	MutatedAminoAcids string
}

// NucleotideMutation holds information for a given nucleotide mutation. Indexing of positions starts at 1, as per convention
type NucleotideMutation struct {
	StartPosition      int
	EndPosition        int
	MutatedNucleotides string
}

// MakeNucleotideMutationPrimers takes a sequence and a given sequence context and generates primers to make a targeted mutation
func MakeNucleotideMutationPrimers(sequence string, sequenceContext string, mutation NucleotideMutation, targetTm float64, overlap int) (string, string, error) {
	// Check that StartPosition < EndPosition
	if mutation.StartPosition > mutation.EndPosition {
		return "", "", errors.New("StartPosition is greater than EndPosition in NucelotideMutation. StartPosition: " + strconv.Itoa(mutation.StartPosition) + " EndPosition: " + strconv.Itoa(mutation.EndPosition))
	}

	// Check that the overall length of the sequence is greater than EndPosition
	if mutation.EndPosition > len(sequence) {
		return "", "", errors.New("EndPosition is greater than length of sequence. EndPosition: " + strconv.Itoa(mutation.EndPosition) + " Sequence length: " + strconv.Itoa(len(sequence)))
	}

	// Check that sequence is within the sequenceContext
	sequenceStart := strings.Index(sequenceContext, sequence)
	if sequenceStart == -1 {
		return "", "", errors.New("sequence is not a substring of sequenceContext")
	}

	// Index the "real" start and end position given the sequenceStart within the sequenceContext.
	// Note: indexing of positions in Mutation notation begins at 1, so 1 is subtracted from the start
	// and end positions
	indexedStartPosition := sequenceStart + (mutation.StartPosition - 1)
	indexedEndPosition := sequenceStart + (mutation.EndPosition - 1)

	// Generate starter primers, in which nucleotide mutations can be added
	seedPrimerFor, _ := MakePrimers(sequenceContext[indexedEndPosition:], targetTm)
	_, seedPrimerRev := MakePrimers(sequenceContext[:indexedStartPosition], targetTm)

	// To assemble properly, the fragments resulting from a PCR with these primers need a certain amount of overlap.
	// This overlap will be the length of the mutation + a little homologous sequence on both ends.
	var additionalForwardOverlap string
	var additionalReverseOverlap string
	additionalOverlap := overlap - len(mutation.MutatedNucleotides)
	if additionalOverlap > 1 {
		// Split overlap between forward and reverse
		targetAdditionalOverlap := additionalOverlap / 2
		additionalForwardOverlap = sequenceContext[indexedStartPosition-targetAdditionalOverlap : indexedStartPosition]
		additionalReverseOverlap = ReverseComplement(sequenceContext[indexedEndPosition : indexedEndPosition+targetAdditionalOverlap])
	}

	// Generate final primers
	primerFor := additionalForwardOverlap + mutation.MutatedNucleotides + seedPrimerFor
	primerRev := additionalReverseOverlap + ReverseComplement(mutation.MutatedNucleotides) + seedPrimerRev

	return primerFor, primerRev, nil
}
