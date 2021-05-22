package poly

import (
	"bytes"
	"math"
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
Start of the De Bruijn stuff
******************************************************************************/

// DeBruijn generates a DNA DeBruijn sequence with alphabet ATGC. DeBruijn sequences are basically a string with all unique substrings of an alphabet represented exactly once.
// https://en.wikipedia.org/wiki/De_Bruijn_sequence
// https://rosettacode.org/wiki/De_Bruijn_sequences#Go
// Pulled and adapted from here
func DeBruijn(n int) string {
	k := 4
	alphabet := "ATGC"
	a := make([]byte, k*n)
	var seq []byte
	var db func(int, int) // recursive closure
	db = func(t, p int) {
		if t > n {
			if n%p == 0 {
				seq = append(seq, a[1:p+1]...)
			}
		} else {
			a[t] = a[t-p]
			db(t+1, p)
			for j := int(a[t-p] + 1); j < k; j++ {
				a[t] = byte(j)
				db(t+1, t)
			}
		}
	}
	db(1, 1)
	var buf bytes.Buffer
	for _, i := range seq {
		buf.WriteByte(alphabet[i])
	}
	b := buf.String()
	return b + b[0:n-1] // as cyclic append first (n-1) digits
}

// CreateBarcodes creates a list of barcodes given a desired barcode length, the maxSubSequence shared in each barcode,
// any banned sequences within the barcode, and any functions which may ban sequences within the barcode.
func CreateBarcodes(length int, maxSubSequence int, banned []string, bannedFunc []func(string) bool) []string {
	var barcodes []string
	var start int
	var end int
	debruijn := DeBruijn(maxSubSequence)
	for i := 0; (i*(length-(maxSubSequence-1)))+length < len(debruijn); {
		start = i * (length - (maxSubSequence - 1))
		end = start + length
		i++
		for _, b := range banned {
			// If the current deBruijn range has the banned sequence, iterate one base pair ahead. If the iteration reaches the end of the deBruijn sequence, close the channel and return the function.
			for strings.Contains(debruijn[start:end], b) {
				if end+1 > len(debruijn) {
					return barcodes
				}
				start++
				end++
				i++
			}
		}
		for _, f := range bannedFunc {
			// If the function returns False for the deBruijn range, iterate one base pair ahead. If the iteration reaches the end of the deBruijn sequence, close the channel and return the function.
			for !f(debruijn[start:end]) {
				if end+1 > len(debruijn) {
					return barcodes
				}
				start++
				end++
				i++
			}
		}
		barcodes = append(barcodes, debruijn[start:end])
	}
	return barcodes
}
