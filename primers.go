package main

import (
	"fmt"
	"math"
	"strings"
)

// For reference: https://www.sigmaaldrich.com/technical-documents/articles/biology/oligos-melting-temp.html

// thermodynamics stores enthalpy (dH, kcal/mol) and entropy (dS, cal/mol-K) values for nucleotide pairs
type thermodynamics struct{ H, S float64 }

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
}                                                       // penalties for nearest neighbor effects
var initThermodynamics = thermodynamics{0.2, -5.7}      // penalty for initiating helix
var symmetryThermodynamics = thermodynamics{0, -1.4}    // penalty for self-complementarity
var terminalATThermodynamics = thermodynamics{2.2, 6.9} // penalty for 3' AT

//SantaLucia calculates the melting point of a short DNA sequence (15-200 bp), using the Nearest Neighbors method [SantaLucia, J. (1998) PNAS, doi:10.1073/pnas.95.4.1460]
func SantaLucia(seq string, cPrimer, cNa, cMg float64) (Tm, dH, dS float64) {
	seq = strings.ToUpper(seq)

	R := 1.9872 // gas constant (cal / mol - K)

	var x float64 // symmetry factor

	// apply initialization penalty
	dH += initThermodynamics.H
	dS += initThermodynamics.S
	// apply symmetry penalty if seq is self-complementary
	if seq == ReverseComplement(seq) {
		fmt.Println("applying symmetry penalty")
		dH += symmetryThermodynamics.H
		dS += symmetryThermodynamics.S
		x = 1
	} else {
		fmt.Println("skipping symmetry penalty")
		x = 4
	}
	// apply penalty if 3' nucleotides are A or T
	if seq[len(seq)-1] == 'A' || seq[len(seq)-1] == 'T' {
		dH += terminalATThermodynamics.H
		dS += terminalATThermodynamics.S
	}
	// apply salt penalty ; von Ahsen et al 1999
	saltEffect := cNa + (cMg * 140)
	dS += (0.368 * float64(len(seq)-1) * math.Log(saltEffect))
	// calculate penalty for nearest neighbor effects
	for i := 0; i+1 < len(seq); i++ {
		dT := nearestNeighborsThermodynamics[seq[i:i+2]]
		dH += dT.H
		dS += dT.S
	}

	fmt.Printf("dH = %f, dS = %f\n", dH, dS)
	Tm = dH*1000/(dS+R*math.Log(cPrimer/x)) - 273.15
	return Tm, dH, dS
}

// MarmurDoty calculates the melting point of an extremely short DNA sequence (<15 bp) using a modified Marmur Doty formula [Marmur J & Doty P (1962). Determination of the base composition of deoxyribonucleic acid from its thermal denaturation temperature. J Mol Biol, 5, 109-118.]
func MarmurDoty(seq string) float64 {
	seq = strings.ToUpper(seq)

	aCount := float64(strings.Count(seq, "A"))
	tCount := float64(strings.Count(seq, "T"))
	cCount := float64(strings.Count(seq, "C"))
	gCount := float64(strings.Count(seq, "G"))

	Tm := 2*(aCount+tCount) + 4*(cCount+gCount) - 7.0
	return Tm
}

// CalcTM calls SantaLucia with default inputs for primer and salt conc.
func CalcTM(seq string) float64 {
	cPrimer := 500e-9 // 500 nM
	cNa := 50e-3      // 50 mM
	cMg := 0.0        // 0 mM

	Tm, _, _ := SantaLucia(seq, cPrimer, cNa, cMg)
	return Tm
}
