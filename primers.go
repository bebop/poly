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
May 23 2021

Start of the De Bruijn stuff

=== Barcode basics ===

We're rapidly getting better at sequencing a lot of DNA. At their core, most
DNA sequencing technologies pool together many samples and sequence them all
at once. For example, let's say we have 2 samples of DNA whose true sequence
is as follows:

DNA-1 := ATGC
DNA-2 := AGGC

If we pooled these two samples together into a single tube, and sequenced
them, we would not be able to tell if ATGC came from DNA-1 or DNA-2. In order
to tell the difference, we would have to go through the process of DNA
barcoding. Let's attach two small barcodes to each DNA fragment separately
in their own tubes and then pool them togehter:

Barcode-1 + DNA-1 = GC + ATGC = GCATGC
Barcode-2 + DNA-2 = AT + AGGC = ATAGGC

When we sequence this pool together, we will end up with two sequences,
GCATGC and ATAGGC. If we correlate the first 2 base pairs with the tube
the sample came from, we can derive DNA-1 is ATGC and DNA-2 is AGGC.



=== Redundancy and start sites ===

Now, let's say we have the need for N number of samples to be pooled
together. The minimal barcode length could be expressed as:

n = number of samples
b = bases required in a minimal barcode

4^b = n
OR
log4n = b

In our perfect case, we would only need 8 base pair barcodes to represent 65536
different samples. Reality is a little different, however.

1. Failure of DNA sequencers to accurately sequence the barcode, leading for
   one barcoed to be mistaken for a different barcode
2. Misalignment of the barcode to the sequence. We cannot guarantee that the
   DNA sequencer will begin sequencing our fragment at an exact base pair.
3. Misreading of sequence as barcode. If our barcode is only 8 base pairs, on
   average, it will occur once within 65536 base pairs, and that occurrence may
   be misread as a barcode.

These challenges force us to build a barcode that has the following features:

1. Any barcode must be different enough from any other barcode that there will
   be no misreading, even with mutated base pairs.
2. Any barcode must be large enough that, on average, it will not occur in a
   natural piece of DNA.

While the second feature is quite easy (use ~20-30 base pair barcodes), the
first can be challenging. When developing a large quantity of barcodes, how
do you guarantee that they are optimally distanced from each other so that
there will be no cross-talk?

=== Our solution to distanced barcodes ===

De Bruijn sequences are an interesting data structure where every possible
substring of length N occurs exactly once as a substring(1). For example, a De
Bruijn sequence of length 3 will only have ATG occur once in the entire
sequence.

By constructing a nucleobase De Bruijn sequence, and selecting barcodes from
within that De Bruijn sequence, we can guarantee that each barcode will never
share any N length substring, since it only occurs once within the whole De
Bruijn sequence.

For example, a nucleobase De Bruijn sequence of substring length 6 is 4101
base pairs long (4^n + (n-1)). You can generate 205 20 base pair barcodes with
each barcode guaranteed to never share any 6 base pairs. This makes it very
easy to unambiguously parse which samples came from where, while maintaining
a guarantee of optimal distancing between your barcodes.

Good luck with barcoding,

Keoni

(1) https://en.wikipedia.org/wiki/De_Bruijn_sequence

******************************************************************************/

// NucleobaseDeBruijnSequence generates a DNA DeBruijn sequence with alphabet ATGC. DeBruijn sequences are basically a string with all unique substrings of an alphabet represented exactly once. Code is adapted from https://rosettacode.org/wiki/De_Bruijn_sequences#Go
func NucleobaseDeBruijnSequence(substringLength int) string {
	alphabet := "ATGC"
	alphabetLength := len(alphabet)
	a := make([]byte, alphabetLength*substringLength)
	var seq []byte
	// The following function is mainly adapted from rosettacode.
	var ConstructDeBruijn func(int, int) // recursive closure
	ConstructDeBruijn = func(t, p int) {
		if t > substringLength {
			if substringLength%p == 0 {
				seq = append(seq, a[1:p+1]...)
			}
		} else {
			a[t] = a[t-p]
			ConstructDeBruijn(t+1, p)
			for j := int(a[t-p] + 1); j < alphabetLength; j++ {
				a[t] = byte(j)
				ConstructDeBruijn(t+1, t)
			}
		}
	}
	ConstructDeBruijn(1, 1)
	var buf bytes.Buffer
	for _, i := range seq {
		buf.WriteByte(alphabet[i])
	}
	b := buf.String()
	return b + b[0:substringLength-1] // as cyclic append first (n-1) digits
}

// CreateBarcodes creates a list of barcodes given a desired barcode length, the maxSubSequence shared in each barcode,
// Sequences may be marked as banned by passing a static list, `bannedSequences`, or, if more flexibility is needed, through a list of `bannedFunctions` that dynamically generates bannedSequences.
// If a sequence is banned, it will not appear within a barcode. The a `bannedFunctions` function can determine if a barcode should be banned or not on the fly. If it is banned, we will continuing iterating until a barcode is found that satisfies the bannedFunction requirement.
func CreateBarcodes(length int, maxSubSequence int, bannedSequences []string, bannedFunctions []func(string) bool) []string {
	var barcodes []string
	var start int
	var end int
	debruijn := NucleobaseDeBruijnSequence(maxSubSequence)
	for barcodeNum := 0; (barcodeNum*(length-(maxSubSequence-1)))+length < len(debruijn); {
		start = barcodeNum * (length - (maxSubSequence - 1))
		end = start + length
		barcodeNum++
		for _, bannedSequence := range bannedSequences {
			// If the current deBruijn range has the banned sequence, iterate one base pair ahead. If the iteration reaches the end of the deBruijn sequence, close the channel and return the function.
			for strings.Contains(debruijn[start:end], bannedSequence) {
				if end+1 > len(debruijn) {
					return barcodes
				}
				start++
				end++
				barcodeNum++
			}
			// Check reverse complement as well for the banned sequence
			for strings.Contains(debruijn[start:end], ReverseComplement(bannedSequence)) {
				if end+1 > len(debruijn) {
					return barcodes
				}
				start++
				end++
				barcodeNum++
			}
		}
		for _, bannedFunction := range bannedFunctions {
			// If the function returns False for the deBruijn range, iterate one base pair ahead. If the iteration reaches the end of the deBruijn sequence, close the channel and return the function.
			for !bannedFunction(debruijn[start:end]) {
				if end+1 > len(debruijn) {
					return barcodes
				}
				start++
				end++
				barcodeNum++
			}
		}
		barcodes = append(barcodes, debruijn[start:end])
	}
	return barcodes
}

// SimpleCreateBarcodes is a simplified version of CreateBarcodes with sane defaults.
func SimpleCreateBarcodes(length int, maxSubSequence int) []string {
	return CreateBarcodes(length, maxSubSequence, []string{}, []func(string) bool{})
}
