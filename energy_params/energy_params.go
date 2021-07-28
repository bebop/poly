package energy_params

import "embed"

/******************************************************************************

This file defines structs needed to contain information about RNA free energy
params, and the types needed to specify which set of energy parameter values
to parse.

For more information on parsing the energy params, please see `parse.go`.

For more information on scaling the energy params, please see `scale.go`.

******************************************************************************/

const (
	// NbDistinguishableBasePairs is the number of distinguishable base pairs:
	// CG, GC, GU, UG, AU, UA, & non-standard
	NbDistinguishableBasePairs int = 7
	// NbDistinguishableNucleotides is the number of distinguishable nucleotides:
	// A, C, G, U
	NbDistinguishableNucleotides int = 4
	// MaxLenLoop is the maximum length of a loop (hairpin, or multi-loop)
	MaxLenLoop int = 30
	// ZeroCelsiusInKelvin is 0 deg Celsius in Kelvin
	ZeroCelsiusInKelvin float64 = 273.15
)

// NewEnergyParams is a wrapper function that calls the required functions
// which parse the specified energy parameter set and scales it by
// `temperatureInCelsius`.
func NewEnergyParams(energyParamsSet EnergyParamsSet, temperatureInCelsius float64) *EnergyParams {
	return newRawEnergyParams(energyParamsSet).scaleByTemperature(temperatureInCelsius)
}

// EnergyParams contains all the energy parameters needed for the free energy
// calculations.
//
// The order of entries to access theses matrices always uses the closing pair or
// pairs as the first indices followed by the unpaired bases in 5' to 3' direction.
// For example, if we have a 2x2 interior loop:
// ```
// 		      5'-GAUA-3'
// 		      3'-CGCU-5'
// ```
// The closing pairs for the loops are GC and UA (not AU!), and the unpaired bases
// are (in 5' to 3' direction, starting at the first pair) A U C G.
// Thus, the energy for this sequence is:
// ```
// 	pairs:                    GC UA A  U  C  G
// 						interior2x2Loop[1][5][1][4][2][3]
// ```
// (See `BasePairEncodedTypeMap` and `NucleotideEncodedIntMap` for more details
// on how pairs and unpaired nucleotides are encoded)
// Note that this sequence is symmetric so the sequence is equivalent to:
// ```
// 					5'-UCGC-3'
// 					3'-AUAG-5'
// ```
// which means the energy of the sequence is equivalent to:
// ```
// 	pairs:                    UA GC C  G  A  U
// 						interior2x2Loop[5][1][1][4][2][3]
// ```
type EnergyParams struct {

	// The matrix of free energies for stacked pairs, indexed by the two encoded closing
	// pairs. The list should be formatted as symmetric a `7*7`
	// (`NbDistinguishableBasePairs`*`NbDistinguishableBasePairs`) matrix,
	// conforming to the order explained above. As an example the stacked pair
	// ```
	// 					5'-GU-3'
	// 					3'-CA-5'
	// ```
	// corresponds to the entry StackingPair[1][4] (GC=1, AU=4) which should be
	// identical to StackingPair[4][1] (AU=4, GC=1).
	// size: [NbDistinguishableBasePairs][NbDistinguishableBasePairs]int
	StackingPair [][]int
	// Free energies of hairpin loops as a function of size. The list should
	// contain 31 (MaxLenLoop + 1) entries. Since the minimum size of a hairpin loop
	// is 3 and we start counting with 0, the first three values should be `inf` to
	// indicate a forbidden value.
	// size: [MaxLenLoop + 1]int
	HairpinLoop []int
	// Free energies of Bulge loops. Should contain 31 (MaxLenLoop + 1) entries,
	// the first one being `inf`.
	// size: [MaxLenLoop + 1]int
	Bulge []int

	// Free energies of interior loops. Should contain 31 (MaxLenLoop + 1) entries,
	// the first 4 being `inf` (since smaller loops are tabulated).
	//
	// This field was previous called internal_loop, but has been renamed to
	// interior loop to remain consistent with the names of other interior loops
	// size: [MaxLenLoop + 1]int
	InteriorLoop []int

	// Free energies for the interaction between the closing pair of an interior
	// loop and the two unpaired bases adjacent to the helix. This is a three
	// dimensional array indexed by the type of the closing pair and the two
	// unpaired bases. Since we distinguish 4 bases (A, C, G, & U) the list
	// contains `7*5*5` entries. The order is such that for example the mismatch
	//
	// ```
	// 							5'-CU-3'
	// 							3'-GC-5'
	// ```
	// corresponds to entry MismatchInteriorLoop[0][3][1] (CG=0, U=3, C=1).
	// Note that the matrix uses dims of length 5 (instead of 4) for the bases
	// since indexing the energy params uses the `NucleotideEncodedIntMap` map
	// which starts at 1 instead of 0.
	//
	// More information about mismatch energy: https://rna.urmc.rochester.edu/NNDB/turner04/tm.html
	// size: [NbDistinguishableBasePairs][NbDistinguishableNucleotides + 1][NbDistinguishableNucleotides + 1]int
	MismatchInteriorLoop [][][]int
	// size: [NbDistinguishableBasePairs][NbDistinguishableNucleotides + 1][NbDistinguishableNucleotides + 1]int
	Mismatch1xnInteriorLoop [][][]int
	// size: [NbDistinguishableBasePairs][NbDistinguishableNucleotides + 1][NbDistinguishableNucleotides + 1]int
	Mismatch2x3InteriorLoop [][][]int
	// size: [NbDistinguishableBasePairs][NbDistinguishableNucleotides + 1][NbDistinguishableNucleotides + 1]int
	MismatchExteriorLoop [][][]int
	// size: [NbDistinguishableBasePairs][NbDistinguishableNucleotides + 1][NbDistinguishableNucleotides + 1]int
	MismatchHairpinLoop [][][]int
	// size: [NbDistinguishableBasePairs][NbDistinguishableNucleotides + 1][NbDistinguishableNucleotides + 1]int
	MismatchMultiLoop [][][]int

	// Energies for the interaction of an unpaired base on the 5' side and
	// adjacent to a helix in multiloops and free ends (the equivalent of mismatch
	// energies in interior and hairpin loops). The array is indexed by the type
	// of pair closing the helix and the unpaired base and, therefore, forms a `7*5`
	// matrix. For example the dangling base in
	// ```
	// 						5'-GC-3'
	// 						3'- G-5'
	// ```
	// corresponds to entry DanglingEndsFivePrime[0][3] (CG=0, G=3).
	//
	// More information about dangling ends: https://rna.urmc.rochester.edu/NNDB/turner04/de.html
	// size: [NbDistinguishableBasePairs][NbDistinguishableNucleotides + 1]int
	DanglingEndsFivePrime [][]int

	// Same as above for bases on the 3' side of a helix.
	// ```
	// 			       5'- A-3'
	// 			       3'-AU-5'
	// ```
	// corresponds to entry DanglingEndsThreePrime[4][1] (AU=4, A=1).
	// size: [NbDistinguishableBasePairs][NbDistinguishableNucleotides + 1]int
	DanglingEndsThreePrime [][]int

	// Free energies for symmetric size 2 interior loops. `7*7*5*5` entries.
	// Example:
	// ```
	// 						5'-CUU-3'
	// 						3'-GCA-5'
	// ```
	// corresponds to entry Interior1x1Loop[0][4][4][2] (CG=0, AU=4, U=4, C=2),
	// which should be identical to Interior1x1Loop[4][0][2][4] (AU=4, CG=0, C=2, U=4).
	// size: [NbDistinguishableBasePairs][NbDistinguishableBasePairs][NbDistinguishableNucleotides + 1][NbDistinguishableNucleotides + 1]int
	Interior1x1Loop [][][][]int

	// Free energies for 2x1 interior loops, where 2 is the number of unpaired
	// nucleotides on the larger 'side' of the interior loop and 1 is the number of
	// unpaired nucleotides on the smaller 'side' of the interior loop.
	// `7*7*5*5*5` entries.
	// Example:
	// ```
	// 						5'-CUUU-3'
	// 						3'-GC A-5'
	// ```
	// corresponds to entry Interior2x1Loop[0][4][4][4][2] (CG=0, AU=4, U=4, U=4, C=2).
	// Note that this matrix is always accessed in the 5' to 3' direction with the
	// larger number of unpaired nucleotides first.
	// size: [NbDistinguishableBasePairs][NbDistinguishableBasePairs][NbDistinguishableNucleotides + 1][NbDistinguishableNucleotides + 1][NbDistinguishableNucleotides + 1]int
	Interior2x1Loop [][][][][]int

	// Free energies for symmetric size 4 interior loops. `7*7*5*5*5*5` entries.
	// Example:
	// ```
	// 						5'-CUAU-3'
	// 						3'-GCCA-5'
	// ```
	// corresponds to entry Interior2x2Loop[0][4][4][1][2][2] (CG=0, AU=4, U=4, A=1, C=2, C=2),
	// which should be identical to Interior2x2Loop[4][0][2][2][1][4] (AU=4, CG=0, C=2, C=2, A=1, U=4).
	// size: [NbDistinguishableBasePairs][NbDistinguishableBasePairs][NbDistinguishableNucleotides + 1][NbDistinguishableNucleotides + 1][NbDistinguishableNucleotides + 1][NbDistinguishableNucleotides + 1]int
	Interior2x2Loop [][][][][][]int

	// LogExtrapolationConstant is used to scale energy parameters for hairpin,
	// bulge and interior loops when the length of the loop is greather than
	// `MaxLenLoop`.
	LogExtrapolationConstant                                                            float64
	MultiLoopUnpairedNucleotideBonus, MultiLoopClosingPenalty, TerminalAUPenalty, Ninio int
	// size: [NbDistinguishableBasePairs]int
	MultiLoopIntern []int

	// Some tetraloops particularly stable tetraloops are assigned an energy
	// bonus. For example:
	// ```
	// 	GAAA    -200
	// ```
	// assigns a bonus energy of -2 kcal/mol to TetraLoop containing
	// the sequence GAAA.
	TetraLoop map[string]int
	TriLoop   map[string]int
	HexaLoop  map[string]int

	MaxNinio int
}

// BasePairType is a type to hold information of the type of a base pair.
// The chosen numbers denote where the energy paramater values can be found
// for the base pair type in the `EnergyParams` energy parameter matrices.
type BasePairType int

const (
	// CG occurs when the base C (on the five prime end) binds to the base G
	// (on the three prime end)
	CG BasePairType = 0
	// GC occurs when the base G (on the five prime end) binds to the base C
	// (on the three prime end)
	GC = 1
	// GU occurs when the base G (on the five prime end) binds to the base U
	// (on the three prime end)
	GU = 2
	// UG occurs when the base U (on the five prime end) binds to the base G
	// (on the three prime end)
	UG = 3
	// Au occurs when the base A (on the five prime end) binds to the base U
	// (on the three prime end)
	AU = 4
	// UA occurs when the base U (on the five prime end) binds to the base A
	// (on the three prime end)
	UA = 5
	// NoPair denotes that two bases don't pair
	NoPair = -1
)

var (
	nucleotideAEncodedTypeMap = map[byte]BasePairType{'U': AU}
	nucleotideCEncodedTypeMap = map[byte]BasePairType{'G': CG}
	nucleotideGEncodedTypeMap = map[byte]BasePairType{'C': GC, 'U': GU}
	nucleotideUEncodedTypeMap = map[byte]BasePairType{'A': UA, 'G': UG}

	// BasePairEncodedTypeMap is a map that encodes a base pair to its numerical
	// representation that is used to access the values of the energy parameters
	// in the `EnergyParams` struct.
	//
	// Various loop energy parameters in the `EnergyParams` struct depend on
	// the base pair closing the loop.
	// The numerical representation for the 7 distinguishable types of pairs
	// (-1 = no pair, CG=0, GC=1, GU=2, UG=3, AU=4, UA=5, nonstandard=6)
	// The map is:
	// 	   _   A   C   G   U
	// _ {-1, -1, -1, -1, -1}
	// A {-1, -1, -1, -1,  4}
	// C {-1, -1, -1,  0, -1}
	// G {-1, -1,  1, -1,  2}
	// U {-1,  5, -1,  3, -1}
	// For example, BasePairEncodedTypeMap['A']['U'] is 4.
	// The encoded numerical representation of a base pair carries no meaning in
	// itself except for where to find the relevent energy contributions of the
	// base pair in the energy paramaters matrices.
	// Thus, any change to this map must be reflected in the energy
	// parameters matrices, and vice versa.
	//
	// Note that in Go maps that return `int` return 0 by default when a value
	// isn't found. This would cause bases that don't pair to seem like a CG pair.
	// To rectify this, the exported func `BasePairTypeEncodedInt` returns `-1`
	// when a value in this map isn't found. Please consider using that function
	// to access this map.
	BasePairEncodedTypeMap = map[byte]map[byte]BasePairType{
		'A': nucleotideAEncodedTypeMap,
		'C': nucleotideCEncodedTypeMap,
		'G': nucleotideGEncodedTypeMap,
		'U': nucleotideUEncodedTypeMap,
	}

	// NucleotideEncodedIntMap is a map that encodes a nucleotide to its numerical
	// representation that is used to access the values of the energy parameters
	// in the `EnergyParams` struct.
	//
	// It would be ideal to have this map start from 0, but the
	// `RNAfold parameter file v2.0` file format has a flaw in it. Please read
	// the documentation of `newRawEnergyParams()` to understand why this occurs,
	// and how it could be rectified in the future.
	NucleotideEncodedIntMap map[byte]int = map[byte]int{
		'A': 1,
		'C': 2,
		'G': 3,
		'U': 4,
	}
)

// EncodeSequence encodes a sequence into its numerical representation based on
// `NucleotideEncodedIntMap`.
func EncodeSequence(sequence string) (encodedSequence []int) {
	lenSequence := len(sequence)
	encodedSequence = make([]int, lenSequence)

	for i := 0; i < lenSequence; i++ {
		encodedSequence[i] = NucleotideEncodedIntMap[sequence[i]]
	}

	return encodedSequence
}

// EncodeBasePair returns the type of a base pair encoded as an `int`,
// which is used to access energy paramater values in the `EnergyParams` struct.
// See `basePairTypeEncodedIntMap` for a detailed explanation of the encoding.
func EncodeBasePair(fivePrimeBase, threePrimeBase byte) BasePairType {
	if val, ok := BasePairEncodedTypeMap[fivePrimeBase][threePrimeBase]; ok {
		return val
	} else {
		return NoPair
	}
}

// EnergyParamsSet is used to specify the set of RNA free energy paramaters to
// parse.
type EnergyParamsSet int

const (
	// Langdon2018 specifies the set of RNA energy parameters obtained from
	// Grow and Graft Genetic Programming (GGGP) as published
	// in Langdon et al. 2018, "Evolving Better RNAfold
	// Structure Prediction", EuroGP-2018, M. Castelli,
	// L. Sekanina, M. Zhang Eds., Parma. 4-6 April 2018
	Langdon2018 EnergyParamsSet = iota

	// Andronescu2007 specifies the set of RNA energy parameters obtained from
	// Andronescu M, Condon A, Hoos HH, Mathews DH, Murphy KP. Efficient
	// parameter estimation for RNA secondary structure prediction.
	// Bioinformatics. 2007 Jul 1;23(13):i19-28
	Andronescu2007

	// Turner2004 specifies the set of RNA energy parameters obtained from
	// Mathews DH, Disney MD, Childs JL, Schroeder SJ, Zuker M, Turner DH.
	// Incorporating chemical modification constraints into a dynamic
	// programming algorithm for prediction of RNA secondary structure.
	// Proc Natl Acad Sci U S A. 2004;101(19):7287-7292.
	Turner2004

	// Turner1999 specifies the set of RNA energy parameters obtained from
	// Mathews DH, Sabina J, Zuker M, Turner DH. Expanded sequence dependence
	// of thermodynamic parameters improves prediction of RNA secondary
	// structure. J Mol Biol. 1999 May 21;288(5):911-40.
	Turner1999
)

//go:embed param_files/*
var embeddedEnergyParamsDirectory embed.FS
var energyParamsDirectory = "param_files"

var energyParamFileNames map[EnergyParamsSet]string = map[EnergyParamsSet]string{
	Langdon2018:    "rna_langdon2018.par",
	Andronescu2007: "rna_andronescu2007.par",
	Turner2004:     "rna_turner2004.par",
	Turner1999:     "rna_turner1999.par",
}
