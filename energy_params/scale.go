package energy_params

/******************************************************************************

This file defines structs needed to contain information about RNA free energy
params and functions needed to scale free energy params by a specified
temperature. For more information of parsing the energy params, please see
`parse.go`.

Exported structs:
* EnergyParams - Contains all the energy parameters needed for free energy
	calculations. See the struct's declaration for more information on the types
	of free energy parameters available, and on how to access the energy params.

Exported funcs:
* NewEnergyParams - A wrapper function that parses the specified set of free
	energy params and scales the params by the required temperature.
* BasePairTypeEncodedInt -
* EncodeSequence -

Exported variables:
* NucleotideEncodedIntMap

******************************************************************************/

const (
	// The number of distinguishable base pairs:
	// CG, GC, GU, UG, AU, UA, & non-standard
	NbDistinguishableBasePairs int = 7
	// The number of distinguishable nucleotides: A, C, G, U
	NbDistinguishableNucleotides int = 4
	// The maximum loop length
	MaxLenLoop int = 30
	// ZeroCelsiusInKelvin is 0 deg Celsius in Kelvin
	ZeroCelsiusInKelvin float64 = 273.15

	// The temperature at which the energy parameters have been measured at
	measurementTemperatureInCelsius float64 = 37.0
)

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

var (
	nucleotideAEncodedTypeMap = map[byte]int{'U': 4}
	nucleotideCEncodedTypeMap = map[byte]int{'G': 0}
	nucleotideGEncodedTypeMap = map[byte]int{'C': 1, 'U': 2}
	nucleotideUEncodedTypeMap = map[byte]int{'A': 5, 'G': 3}

	// basePairTypeEncodedIntMap is a map that encodes a base pair to its numerical
	// representation that is used to access the values of the energy parameters
	// in the `EnergyParams` struct.
	//
	// Various loop energy parameters in the `EnergyParams` struct depend in on
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
	// For example, BasePairTypeEncodedIntMap['A']['U'] is 4.
	// The encoded numerical representation of a base pair carries no meaning in
	// itself except for where to find the relevent energy contributions of the
	// base pair in the energy paramaters matrices.
	// Thus, any change to this map must be reflected in the energy
	// parameters matrices, and vice versa.
	//
	// Note that this map is un-exported on purpose. Since maps that return `int`
	// in Go return `0` by default when a value isn't found, this map is wrapped
	// by the exported func `BasePairTypeEncodedInt` which returns `-1` when
	// a value in this map isn't found.
	basePairTypeEncodedIntMap = map[byte]map[byte]int{
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

// BasePairTypeEncodedInt returns the type of a base pair encoded as an `int`,
// which is used to access energy paramater values in the `EnergyParams` struct.
// See `basePairTypeEncodedIntMap` for a detailed explanation of the encoding.
func BasePairTypeEncodedInt(fivePrimeBase, threePrimeBase byte) int {
	if val, ok := basePairTypeEncodedIntMap[fivePrimeBase][threePrimeBase]; ok {
		return val
	} else {
		return -1
	}
}

// scaleByTemperature scales energy paramaters according to the specificed temperatue.
// See `rescaleDg` for more information of how energy values are rescaled.
func (rawEnergyParams rawEnergyParams) scaleByTemperature(temperatureInCelsius float64) *EnergyParams {

	// set the non-matrix energy parameters
	var params *EnergyParams = &EnergyParams{
		LogExtrapolationConstant:         rescaleDgFloat64(rawEnergyParams.logExtrapolationConstant, 0, temperatureInCelsius),
		TerminalAUPenalty:                rescaleDg(rawEnergyParams.terminalAU37C, rawEnergyParams.terminalAUEnthalpy, temperatureInCelsius),
		MultiLoopUnpairedNucleotideBonus: rescaleDg(rawEnergyParams.multiLoopBase37C, rawEnergyParams.multiLoopBaseEnthalpy, temperatureInCelsius),
		MultiLoopClosingPenalty:          rescaleDg(rawEnergyParams.multiLoopClosing37C, rawEnergyParams.multiLoopClosingEnthalpy, temperatureInCelsius),
		Ninio:                            rescaleDg(rawEnergyParams.ninio37C, rawEnergyParams.ninioEnthalpy, temperatureInCelsius),
		MaxNinio:                         rawEnergyParams.maxNinio,
	}

	params.HairpinLoop = rescaleDgMatrix(rawEnergyParams.hairpinLoopEnergy37C, rawEnergyParams.hairpinLoopEnthalpy, temperatureInCelsius, identityInt).([]int)
	params.Bulge = rescaleDgMatrix(rawEnergyParams.bulgeEnergy37C, rawEnergyParams.bulgeEnthalpy, temperatureInCelsius, identityInt).([]int)
	params.InteriorLoop = rescaleDgMatrix(rawEnergyParams.interiorLoopEnergy37C, rawEnergyParams.interiorLoopEnthalpy, temperatureInCelsius, identityInt).([]int)

	params.MultiLoopIntern = make([]int, MaxLenLoop+1)
	for i := 0; i <= MaxLenLoop; i++ {
		params.MultiLoopIntern[i] = rescaleDg(rawEnergyParams.multiLoopIntern37C, rawEnergyParams.multiLoopInternEnthalpy, temperatureInCelsius)
	}

	params.TetraLoop = make(map[string]int)
	for loop := range rawEnergyParams.tetraLoopEnergy37C {
		params.TetraLoop[loop] = rescaleDg(rawEnergyParams.tetraLoopEnergy37C[loop], rawEnergyParams.tetraLoopEnthalpy[loop], temperatureInCelsius)
	}

	params.TriLoop = make(map[string]int)
	for loop := range rawEnergyParams.triLoopEnergy37C {
		params.TriLoop[loop] = rescaleDg(rawEnergyParams.triLoopEnergy37C[loop], rawEnergyParams.triLoopEnthalpy[loop], temperatureInCelsius)
	}

	params.HexaLoop = make(map[string]int)
	for loop := range rawEnergyParams.hexaLoopEnergy37C {
		params.HexaLoop[loop] = rescaleDg(rawEnergyParams.hexaLoopEnergy37C[loop], rawEnergyParams.hexaLoopEnthalpy[loop], temperatureInCelsius)
	}

	params.StackingPair = rescaleDgMatrix(rawEnergyParams.stackingPairEnergy37C, rawEnergyParams.stackingPairEnthalpy, temperatureInCelsius, identityInt).([][]int)

	/* mismatches */
	params.MismatchInteriorLoop = rescaleDgMatrix(rawEnergyParams.mismatchInteriorLoopEnergy37C, rawEnergyParams.mismatchInteriorLoopEnthalpy, temperatureInCelsius, identityInt).([][][]int)
	params.MismatchHairpinLoop = rescaleDgMatrix(rawEnergyParams.mismatchHairpinLoopEnergy37C, rawEnergyParams.mismatchHairpinLoopEnthalpy, temperatureInCelsius, identityInt).([][][]int)
	params.Mismatch1xnInteriorLoop = rescaleDgMatrix(rawEnergyParams.mismatch1xnInteriorLoopEnergy37C, rawEnergyParams.mismatch1xnInteriorLoopEnthalpy, temperatureInCelsius, identityInt).([][][]int)
	params.Mismatch2x3InteriorLoop = rescaleDgMatrix(rawEnergyParams.mismatch2x3InteriorLoopEnergy37C, rawEnergyParams.mismatch2x3InteriorLoopEnthalpy, temperatureInCelsius, identityInt).([][][]int)

	params.MismatchMultiLoop = rescaleDgMatrix(rawEnergyParams.mismatchMultiLoopEnergy37C, rawEnergyParams.mismatchMultiLoopEnthalpy, temperatureInCelsius, onlyLessThanOrEqualToZero).([][][]int)
	params.MismatchExteriorLoop = rescaleDgMatrix(rawEnergyParams.mismatchExteriorLoopEnergy37C, rawEnergyParams.mismatchExteriorLoopEnthalpy, temperatureInCelsius, onlyLessThanOrEqualToZero).([][][]int)

	/* dangling ends energies */
	params.DanglingEndsFivePrime = rescaleDgMatrix(rawEnergyParams.danglingEndsFivePrimeEnergy37C, rawEnergyParams.danglingEndsFivePrimeEnthalpy, temperatureInCelsius, onlyLessThanOrEqualToZero).([][]int)
	params.DanglingEndsThreePrime = rescaleDgMatrix(rawEnergyParams.danglingEndsThreePrimeEnergy37C, rawEnergyParams.danglingEndsThreePrimeEnthalpy, temperatureInCelsius, onlyLessThanOrEqualToZero).([][]int)

	/* interior 1x1 loops */
	params.Interior1x1Loop = rescaleDgMatrix(rawEnergyParams.interior1x1LoopEnergy37C, rawEnergyParams.interior1x1LoopEnthalpy, temperatureInCelsius, identityInt).([][][][]int)

	/* interior 2x1 loops */
	params.Interior2x1Loop = rescaleDgMatrix(rawEnergyParams.interior2x1LoopEnergy37C, rawEnergyParams.interior2x1LoopEnthalpy, temperatureInCelsius, identityInt).([][][][][]int)

	/* interior 2x2 loops */
	params.Interior2x2Loop = rescaleDgMatrix(rawEnergyParams.interior2x2LoopEnergy37C, rawEnergyParams.interior2x2LoopEnthalpy, temperatureInCelsius, identityInt).([][][][][][]int)

	return params
}

// Rescale Gibbs free energy according to the equation dG = dH - T * dS
// where dG is the change in Gibbs free energy
// 			dH is the change in enthalpy
// 			dS is the change in entropy
// 			T is the temperature
// more information: https://chemed.chem.purdue.edu/genchem/topicreview/bp/ch21/gibbs.php
func rescaleDg(dG, dH int, temperatureInCelsius float64) int {
	// if temperate == measurementTemperatureInCelsius then below calculation will
	// always return dG. So we save some computation with this check.
	if temperatureInCelsius == measurementTemperatureInCelsius {
		return dG
	}

	measurementTemperatureInKelvin := measurementTemperatureInCelsius + ZeroCelsiusInKelvin
	temperatureInKelvin := temperatureInCelsius + ZeroCelsiusInKelvin
	var T float64 = float64(temperatureInKelvin / measurementTemperatureInKelvin)

	dGFloat64 := float64(dG)
	dHFloat64 := float64(dH)

	dSFloat64 := dHFloat64 - dGFloat64

	return int(dHFloat64 - dSFloat64*T)
}

// rescaleDgFloat64 is the same as rescaleDg, but for float64
func rescaleDgFloat64(dG, dH, temperatureInCelsius float64) float64 {
	// if temperate == energyParamsTemperature then below calculation will always
	// return dG. So we save some computation with this check.
	if temperatureInCelsius == measurementTemperatureInCelsius {
		return dG
	}

	defaultEnergyParamsTemperatureKelvin := measurementTemperatureInCelsius + ZeroCelsiusInKelvin
	temperatureKelvin := temperatureInCelsius + ZeroCelsiusInKelvin
	var T float64 = temperatureKelvin / defaultEnergyParamsTemperatureKelvin

	dS := dH - dG
	return dH - dS*T
}

// rescaleDgMatrix maps the function `rescaleDg` onto a given `energy` and
// `enthalpy` matrix, and returns a matrix with the scaled energy param values
// (which is the same size as `energy` and `enthalpy`).
//
// Since some energy param values (such as `MismatchMultiLoop` and
// `DanglingEndsThreePrime`) require that energy param values be less than
// or equal to 0, this function includes a `fn` argument which allows the
// caller to specify how to transfrom the computed free energy value (i.e. the
// value after scaling). In most cases `fn` will default to `idInt` (which is
// the identify function for type `int`), but in some cases (as with
// `MismatchMultiLoop`), `onlyLessThanOrEqualToZero` will be used as `fn`'s
// argument.
func rescaleDgMatrix(energy interface{}, enthalpy interface{}, temperatureInCelsius float64, fn func(int) int) (ret interface{}) {
	switch energy := energy.(type) {
	case []int:
		return rescaleDg1Dim(energy, enthalpy.([]int), temperatureInCelsius, fn)
	case [][]int:
		return rescaleDg2Dim(energy, enthalpy.([][]int), temperatureInCelsius, fn)
	case [][][]int:
		return rescaleDg3Dim(energy, enthalpy.([][][]int), temperatureInCelsius, fn)
	case [][][][]int:
		return rescaleDg4Dim(energy, enthalpy.([][][][]int), temperatureInCelsius, fn)
	case [][][][][]int:
		return rescaleDg5Dim(energy, enthalpy.([][][][][]int), temperatureInCelsius, fn)
	case [][][][][][]int:
		return rescaleDg6Dim(energy, enthalpy.([][][][][][]int), temperatureInCelsius, fn)
	}
	return
}

// identityInt is the identify func for int values
func identityInt(x int) int {
	return x
}

// onlyLessThanOrEqualToZero returns x if x <= 0, else 0
func onlyLessThanOrEqualToZero(x int) int {
	return minInt(0, x)
}

func rescaleDg1Dim(energy []int, enthalpy []int, temperatureInCelsius float64, fn func(int) int) (ret []int) {
	lenEnergy := len(energy)
	ret = make([]int, lenEnergy)
	for i := 0; i < lenEnergy; i++ {
		val := rescaleDg(energy[i], enthalpy[i], temperatureInCelsius)
		ret[i] = fn(val)
	}
	return
}

func rescaleDg2Dim(energy [][]int, enthalpy [][]int, temperatureInCelsius float64, fn func(int) int) (ret [][]int) {
	lenEnergy := len(energy)
	ret = make([][]int, lenEnergy)
	for i := 0; i < lenEnergy; i++ {
		ret[i] = rescaleDgMatrix(energy[i], enthalpy[i], temperatureInCelsius, fn).([]int)
	}
	return
}

func rescaleDg3Dim(energy [][][]int, enthalpy [][][]int, temperatureInCelsius float64, fn func(int) int) (ret [][][]int) {
	lenEnergy := len(energy)
	ret = make([][][]int, lenEnergy)
	for i := 0; i < lenEnergy; i++ {
		ret[i] = rescaleDgMatrix(energy[i], enthalpy[i], temperatureInCelsius, fn).([][]int)
	}
	return
}

func rescaleDg4Dim(energy [][][][]int, enthalpy [][][][]int, temperatureInCelsius float64, fn func(int) int) (ret [][][][]int) {
	lenEnergy := len(energy)
	ret = make([][][][]int, lenEnergy)
	for i := 0; i < lenEnergy; i++ {
		ret[i] = rescaleDgMatrix(energy[i], enthalpy[i], temperatureInCelsius, fn).([][][]int)
	}
	return
}

func rescaleDg5Dim(energy [][][][][]int, enthalpy [][][][][]int, temperatureInCelsius float64, fn func(int) int) (ret [][][][][]int) {
	lenEnergy := len(energy)
	ret = make([][][][][]int, lenEnergy)
	for i := 0; i < lenEnergy; i++ {
		ret[i] = rescaleDgMatrix(energy[i], enthalpy[i], temperatureInCelsius, fn).([][][][]int)
	}
	return
}

func rescaleDg6Dim(energy [][][][][][]int, enthalpy [][][][][][]int, temperatureInCelsius float64, fn func(int) int) (ret [][][][][][]int) {
	lenEnergy := len(energy)
	ret = make([][][][][][]int, lenEnergy)
	for i := 0; i < lenEnergy; i++ {
		ret[i] = rescaleDgMatrix(energy[i], enthalpy[i], temperatureInCelsius, fn).([][][][][]int)
	}
	return
}

// minInt returns the minimum of two ints
func minInt(a, b int) int {
	if a < b {
		return a
	}
	return b
}
