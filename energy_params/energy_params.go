package energy_params

import (
	"bufio"
	"embed"
	"fmt"
	"log"
	"strconv"
	"strings"
)

/******************************************************************************

This file defines functions and structs needed to parse RNA energy parameters
as specified by the `RNAfold parameter file v2.0` file format from
[ViennaRNA](https://github.com/ViennaRNA/ViennaRNA#energy-parameters).

Exported structs:
* EnergyParams - Contains all the energy parameters needed for free energy
	calculations.



******************************************************************************/

/**
* These energy parameters are nearest neighbor parameters for RNA folding
* compiled by the Turner group. The names of variables from the original repo
* (ViennaRNA) have been changed to be more descriptive and the original name
* of the variable has been commented out in the line above the variable
* declaration (to allow for easy comparison and updating of the variables).
* `mfe.go` includes an explanation of these variables in the declaration of the
* `energyParams` struct.
 */

const (
	// The number of distinguishable base pairs:
	// CG, GC, GU, UG, AU, UA, & non-standard (see `energy_params.md`)
	NbDistinguishableBasePairs int = 7
	// The number of distinguishable nucleotides: A, C, G, U
	NbDistinguishableNucleotides int = 4
	// The maximum loop length
	MaxLenLoop int = 30

	// ZeroCKelvin is 0 deg Celsius in Kelvin
	ZeroCKelvin float64 = 273.15

	// The constant used to extrapolate energy values when length of loop > `MaxLenLoop`
	defaultLogExtrapolationConstantAt37C float64 = 107.856

	// The temperature at which the energy parameters have been measured at
	measurementTemperatureInCelsius float64 = 37.0

	// inf is infinity as used in minimization routines (INT_MAX/10)
	inf int = 10000000
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

//go:embed param_files/*
var embeddedEnergyParamsDirectory embed.FS
var energyParamsDirectory = "param_files"

var energyParamFileNames map[EnergyParamsSet]string = map[EnergyParamsSet]string{
	Langdon2018:    "rna_langdon2018.par",
	Andronescu2007: "rna_andronescu2007.par",
	Turner2004:     "rna_turner2004.par",
	Turner1999:     "rna_turner1999.par",
}

// rawEnergyParams is an un-exported intermediate struct used to store parsed
// energy param values. It is converted into usable energy params by the
// function `scaleByTemperature()`. For more information on the fields of
// this struct, read the documentation of the fields of the `EnergyParams`
// struct
type rawEnergyParams struct {
	interior2x2LoopEnergy37C [][][][][][]int
	interior2x2LoopEnthalpy  [][][][][][]int

	interior2x1LoopEnergy37C [][][][][]int
	interior2x1LoopEnthalpy  [][][][][]int

	interior1x1LoopEnergy37C [][][][]int
	interior1x1LoopEnthalpy  [][][][]int

	mismatchInteriorLoopEnergy37C    [][][]int
	mismatchInteriorLoopEnthalpy     [][][]int
	mismatchHairpinLoopEnergy37C     [][][]int
	mismatchHairpinLoopEnthalpy      [][][]int
	mismatchMultiLoopEnergy37C       [][][]int
	mismatchMultiLoopEnthalpy        [][][]int
	mismatch1xnInteriorLoopEnergy37C [][][]int
	mismatch1xnInteriorLoopEnthalpy  [][][]int
	mismatch2x3InteriorLoopEnergy37C [][][]int
	mismatch2x3InteriorLoopEnthalpy  [][][]int
	mismatchExteriorLoopEnergy37C    [][][]int
	mismatchExteriorLoopEnthalpy     [][][]int

	dangle5Energy37C      [][]int
	dangle5Enthalpy       [][]int
	dangle3Energy37C      [][]int
	dangle3Enthalpy       [][]int
	stackingPairEnergy37C [][]int
	stackingPairEnthalpy  [][]int

	hairpinLoopEnergy37C  []int
	hairpinLoopEnthalpy   []int
	bulgeEnergy37C        []int
	bulgeEnthalpy         []int
	interiorLoopEnergy37C []int
	interiorLoopEnthalpy  []int

	multiLoopIntern37C       int
	multiLoopInternEnthalpy  int
	multiLoopClosing37C      int
	multiLoopClosingEnthalpy int
	multiLoopBase37C         int
	multiLoopBaseEnthalpy    int
	maxNinio                 int
	ninio37C                 int
	ninioEnthalpy            int
	terminalAU37C            int
	terminalAUEnthalpy       int
	logExtrapolationConstant float64

	tetraLoops         map[string]bool
	tetraLoopEnergy37C map[string]int
	tetraLoopEnthalpy  map[string]int

	triLoops         map[string]bool
	triLoopEnergy37C map[string]int
	triLoopEnthalpy  map[string]int

	hexaLoops         map[string]bool
	hexaLoopEnergy37C map[string]int
	hexaLoopEnthalpy  map[string]int
}

// Contains all the energy parameters needed for the free energy calculations.
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
// 						interior2x2Loop[2][6][1][4][2][3]
// ```
// (See `basePairEncodedTypeMap()` and `encodeSequence()` for more details on how
// pairs and unpaired nucleotides are encoded)
// Note that this sequence is symmetric so the sequence is equivalent to:
// ```
// 					5'-UCGC-3'
// 					3'-AUAG-5'
// ```
// which means the energy of the sequence is equivalent to:
// ```
// 	pairs:                    UA GC C  G  A  U
// 						interior2x2Loop[2][6][1][4][2][3]
// ```
type EnergyParams struct {

	// The matrix of free energies for stacked pairs, indexed by the two encoded closing
	// pairs. The list should be formatted as symmetric a `7*7` matrix, conforming
	// to the order explained above. As an example the stacked pair
	// ```
	// 					5'-GU-3'
	// 					3'-CA-5'
	// ```
	// corresponds to the entry StackingPair[2][5] (GC=2, AU=5) which should be
	// identical to StackingPair[5][2] (AU=5, GC=2).
	// size: [NBDistinguishablePairs + 1][NBDistinguishablePairs + 1]int
	StackingPair [][]int
	// Free energies of hairpin loops as a function of size. The list should
	// contain 31 (MaxLenLoop + 1) entries. Since the minimum size of a hairpin loop
	// is 3 and we start counting with 0, the first three values should be `inf` to
	// indicate a forbidden value.
	// size: [MaxLenLoop + 1]int
	HairpinLoop []int
	// Free energies of Bulge loops. Should contain 31 (MaxLenLoop + 1) entries, the
	// first one being `inf`.
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
	// unpaired bases. Since we distinguish 5 bases the list contains
	// `7*5*5` entries. The order is such that for example the mismatch
	//
	// ```
	// 							5'-CU-3'
	// 							3'-GC-5'
	// ```
	// corresponds to entry MismatchInteriorLoop[1][4][2] (CG=1, U=4, C=2).
	//
	// More information about mismatch energy: https://rna.urmc.rochester.edu/NNDB/turner04/tm.html
	// size: [NBDistinguishablePairs + 1][NBDistinguishableNucleotides + 1][NBDistinguishableNucleotides + 1]int
	MismatchInteriorLoop [][][]int
	// size: [NBDistinguishablePairs + 1][NBDistinguishableNucleotides + 1][NBDistinguishableNucleotides + 1]int
	Mismatch1xnInteriorLoop [][][]int
	// size: [NBDistinguishablePairs + 1][NBDistinguishableNucleotides + 1][NBDistinguishableNucleotides + 1]int
	Mismatch2x3InteriorLoop [][][]int
	// size: [NBDistinguishablePairs + 1][NBDistinguishableNucleotides + 1][NBDistinguishableNucleotides + 1]int
	MismatchExteriorLoop [][][]int
	// size: [NBDistinguishablePairs + 1][NBDistinguishableNucleotides + 1][NBDistinguishableNucleotides + 1]int
	MismatchHairpinLoop [][][]int
	// size: [NBDistinguishablePairs + 1][NBDistinguishableNucleotides + 1][NBDistinguishableNucleotides + 1]int
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
	// corresponds to entry DanglingEndsFivePrime[1][3] (CG=1, G=3).
	//
	// More information about dangling ends: https://rna.urmc.rochester.edu/NNDB/turner04/de.html
	// size: [NBDistinguishablePairs + 1][NBDistinguishableNucleotides + 1]int
	DanglingEndsFivePrime [][]int

	// Same as above for bases on the 3' side of a helix.
	// ```
	// 			       5'- A-3'
	// 			       3'-AU-5'
	// ```
	// corresponds to entry DanglingEndsThreePrime[5][1] (AU=5, A=1).
	// size: [NBDistinguishablePairs + 1][NBDistinguishableNucleotides + 1]int
	DanglingEndsThreePrime [][]int

	// Free energies for symmetric size 2 interior loops. `7*7*5*5` entries.
	// Example:
	// ```
	// 						5'-CUU-3'
	// 						3'-GCA-5'
	// ```
	// corresponds to entry Interior1x1Loop[1][5][4][2] (CG=1, AU=5, U=4, C=2),
	// which should be identical to Interior1x1Loop[5][1][2][4] (AU=5, CG=1, C=2, U=4).
	// size: [NBDistinguishablePairs + 1][NBDistinguishablePairs + 1][NBDistinguishableNucleotides + 1][NBDistinguishableNucleotides + 1]int
	Interior1x1Loop [][][][]int

	// Free energies for 2x1 interior loops, where 2 is the number of unpaired
	// nucelobases on the larger 'side' of the interior loop and 1 is the number of
	// unpaired nucelobases on the smaller 'side' of the interior loop.
	// `7*7*5*5*5` entries.
	// Example:
	// ```
	// 						5'-CUUU-3'
	// 						3'-GC A-5'
	// ```
	// corresponds to entry Interior2x1Loop[1][5][4][4][2] (CG=1, AU=5, U=4, U=4, C=2).
	// Note that this matrix is always accessed in the 5' to 3' direction with the
	// larger number of unpaired nucleobases first.
	// size: [NBDistinguishablePairs + 1][NBDistinguishablePairs + 1][NBDistinguishableNucleotides + 1][NBDistinguishableNucleotides + 1][NBDistinguishableNucleotides + 1]int
	Interior2x1Loop [][][][][]int

	// Free energies for symmetric size 4 interior loops. `7*7*5*5*5*5` entries.
	// Example:
	// ```
	// 						5'-CUAU-3'
	// 						3'-GCCA-5'
	// ```
	// corresponds to entry Interior2x2Loop[1][5][4][1][2][2] (CG=1, AU=5, U=4, A=1, C=2, C=2),
	// which should be identical to Interior2x2Loop[5][1][2][2][1][4] (AU=5, CG=1, C=2, C=2, A=1, U=4).
	// size: [NBDistinguishablePairs + 1][NBDistinguishablePairs + 1][NBDistinguishableNucleotides + 1][NBDistinguishableNucleotides + 1][NBDistinguishableNucleotides + 1][NBDistinguishableNucleotides + 1]int
	Interior2x2Loop [][][][][][]int

	// Parameter for logarithmic loop energy extrapolation. Used to scale energy
	// parameters for hairpin, bulge and interior loops when the length of the loop
	// is greather than `MaxLenLoop`.
	LogExtrapolationConstant                                                            float64
	MultiLoopUnpairedNucleotideBonus, MultiLoopClosingPenalty, TerminalAUPenalty, Ninio int
	// size: [NBDistinguishablePairs + 1]int
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

func readLine(scanner *bufio.Scanner) (string, bool) {
	lineAvailable := scanner.Scan()

	if err := scanner.Err(); err != nil {
		log.Fatal(err)
	}

	if lineAvailable {
		line := scanner.Text()
		return line, lineAvailable
	} else {
		return "", false
	}
}

func NewEnergyParams(energyParamsSet EnergyParamsSet, temperatureInCelsius float64) *EnergyParams {
	return newRawEnergyParams(energyParamsSet).scaleByTemperature(temperatureInCelsius)
}

func newRawEnergyParams(energyParamsSet EnergyParamsSet) (rawEnergyParams rawEnergyParams) {
	fileName := energyParamFileNames[energyParamsSet]
	file, err := embeddedEnergyParamsDirectory.Open(energyParamsDirectory + "/" + fileName)
	if err != nil {
		log.Fatal(err)
	}
	defer file.Close()

	scanner := bufio.NewScanner(file)

	headerLine, lineAvailable := readLine(scanner)
	if lineAvailable && headerLine != "## RNAfold parameter file v2.0" {
		panic("Missing header line in file.\nMay be this file has not v2.0 format.\n")
	}

	line, lineAvailable := readLine(scanner)
	for lineAvailable {
		var energyParamSet string
		n, _ := fmt.Sscanf(line, "# %255s", &energyParamSet)

		if n == 1 {
			switch energyParamSet {
			case "END":
				break
			case "stack":
				rawEnergyParams.stackingPairEnergy37C = parseParamValues(scanner, NbDistinguishableBasePairs, NbDistinguishableBasePairs).([][]int)
				rawEnergyParams.stackingPairEnergy37C = addOffset(rawEnergyParams.stackingPairEnergy37C, preOffset, 1, 1).([][]int)

			case "stack_enthalpies":
				rawEnergyParams.stackingPairEnthalpy = parseParamValues(scanner, NbDistinguishableBasePairs, NbDistinguishableBasePairs).([][]int)
				rawEnergyParams.stackingPairEnthalpy = addOffset(rawEnergyParams.stackingPairEnthalpy, preOffset, 1, 1).([][]int)

			case "hairpin":
				rawEnergyParams.hairpinLoopEnergy37C = parseParamValues(scanner, 31).([]int)

			case "hairpin_enthalpies":
				rawEnergyParams.hairpinLoopEnthalpy = parseParamValues(scanner, 31).([]int)

			case "bulge":
				rawEnergyParams.bulgeEnergy37C = parseParamValues(scanner, 31).([]int)

			case "bulge_enthalpies":
				rawEnergyParams.bulgeEnthalpy = parseParamValues(scanner, 31).([]int)

			case "interior":
				rawEnergyParams.interiorLoopEnergy37C = parseParamValues(scanner, 31).([]int)

			case "interior_enthalpies":
				rawEnergyParams.interiorLoopEnthalpy = parseParamValues(scanner, 31).([]int)

			case "mismatch_exterior":
				rawEnergyParams.mismatchExteriorLoopEnergy37C = parseParamValues(scanner, NbDistinguishableBasePairs, NbDistinguishableNucleotides+1, NbDistinguishableNucleotides+1).([][][]int)
				rawEnergyParams.mismatchExteriorLoopEnergy37C = addOffset(rawEnergyParams.mismatchExteriorLoopEnergy37C, preOffset, 1, 0, 0).([][][]int)

			case "mismatch_exterior_enthalpies":
				rawEnergyParams.mismatchExteriorLoopEnthalpy = parseParamValues(scanner, NbDistinguishableBasePairs, NbDistinguishableNucleotides+1, NbDistinguishableNucleotides+1).([][][]int)
				rawEnergyParams.mismatchExteriorLoopEnthalpy = addOffset(rawEnergyParams.mismatchExteriorLoopEnthalpy, preOffset, 1, 0, 0).([][][]int)

			case "mismatch_hairpin":
				rawEnergyParams.mismatchHairpinLoopEnergy37C = parseParamValues(scanner, NbDistinguishableBasePairs, NbDistinguishableNucleotides+1, NbDistinguishableNucleotides+1).([][][]int)
				rawEnergyParams.mismatchHairpinLoopEnergy37C = addOffset(rawEnergyParams.mismatchHairpinLoopEnergy37C, preOffset, 1, 0, 0).([][][]int)

			case "mismatch_hairpin_enthalpies":
				rawEnergyParams.mismatchHairpinLoopEnthalpy = parseParamValues(scanner, NbDistinguishableBasePairs, NbDistinguishableNucleotides+1, NbDistinguishableNucleotides+1).([][][]int)
				rawEnergyParams.mismatchHairpinLoopEnthalpy = addOffset(rawEnergyParams.mismatchHairpinLoopEnthalpy, preOffset, 1, 0, 0).([][][]int)

			case "mismatch_interior":
				rawEnergyParams.mismatchInteriorLoopEnergy37C = parseParamValues(scanner, NbDistinguishableBasePairs, NbDistinguishableNucleotides+1, NbDistinguishableNucleotides+1).([][][]int)
				rawEnergyParams.mismatchInteriorLoopEnergy37C = addOffset(rawEnergyParams.mismatchInteriorLoopEnergy37C, preOffset, 1, 0, 0).([][][]int)

			case "mismatch_interior_enthalpies":
				rawEnergyParams.mismatchInteriorLoopEnthalpy = parseParamValues(scanner, NbDistinguishableBasePairs, NbDistinguishableNucleotides+1, NbDistinguishableNucleotides+1).([][][]int)
				rawEnergyParams.mismatchInteriorLoopEnthalpy = addOffset(rawEnergyParams.mismatchInteriorLoopEnthalpy, preOffset, 1, 0, 0).([][][]int)

			case "mismatch_interior_1n":
				rawEnergyParams.mismatch1xnInteriorLoopEnergy37C = parseParamValues(scanner, NbDistinguishableBasePairs, NbDistinguishableNucleotides+1, NbDistinguishableNucleotides+1).([][][]int)
				rawEnergyParams.mismatch1xnInteriorLoopEnergy37C = addOffset(rawEnergyParams.mismatch1xnInteriorLoopEnergy37C, preOffset, 1, 0, 0).([][][]int)

			case "mismatch_interior_1n_enthalpies":
				rawEnergyParams.mismatch1xnInteriorLoopEnthalpy = parseParamValues(scanner, NbDistinguishableBasePairs, NbDistinguishableNucleotides+1, NbDistinguishableNucleotides+1).([][][]int)
				rawEnergyParams.mismatch1xnInteriorLoopEnthalpy = addOffset(rawEnergyParams.mismatch1xnInteriorLoopEnthalpy, preOffset, 1, 0, 0).([][][]int)

			case "mismatch_interior_23":
				rawEnergyParams.mismatch2x3InteriorLoopEnergy37C = parseParamValues(scanner, NbDistinguishableBasePairs, NbDistinguishableNucleotides+1, NbDistinguishableNucleotides+1).([][][]int)
				rawEnergyParams.mismatch2x3InteriorLoopEnergy37C = addOffset(rawEnergyParams.mismatch2x3InteriorLoopEnergy37C, preOffset, 1, 0, 0).([][][]int)

			case "mismatch_interior_23_enthalpies":
				rawEnergyParams.mismatch2x3InteriorLoopEnthalpy = parseParamValues(scanner, NbDistinguishableBasePairs, NbDistinguishableNucleotides+1, NbDistinguishableNucleotides+1).([][][]int)
				rawEnergyParams.mismatch2x3InteriorLoopEnthalpy = addOffset(rawEnergyParams.mismatch2x3InteriorLoopEnthalpy, preOffset, 1, 0, 0).([][][]int)

			case "mismatch_multi":
				rawEnergyParams.mismatchMultiLoopEnergy37C = parseParamValues(scanner, NbDistinguishableBasePairs, NbDistinguishableNucleotides+1, NbDistinguishableNucleotides+1).([][][]int)
				rawEnergyParams.mismatchMultiLoopEnergy37C = addOffset(rawEnergyParams.mismatchMultiLoopEnergy37C, preOffset, 1, 0, 0).([][][]int)

			case "mismatch_multi_enthalpies":
				rawEnergyParams.mismatchMultiLoopEnthalpy = parseParamValues(scanner, NbDistinguishableBasePairs, NbDistinguishableNucleotides+1, NbDistinguishableNucleotides+1).([][][]int)
				rawEnergyParams.mismatchMultiLoopEnthalpy = addOffset(rawEnergyParams.mismatchMultiLoopEnthalpy, preOffset, 1, 0, 0).([][][]int)

			case "int11":
				rawEnergyParams.interior1x1LoopEnergy37C = parseParamValues(scanner, NbDistinguishableBasePairs, NbDistinguishableBasePairs, NbDistinguishableNucleotides+1, NbDistinguishableNucleotides+1).([][][][]int)
				rawEnergyParams.interior1x1LoopEnergy37C = addOffset(rawEnergyParams.interior1x1LoopEnergy37C, preOffset, 1, 1, 0, 0).([][][][]int)

			case "int11_enthalpies":
				rawEnergyParams.interior1x1LoopEnthalpy = parseParamValues(scanner, NbDistinguishableBasePairs, NbDistinguishableBasePairs, NbDistinguishableNucleotides+1, NbDistinguishableNucleotides+1).([][][][]int)
				rawEnergyParams.interior1x1LoopEnthalpy = addOffset(rawEnergyParams.interior1x1LoopEnthalpy, preOffset, 1, 1, 0, 0).([][][][]int)

			case "int21":
				rawEnergyParams.interior2x1LoopEnergy37C = parseParamValues(scanner, NbDistinguishableBasePairs, NbDistinguishableBasePairs, NbDistinguishableNucleotides+1, NbDistinguishableNucleotides+1, NbDistinguishableNucleotides+1).([][][][][]int)
				rawEnergyParams.interior2x1LoopEnergy37C = addOffset(rawEnergyParams.interior2x1LoopEnergy37C, preOffset, 1, 1, 0, 0, 0).([][][][][]int)

			case "int21_enthalpies":
				rawEnergyParams.interior2x1LoopEnthalpy = parseParamValues(scanner, NbDistinguishableBasePairs, NbDistinguishableBasePairs, NbDistinguishableNucleotides+1, NbDistinguishableNucleotides+1, NbDistinguishableNucleotides+1).([][][][][]int)
				rawEnergyParams.interior2x1LoopEnthalpy = addOffset(rawEnergyParams.interior2x1LoopEnthalpy, preOffset, 1, 1, 0, 0, 0).([][][][][]int)

			case "int22":
				rawEnergyParams.interior2x2LoopEnergy37C = parseParamValues(scanner, NbDistinguishableBasePairs-1, NbDistinguishableBasePairs-1, NbDistinguishableNucleotides, NbDistinguishableNucleotides, NbDistinguishableNucleotides, NbDistinguishableNucleotides).([][][][][][]int)
				rawEnergyParams.interior2x2LoopEnergy37C = addOffset(rawEnergyParams.interior2x2LoopEnergy37C, preOffset, 1, 1, 1, 1, 1, 1).([][][][][][]int)
				rawEnergyParams.interior2x2LoopEnergy37C = addOffset(rawEnergyParams.interior2x2LoopEnergy37C, postOffset, 1, 1, 0, 0, 0, 0).([][][][][][]int)

			case "int22_enthalpies":
				rawEnergyParams.interior2x2LoopEnthalpy = parseParamValues(scanner, NbDistinguishableBasePairs-1, NbDistinguishableBasePairs-1, NbDistinguishableNucleotides, NbDistinguishableNucleotides, NbDistinguishableNucleotides, NbDistinguishableNucleotides).([][][][][][]int)
				rawEnergyParams.interior2x2LoopEnthalpy = addOffset(rawEnergyParams.interior2x2LoopEnthalpy, preOffset, 1, 1, 1, 1, 1, 1).([][][][][][]int)
				rawEnergyParams.interior2x2LoopEnthalpy = addOffset(rawEnergyParams.interior2x2LoopEnthalpy, postOffset, 1, 1, 0, 0, 0, 0).([][][][][][]int)

			case "dangle5":
				rawEnergyParams.dangle5Energy37C = parseParamValues(scanner, NbDistinguishableBasePairs, NbDistinguishableNucleotides+1).([][]int)
				rawEnergyParams.dangle5Energy37C = addOffset(rawEnergyParams.dangle5Energy37C, preOffset, 1, 0).([][]int)

			case "dangle5_enthalpies":
				rawEnergyParams.dangle5Enthalpy = parseParamValues(scanner, NbDistinguishableBasePairs, NbDistinguishableNucleotides+1).([][]int)
				rawEnergyParams.dangle5Enthalpy = addOffset(rawEnergyParams.dangle5Enthalpy, preOffset, 1, 0).([][]int)

			case "dangle3":
				rawEnergyParams.dangle3Energy37C = parseParamValues(scanner, NbDistinguishableBasePairs, NbDistinguishableNucleotides+1).([][]int)
				rawEnergyParams.dangle3Energy37C = addOffset(rawEnergyParams.dangle3Energy37C, preOffset, 1, 0).([][]int)

			case "dangle3_enthalpies":
				rawEnergyParams.dangle3Enthalpy = parseParamValues(scanner, NbDistinguishableBasePairs, NbDistinguishableNucleotides+1).([][]int)
				rawEnergyParams.dangle3Enthalpy = addOffset(rawEnergyParams.dangle3Enthalpy, preOffset, 1, 0).([][]int)

			case "ML_params":
				mlParams := parseParamValues(scanner, 6).([]int)
				rawEnergyParams.setMultiLoopParams(mlParams)

			case "NINIO":
				ninioParams := parseParamValues(scanner, 3).([]int)
				rawEnergyParams.setNinioParams(ninioParams)

			case "Triloops":
				rawEnergyParams.triLoopEnergy37C,
					rawEnergyParams.triLoopEnthalpy = parseTriTetraHexaLoopParams(scanner)

			case "Tetraloops":
				rawEnergyParams.tetraLoopEnergy37C,
					rawEnergyParams.tetraLoopEnthalpy = parseTriTetraHexaLoopParams(scanner)

			case "Hexaloops":
				rawEnergyParams.hexaLoopEnergy37C,
					rawEnergyParams.hexaLoopEnthalpy = parseTriTetraHexaLoopParams(scanner)

			case "Misc":
				paramLine := parseLineIntoSlice(scanner)
				miscParams := convertToFloat64(paramLine)
				rawEnergyParams.setMiscParams(miscParams)
			}
		}

		line, lineAvailable = readLine(scanner)
	}

	return
}

func parseTriTetraHexaLoopParams(scanner *bufio.Scanner) (energies map[string]int, enthalpies map[string]int) {
	energies, enthalpies = make(map[string]int), make(map[string]int)
	line, lineAvailable := readLine(scanner)
	for lineAvailable {
		if len(strings.TrimSpace(line)) == 0 {
			// blank line encountered
			return
		} else {
			line = removeComments(line)

			if len(strings.TrimSpace(line)) != 0 {
				values := strings.Fields(line)
				if len(values) != 3 {
					panic(fmt.Sprintf("encountered incorrect number of values. expected 3, got %v", len(values)))
				}
				loop, energy, enthalpy := values[0], parseInt(values[1]), parseInt(values[2])
				energies[loop] = energy
				enthalpies[loop] = enthalpy
			} // else line only contain comments so continue reading

			line, lineAvailable = readLine(scanner)
		}
	}
	return
}

// returns a slice of length `length` with all values set to `value`
func newIntSlice(value, length int) (ret []int) {
	ret = make([]int, length)
	for i := 0; i < length; i++ {
		ret[i] = value
	}
	return
}

// adds `length` `inf`s to the front of a slice
func prependInfsToSlice(slice []int, length int) []int {
	return append(newIntSlice(inf, length), slice...)
}

func appendInfsToSlice(slice []int, length int) []int {
	return append(slice, newIntSlice(inf, length)...)
}

type offsetType int

const (
	preOffset offsetType = iota
	postOffset
)

func addOffset(values interface{}, offsetType offsetType, dims ...int) interface{} {
	switch values := values.(type) {
	case []int:
		return addOffset1Dim(values, dims[0], offsetType)
	case [][]int:
		return addOffset2Dim(values, dims[0], dims[1], offsetType)
	case [][][]int:
		return addOffset3Dim(values, dims[0], dims[1], dims[2], offsetType)
	case [][][][]int:
		return addOffset4Dim(values, dims[0], dims[1], dims[2], dims[3], offsetType)
	case [][][][][]int:
		return addOffset5Dim(values, dims[0], dims[1], dims[2], dims[3], dims[4], offsetType)
	case [][][][][][]int:
		return addOffset6Dim(values, dims[0], dims[1], dims[2], dims[3], dims[4], dims[5], offsetType)
	}
	return nil
}

func addOffset1Dim(values []int, dim1Offset int, offsetType offsetType) (ret []int) {
	switch offsetType {
	case preOffset:
		ret = prependInfsToSlice(values, dim1Offset)
	case postOffset:
		ret = appendInfsToSlice(values, dim1Offset)
	}
	return
}

func addOffset2Dim(values [][]int, dim1Offset, dim2Offset int, offsetType offsetType) (ret [][]int) {
	valuesDims := getSliceDims(values)
	currLenDim1, currLenDim2 := valuesDims[0], valuesDims[1]

	newLenDim1 := currLenDim1 + dim1Offset
	ret = make([][]int, 0, newLenDim1)

	var newLenDim2 int = dim2Offset + currLenDim2

	switch offsetType {
	case preOffset:
		for i := 0; i < dim1Offset; i++ {
			infSlice := newIntSlice(inf, newLenDim2)
			ret = append(ret, infSlice)
		}

		for i := 0; i < currLenDim1; i++ {
			ret = append(ret, prependInfsToSlice(values[i], dim2Offset))
		}
	case postOffset:
		for i := 0; i < currLenDim1; i++ {
			ret = append(ret, appendInfsToSlice(values[i], dim2Offset))
		}

		for i := 0; i < dim1Offset; i++ {
			infSlice := newIntSlice(inf, newLenDim2)
			ret = append(ret, infSlice)
		}
	}
	return
}

func addOffset3Dim(values [][][]int, dim1Offset, dim2Offset, dim3Offset int, offsetType offsetType) (ret [][][]int) {
	valuesDims := getSliceDims(values)
	currLenDim1, currLenDim2, currLenDim3 := valuesDims[0], valuesDims[1], valuesDims[2]

	newLenDim1 := currLenDim1 + dim1Offset
	ret = make([][][]int, 0, newLenDim1)

	var newLenDim2 int = dim2Offset + currLenDim2
	var newLenDim3 int = dim3Offset + currLenDim3

	switch offsetType {
	case preOffset:
		for i := 0; i < dim1Offset; i++ {
			var infMatrix [][]int
			infMatrix = addOffset2Dim(infMatrix, newLenDim2, newLenDim3, offsetType)
			ret = append(ret, infMatrix)
		}

		for i := 0; i < currLenDim1; i++ {
			ret = append(ret, addOffset2Dim(values[i], dim2Offset, dim3Offset, offsetType))
		}
	case postOffset:
		for i := 0; i < currLenDim1; i++ {
			ret = append(ret, addOffset2Dim(values[i], dim2Offset, dim3Offset, offsetType))
		}

		for i := 0; i < dim1Offset; i++ {
			var infMatrix [][]int
			infMatrix = addOffset2Dim(infMatrix, newLenDim2, newLenDim3, offsetType)
			ret = append(ret, infMatrix)
		}
	}
	return
}

func addOffset4Dim(values [][][][]int, dim1Offset, dim2Offset, dim3Offset, dim4Offset int, offsetType offsetType) (ret [][][][]int) {
	valuesDims := getSliceDims(values)
	currLenDim1, currLenDim2, currLenDim3, currLenDim4 := valuesDims[0], valuesDims[1], valuesDims[2], valuesDims[3]

	newLenDim1 := currLenDim1 + dim1Offset
	ret = make([][][][]int, 0, newLenDim1)

	var newLenDim2 int = dim2Offset + currLenDim2
	var newLenDim3 int = dim3Offset + currLenDim3
	var newLenDim4 int = dim4Offset + currLenDim4

	switch offsetType {
	case preOffset:
		for i := 0; i < dim1Offset; i++ {
			var infMatrix [][][]int
			infMatrix = addOffset3Dim(infMatrix, newLenDim2, newLenDim3, newLenDim4, offsetType)
			ret = append(ret, infMatrix)
		}

		for i := 0; i < currLenDim1; i++ {
			ret = append(ret, addOffset3Dim(values[i], dim2Offset, dim3Offset, dim4Offset, offsetType))
		}
	case postOffset:
		for i := 0; i < currLenDim1; i++ {
			ret = append(ret, addOffset3Dim(values[i], dim2Offset, dim3Offset, dim4Offset, offsetType))
		}

		for i := 0; i < dim1Offset; i++ {
			var infMatrix [][][]int
			infMatrix = addOffset3Dim(infMatrix, newLenDim2, newLenDim3, newLenDim4, offsetType)
			ret = append(ret, infMatrix)
		}
	}

	return
}

func addOffset5Dim(values [][][][][]int, dim1Offset, dim2Offset, dim3Offset, dim4Offset, dim5Offset int, offsetType offsetType) (ret [][][][][]int) {
	valuesDims := getSliceDims(values)
	currLenDim1, currLenDim2, currLenDim3, currLenDim4, currLenDim5 := valuesDims[0], valuesDims[1], valuesDims[2], valuesDims[3], valuesDims[4]

	newLenDim1 := currLenDim1 + dim1Offset
	ret = make([][][][][]int, 0, newLenDim1)

	var newLenDim2 int = dim2Offset + currLenDim2
	var newLenDim3 int = dim3Offset + currLenDim3
	var newLenDim4 int = dim4Offset + currLenDim4
	var newLenDim5 int = dim5Offset + currLenDim5

	switch offsetType {
	case preOffset:
		for i := 0; i < dim1Offset; i++ {
			var infMatrix [][][][]int
			infMatrix = addOffset4Dim(infMatrix, newLenDim2, newLenDim3, newLenDim4, newLenDim5, offsetType)
			ret = append(ret, infMatrix)
		}

		for i := 0; i < currLenDim1; i++ {
			ret = append(ret, addOffset4Dim(values[i], dim2Offset, dim3Offset, dim4Offset, dim5Offset, offsetType))
		}
	case postOffset:
		for i := 0; i < currLenDim1; i++ {
			ret = append(ret, addOffset4Dim(values[i], dim2Offset, dim3Offset, dim4Offset, dim5Offset, offsetType))
		}

		for i := 0; i < dim1Offset; i++ {
			var infMatrix [][][][]int
			infMatrix = addOffset4Dim(infMatrix, newLenDim2, newLenDim3, newLenDim4, newLenDim5, offsetType)
			ret = append(ret, infMatrix)
		}
	}
	return
}

func addOffset6Dim(values [][][][][][]int, dim1Offset, dim2Offset, dim3Offset, dim4Offset, dim5Offset, dim6Offset int, offsetType offsetType) (ret [][][][][][]int) {
	valuesDims := getSliceDims(values)
	currLenDim1, currLenDim2, currLenDim3, currLenDim4, currLenDim5, currLenDim6 := valuesDims[0], valuesDims[1], valuesDims[2], valuesDims[3], valuesDims[4], valuesDims[5]

	newLenDim1 := currLenDim1 + dim1Offset
	ret = make([][][][][][]int, 0, newLenDim1)

	var newLenDim2 int = dim2Offset + currLenDim2
	var newLenDim3 int = dim3Offset + currLenDim3
	var newLenDim4 int = dim4Offset + currLenDim4
	var newLenDim5 int = dim5Offset + currLenDim5
	var newLenDim6 int = dim6Offset + currLenDim6

	switch offsetType {
	case preOffset:
		for i := 0; i < dim1Offset; i++ {
			var infMatrix [][][][][]int
			infMatrix = addOffset5Dim(infMatrix, newLenDim2, newLenDim3, newLenDim4, newLenDim5, newLenDim6, offsetType)
			ret = append(ret, infMatrix)
		}

		for i := 0; i < currLenDim1; i++ {
			ret = append(ret, addOffset5Dim(values[i], dim2Offset, dim3Offset, dim4Offset, dim5Offset, dim6Offset, offsetType))
		}
	case postOffset:
		for i := 0; i < currLenDim1; i++ {
			ret = append(ret, addOffset5Dim(values[i], dim2Offset, dim3Offset, dim4Offset, dim5Offset, dim6Offset, offsetType))
		}

		for i := 0; i < dim1Offset; i++ {
			var infMatrix [][][][][]int
			infMatrix = addOffset5Dim(infMatrix, newLenDim2, newLenDim3, newLenDim4, newLenDim5, newLenDim6, offsetType)
			ret = append(ret, infMatrix)
		}
	}

	return
}

func getSliceDims(values interface{}) (ret []int) {
	switch values := values.(type) {
	case []int:
		return []int{len(values)}
	case [][]int:
		lenCurrDim := len(values)
		ret = []int{lenCurrDim}
		if lenCurrDim > 0 {
			ret = append(ret, getSliceDims(values[0])...)
		} else {
			ret = append(ret, newIntSlice(0, 1)...)
		}
	case [][][]int:
		lenCurrDim := len(values)
		ret = []int{lenCurrDim}
		if lenCurrDim > 0 {
			ret = append(ret, getSliceDims(values[0])...)
		} else {
			ret = append(ret, newIntSlice(0, 2)...)
		}
	case [][][][]int:
		lenCurrDim := len(values)
		ret = []int{lenCurrDim}
		if lenCurrDim > 0 {
			ret = append(ret, getSliceDims(values[0])...)
		} else {
			ret = append(ret, newIntSlice(0, 3)...)
		}
	case [][][][][]int:
		lenCurrDim := len(values)
		ret = []int{lenCurrDim}
		if lenCurrDim > 0 {
			ret = append(ret, getSliceDims(values[0])...)
		} else {
			ret = append(ret, newIntSlice(0, 4)...)
		}
	case [][][][][][]int:
		lenCurrDim := len(values)
		ret = []int{lenCurrDim}
		if lenCurrDim > 0 {
			ret = append(ret, getSliceDims(values[0])...)
		} else {
			ret = append(ret, newIntSlice(0, 5)...)
		}
	}
	return
}

func parseParamValues(scanner *bufio.Scanner, dims ...int) interface{} {
	switch len(dims) {
	case 0:
		panic("invalid number of dims passed to parseParamValues")
	case 1:
		return parseUntilEnoughValuesIntoSlice(scanner, dims[0])
	case 2:
		return parseParamValuesInto2DimSlice(scanner, dims[0], dims[1])
	case 3:
		return parseParamValuesInto3DimSlice(scanner, dims[0], dims[1], dims[2])
	case 4:
		return parseParamValuesInto4DimSlice(scanner, dims[0], dims[1], dims[2], dims[3])
	case 5:
		return parseParamValuesInto5DimSlice(scanner, dims[0], dims[1], dims[2], dims[3], dims[4])
	case 6:
		return parseParamValuesInto6DimSlice(scanner, dims[0], dims[1], dims[2], dims[3], dims[4], dims[5])
	}
	return nil
}

func parseFloatParamLineIntoSlice(scanner *bufio.Scanner) (ret []float64) {
	line := readUncommentedLine(scanner)
	values := strings.Fields(line)

	for _, value := range values {
		if value == "INF" {
			ret = append(ret, float64(inf))
		} else {
			ret = append(ret, parseFloat64(value))
		}
	}
	return
}

func convertToInt(values []string) (ret []int) {
	for _, value := range values {
		ret = append(ret, parseInt(value))
	}
	return
}

func convertToFloat64(values []string) (ret []float64) {
	for _, value := range values {
		ret = append(ret, parseFloat64(value))
	}
	return
}

func parseLineIntoSlice(scanner *bufio.Scanner) (values []string) {
	line := readUncommentedLine(scanner)
	values = strings.Fields(line)
	return
}

func parseUntilEnoughValuesIntoSlice(scanner *bufio.Scanner, numValuesToParse int) (ret []int) {
	totalValuesParsed := 0
	ret = make([]int, 0, numValuesToParse)
	for totalValuesParsed < numValuesToParse {
		// parsedParamsSlice := parseLineIntoSlice(scanner)
		parsedValues := parseLineIntoSlice(scanner)
		parsedIntValues := convertToInt(parsedValues)
		ret = append(ret, parsedIntValues...)
		totalValuesParsed += len(parsedIntValues)
	}
	if totalValuesParsed > numValuesToParse {
		panic("parsed too many values")
	}
	return
}

func parseInt(token string) int {
	if token == "INF" {
		return inf
	}

	valueInt64, err := strconv.ParseInt(token, 10, 0)
	if err != nil {
		panic(err)
	}
	return int(valueInt64)
}

func parseFloat64(token string) float64 {
	if token == "INF" {
		return float64(inf)
	}

	ret, err := strconv.ParseFloat(token, 64)
	if err != nil {
		panic(err)
	}
	return ret
}

func parseParamValuesInto2DimSlice(scanner *bufio.Scanner, lenDim1, lenDim2 int) (ret [][]int) {
	ret = make([][]int, lenDim1)
	for i := 0; i < lenDim1; i++ {
		ret[i] = parseUntilEnoughValuesIntoSlice(scanner, lenDim2)
	}

	return ret
}

func parseParamValuesInto3DimSlice(scanner *bufio.Scanner, lenDim1, lenDim2, lenDim3 int) (ret [][][]int) {
	ret = make([][][]int, lenDim1)
	for i := 0; i < lenDim1; i++ {
		ret[i] = parseParamValuesInto2DimSlice(scanner, lenDim2, lenDim3)
	}
	return
}

func parseParamValuesInto4DimSlice(scanner *bufio.Scanner, lenDim1, lenDim2, lenDim3, lenDim4 int) (ret [][][][]int) {
	ret = make([][][][]int, lenDim1)
	for i := 0; i < lenDim1; i++ {
		ret[i] = parseParamValuesInto3DimSlice(scanner, lenDim2, lenDim3, lenDim4)
	}
	return
}

func parseParamValuesInto5DimSlice(scanner *bufio.Scanner, lenDim1, lenDim2, lenDim3, lenDim4, lenDim5 int) (ret [][][][][]int) {
	ret = make([][][][][]int, lenDim1)
	for i := 0; i < lenDim1; i++ {
		ret[i] = parseParamValuesInto4DimSlice(scanner, lenDim2, lenDim3, lenDim4, lenDim5)
	}
	return
}

func parseParamValuesInto6DimSlice(scanner *bufio.Scanner, lenDim1, lenDim2, lenDim3, lenDim4, lenDim5, lenDim6 int) (ret [][][][][][]int) {
	ret = make([][][][][][]int, lenDim1)
	for i := 0; i < lenDim1; i++ {
		ret[i] = parseParamValuesInto5DimSlice(scanner, lenDim2, lenDim3, lenDim4, lenDim5, lenDim6)
	}
	return
}

func readUncommentedLine(scanner *bufio.Scanner) (ret string) {
	line, lineAvailable := readLine(scanner)
	for lineAvailable {
		line = removeComments(line)
		if len(strings.TrimSpace(line)) == 0 {
			// line only contain comments so continue reading
			line, lineAvailable = readLine(scanner)
		} else {
			ret = line
			return
		}
	}
	panic("no more lines to read")
}

func removeComments(source string) string {
	commentStartIdx := strings.Index(source, "/*")
	for commentStartIdx != -1 {
		commentEndIdx := strings.Index(source, "*/")
		if commentEndIdx == -1 {
			panic("unclosed comment in parameter file")
		}

		source = source[:commentStartIdx] + source[commentEndIdx+2:]

		commentStartIdx = strings.Index(source, "/*")
	}
	return source
}

func (rawEnergyParams *rawEnergyParams) setMultiLoopParams(parsedMultiLoopParams []int) {
	rawEnergyParams.multiLoopBase37C = parsedMultiLoopParams[0]
	rawEnergyParams.multiLoopBaseEnthalpy = parsedMultiLoopParams[1]
	rawEnergyParams.multiLoopClosing37C = parsedMultiLoopParams[2]
	rawEnergyParams.multiLoopClosingEnthalpy = parsedMultiLoopParams[3]
	rawEnergyParams.multiLoopIntern37C = parsedMultiLoopParams[4]
	rawEnergyParams.multiLoopInternEnthalpy = parsedMultiLoopParams[5]
}

func (rawEnergyParams *rawEnergyParams) setNinioParams(parsedNinioParams []int) {
	rawEnergyParams.ninio37C = parsedNinioParams[0]
	rawEnergyParams.ninioEnthalpy = parsedNinioParams[1]
	rawEnergyParams.maxNinio = parsedNinioParams[2]
}

func (rawEnergyParams *rawEnergyParams) setMiscParams(parsedMiscParams []float64) {
	rawEnergyParams.terminalAU37C = int(parsedMiscParams[2])
	rawEnergyParams.terminalAUEnthalpy = int(parsedMiscParams[3])
	// parameter files may or may not include the log extrapolation constant
	if len(parsedMiscParams) > 4 {
		rawEnergyParams.logExtrapolationConstant = parsedMiscParams[5]
	} else {
		// no log extrapolation constant so set to default
		rawEnergyParams.logExtrapolationConstant = defaultLogExtrapolationConstantAt37C
	}
}

type intFunc = func(int) int

// idInt is the identify func for int values
func idInt(x int) int {
	return x
}

// onlyLessThanOrEqualToZero returns x if x <= 0, else 0
func onlyLessThanOrEqualToZero(x int) int {
	return minInt(0, x)
}

func rescaleDgSlice(energy interface{}, enthalpy interface{}, temperatureInCelsius float64, fn intFunc) (ret interface{}) {
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

func rescaleDg1Dim(energy []int, enthalpy []int, temperatureInCelsius float64, fn intFunc) (ret []int) {
	lenEnergy := len(energy)
	ret = make([]int, lenEnergy)
	for i := 0; i < lenEnergy; i++ {
		val := rescaleDg(energy[i], enthalpy[i], temperatureInCelsius)
		ret[i] = fn(val)
	}
	return
}

func rescaleDg2Dim(energy [][]int, enthalpy [][]int, temperatureInCelsius float64, fn intFunc) (ret [][]int) {
	lenEnergy := len(energy)
	ret = make([][]int, lenEnergy)
	for i := 0; i < lenEnergy; i++ {
		ret[i] = rescaleDgSlice(energy[i], enthalpy[i], temperatureInCelsius, fn).([]int)
	}
	return
}

func rescaleDg3Dim(energy [][][]int, enthalpy [][][]int, temperatureInCelsius float64, fn intFunc) (ret [][][]int) {
	lenEnergy := len(energy)
	ret = make([][][]int, lenEnergy)
	for i := 0; i < lenEnergy; i++ {
		ret[i] = rescaleDgSlice(energy[i], enthalpy[i], temperatureInCelsius, fn).([][]int)
	}
	return
}

func rescaleDg4Dim(energy [][][][]int, enthalpy [][][][]int, temperatureInCelsius float64, fn intFunc) (ret [][][][]int) {
	lenEnergy := len(energy)
	ret = make([][][][]int, lenEnergy)
	for i := 0; i < lenEnergy; i++ {
		ret[i] = rescaleDgSlice(energy[i], enthalpy[i], temperatureInCelsius, fn).([][][]int)
	}
	return
}

func rescaleDg5Dim(energy [][][][][]int, enthalpy [][][][][]int, temperatureInCelsius float64, fn intFunc) (ret [][][][][]int) {
	lenEnergy := len(energy)
	ret = make([][][][][]int, lenEnergy)
	for i := 0; i < lenEnergy; i++ {
		ret[i] = rescaleDgSlice(energy[i], enthalpy[i], temperatureInCelsius, fn).([][][][]int)
	}
	return
}

func rescaleDg6Dim(energy [][][][][][]int, enthalpy [][][][][][]int, temperatureInCelsius float64, fn intFunc) (ret [][][][][][]int) {
	lenEnergy := len(energy)
	ret = make([][][][][][]int, lenEnergy)
	for i := 0; i < lenEnergy; i++ {
		ret[i] = rescaleDgSlice(energy[i], enthalpy[i], temperatureInCelsius, fn).([][][][][]int)
	}
	return
}

// scales energy paramaters according to the specificed temperatue
func (rawEnergyParams rawEnergyParams) scaleByTemperature(temperature float64) *EnergyParams {

	// set the non-matrix energy parameters
	var params *EnergyParams = &EnergyParams{
		LogExtrapolationConstant:         rescaleDgFloat64(rawEnergyParams.logExtrapolationConstant, 0, temperature),
		TerminalAUPenalty:                rescaleDg(rawEnergyParams.terminalAU37C, rawEnergyParams.terminalAUEnthalpy, temperature),
		MultiLoopUnpairedNucleotideBonus: rescaleDg(rawEnergyParams.multiLoopBase37C, rawEnergyParams.multiLoopBaseEnthalpy, temperature),
		MultiLoopClosingPenalty:          rescaleDg(rawEnergyParams.multiLoopClosing37C, rawEnergyParams.multiLoopClosingEnthalpy, temperature),
		Ninio:                            rescaleDg(rawEnergyParams.ninio37C, rawEnergyParams.ninioEnthalpy, temperature),
		MaxNinio:                         rawEnergyParams.maxNinio,
	}

	params.HairpinLoop = rescaleDgSlice(rawEnergyParams.hairpinLoopEnergy37C, rawEnergyParams.hairpinLoopEnthalpy, temperature, idInt).([]int)
	params.Bulge = rescaleDgSlice(rawEnergyParams.bulgeEnergy37C, rawEnergyParams.bulgeEnthalpy, temperature, idInt).([]int)
	params.InteriorLoop = rescaleDgSlice(rawEnergyParams.interiorLoopEnergy37C, rawEnergyParams.interiorLoopEnthalpy, temperature, idInt).([]int)

	params.MultiLoopIntern = make([]int, MaxLenLoop+1)
	for i := 0; i <= MaxLenLoop; i++ {
		params.MultiLoopIntern[i] = rescaleDg(rawEnergyParams.multiLoopIntern37C, rawEnergyParams.multiLoopInternEnthalpy, temperature)
	}

	params.TetraLoop = make(map[string]int)
	for loop := range rawEnergyParams.tetraLoopEnergy37C {
		params.TetraLoop[loop] = rescaleDg(rawEnergyParams.tetraLoopEnergy37C[loop], rawEnergyParams.tetraLoopEnthalpy[loop], temperature)
	}

	params.TriLoop = make(map[string]int)
	for loop := range rawEnergyParams.triLoopEnergy37C {
		params.TriLoop[loop] = rescaleDg(rawEnergyParams.triLoopEnergy37C[loop], rawEnergyParams.triLoopEnthalpy[loop], temperature)
	}

	params.HexaLoop = make(map[string]int)
	for loop := range rawEnergyParams.hexaLoopEnergy37C {
		params.HexaLoop[loop] = rescaleDg(rawEnergyParams.hexaLoopEnergy37C[loop], rawEnergyParams.hexaLoopEnthalpy[loop], temperature)
	}

	/* stacks    G(T) = H - [H - G(T0)]*T/T0 */
	params.StackingPair = rescaleDgSlice(rawEnergyParams.stackingPairEnergy37C, rawEnergyParams.stackingPairEnthalpy, temperature, idInt).([][]int)

	/* mismatches */
	params.MismatchInteriorLoop = rescaleDgSlice(rawEnergyParams.mismatchInteriorLoopEnergy37C, rawEnergyParams.mismatchInteriorLoopEnthalpy, temperature, idInt).([][][]int)
	params.MismatchHairpinLoop = rescaleDgSlice(rawEnergyParams.mismatchHairpinLoopEnergy37C, rawEnergyParams.mismatchHairpinLoopEnthalpy, temperature, idInt).([][][]int)
	params.Mismatch1xnInteriorLoop = rescaleDgSlice(rawEnergyParams.mismatch1xnInteriorLoopEnergy37C, rawEnergyParams.mismatch1xnInteriorLoopEnthalpy, temperature, idInt).([][][]int)
	params.Mismatch2x3InteriorLoop = rescaleDgSlice(rawEnergyParams.mismatch2x3InteriorLoopEnergy37C, rawEnergyParams.mismatch2x3InteriorLoopEnthalpy, temperature, idInt).([][][]int)

	params.MismatchMultiLoop = rescaleDgSlice(rawEnergyParams.mismatchMultiLoopEnergy37C, rawEnergyParams.mismatchMultiLoopEnthalpy, temperature, onlyLessThanOrEqualToZero).([][][]int)
	params.MismatchExteriorLoop = rescaleDgSlice(rawEnergyParams.mismatchExteriorLoopEnergy37C, rawEnergyParams.mismatchExteriorLoopEnthalpy, temperature, onlyLessThanOrEqualToZero).([][][]int)

	/* dangling ends energies */
	params.DanglingEndsFivePrime = rescaleDgSlice(rawEnergyParams.dangle5Energy37C, rawEnergyParams.dangle5Enthalpy, temperature, onlyLessThanOrEqualToZero).([][]int)
	params.DanglingEndsThreePrime = rescaleDgSlice(rawEnergyParams.dangle3Energy37C, rawEnergyParams.dangle3Enthalpy, temperature, onlyLessThanOrEqualToZero).([][]int)

	/* interior 1x1 loops */
	params.Interior1x1Loop = rescaleDgSlice(rawEnergyParams.interior1x1LoopEnergy37C, rawEnergyParams.interior1x1LoopEnthalpy, temperature, idInt).([][][][]int)

	/* interior 2x1 loops */
	params.Interior2x1Loop = rescaleDgSlice(rawEnergyParams.interior2x1LoopEnergy37C, rawEnergyParams.interior2x1LoopEnthalpy, temperature, idInt).([][][][][]int)

	/* interior 2x2 loops */
	params.Interior2x2Loop = rescaleDgSlice(rawEnergyParams.interior2x2LoopEnergy37C, rawEnergyParams.interior2x2LoopEnthalpy, temperature, idInt).([][][][][][]int)

	return params
}

/**
Rescale Gibbs free energy according to the equation dG = dH - T * dS
where dG is the change in Gibbs free energy
			dH is the change in enthalpy
			dS is the change in entrophy
			T is the temperature
*/
func rescaleDg(dG, dH int, temperature float64) int {
	// if temperate == measurementTemperatureInCelsius then below calculation will
	// always return dG. So we save some computation with this check.
	if temperature == measurementTemperatureInCelsius {
		return dG
	}

	measurementTemperatureInKelvin := measurementTemperatureInCelsius + ZeroCKelvin
	temperatureInKelvin := temperature + ZeroCKelvin
	var T float64 = float64(temperatureInKelvin / measurementTemperatureInKelvin)

	dGFloat64 := float64(dG)
	dHFloat64 := float64(dH)

	dSFloat64 := dHFloat64 - dGFloat64

	return int(dHFloat64 - dSFloat64*T)
}

// same as rescaleDg, but for floats
func rescaleDgFloat64(dG, dH, temperature float64) float64 {
	// if temperate == energyParamsTemperature then below calculation will always
	// return dG. So we save some computation with this check.
	if temperature == measurementTemperatureInCelsius {
		return dG
	}

	defaultEnergyParamsTemperatureKelvin := measurementTemperatureInCelsius + ZeroCKelvin
	temperatureKelvin := temperature + ZeroCKelvin
	var T float64 = temperatureKelvin / defaultEnergyParamsTemperatureKelvin

	dS := dH - dG
	return dH - dS*T
}

// Returns the minimum of two ints
func minInt(a, b int) int {
	if a < b {
		return a
	}
	return b
}

func maxInt(a, b int) int {
	if a > b {
		return a
	}
	return b
}

func maxIntSlice(s []int) (max int) {
	max = -inf
	for _, elem := range s {
		max = maxInt(max, elem)
	}
	return
}
