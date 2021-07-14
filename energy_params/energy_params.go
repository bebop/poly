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
	NBDistinguishableBasePairs int = 7
	// The number of distinguishable nucleotides: A, C, G, U
	NBDistinguishableNucleotides int = 4
	// The maximum loop length
	MaxLenLoop int = 30

	// ZeroCKelvin is 0 deg Celsius in Kelvin
	ZeroCKelvin float64 = 273.15
	// lxc37
	logExtrapolationConstantAt37C float64 = 107.856

	// The temperature at which the energy parameters have been measured at
	measurementTemperatureInCelcius float64 = 37.0

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
var embededEnergyParamsDirectory embed.FS
var energyParamsDirectory = "param_files"

var energyParamFileNames map[EnergyParamsSet]string = map[EnergyParamsSet]string{
	Langdon2018:    "rna_langdon2018.par",
	Andronescu2007: "rna_andronescu2007.par",
	Turner2004:     "rna_turner2004.par",
	Turner1999:     "rna_turner1999.par",
}

// rawEnergyParams is an un-exported intermediate struct used to store parsed
// energy param values. It is converted into usable energy params by the
// function `scaleByTemperature()`.
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

/**
Contains all the energy parameters needed for the free energy calculations.

The order of entries to access theses matrices always uses the closing pair or
pairs as the first indices followed by the unpaired bases in 5' to 3' direction.
For example, if we have a 2x2 interior loop:
```
		      5'-GAUA-3'
		      3'-CGCU-5'
```
The closing pairs for the loops are GC and UA (not AU!), and the unpaired bases
are (in 5' to 3' direction, starting at the first pair) A U C G.
Thus, the energy for this sequence is:
```
	pairs:                    GC UA A  U  C  G
						interior2x2Loop[2][6][1][4][2][3]
```
(See `basePairEncodedTypeMap()` and `encodeSequence()` for more details on how
pairs and unpaired nucleotides are encoded)
Note that this sequence is symmetric so the sequence is equivalent to:
```
					5'-UCGC-3'
					3'-AUAG-5'
```
which means the energy of the sequence is equivalent to:
```
	pairs:                    UA GC C  G  A  U
						interior2x2Loop[2][6][1][4][2][3]
```
*/
type EnergyParams struct {
	/**
	The matrix of free energies for stacked pairs, indexed by the two encoded closing
	pairs. The list should be formatted as symmetric a `7*7` matrix, conforming
	to the order explained above. As an example the stacked pair
	```
						5'-GU-3'
						3'-CA-5'
	```
	corresponds to the entry StackingPair[2][5] (GC=2, AU=5) which should be
	identical to StackingPair[5][2] (AU=5, GC=2).
	size: [NBDistinguishablePairs + 1][NBDistinguishablePairs + 1]int
	*/
	StackingPair [][]int
	/**
	Free energies of hairpin loops as a function of size. The list should
	contain 31 (MaxLenLoop + 1) entries. Since the minimum size of a hairpin loop
	is 3 and we start counting with 0, the first three values should be `inf` to
	indicate a forbidden value.
	size: [MaxLenLoop + 1]int
	*/
	HairpinLoop []int
	/**
	Free energies of Bulge loops. Should contain 31 (MaxLenLoop + 1) entries, the
	first one being `inf`.
	size: [MaxLenLoop + 1]int
	*/
	Bulge []int
	/**
	Free energies of interior loops. Should contain 31 (MaxLenLoop + 1) entries,
	the first 4 being `inf` (since smaller loops are tabulated).

	This field was previous called internal_loop, but has been renamed to
	interior loop to remain consistent with the names of other interior loops
	size: [MaxLenLoop + 1]int
	*/
	InteriorLoop []int
	/**
		Free energies for the interaction between the closing pair of an interior
	  loop and the two unpaired bases adjacent to the helix. This is a three
	  dimensional array indexed by the type of the closing pair and the two
		unpaired bases. Since we distinguish 5 bases the list contains
	  `7*5*5` entries. The order is such that for example the mismatch

	  ```
	  			       5'-CU-3'
	  			       3'-GC-5'
	  ```
	  corresponds to entry MismatchInteriorLoop[1][4][2] (CG=1, U=4, C=2).

		More information about mismatch energy: https://rna.urmc.rochester.edu/NNDB/turner04/tm.html
		size: [NBDistinguishablePairs + 1][NBDistinguishableNucleotides + 1][NBDistinguishableNucleotides + 1]int
	*/
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
	/**
	Energies for the interaction of an unpaired base on the 5' side and
	adjacent to a helix in multiloops and free ends (the equivalent of mismatch
	energies in interior and hairpin loops). The array is indexed by the type
	of pair closing the helix and the unpaired base and, therefore, forms a `7*5`
	matrix. For example the dangling base in
	```
							5'-GC-3'
							3'- G-5'
	```

	corresponds to entry Dangle5[1][3] (CG=1, G=3).

	More information about dangling ends: https://rna.urmc.rochester.edu/NNDB/turner04/de.html
	size: [NBDistinguishablePairs + 1][NBDistinguishableNucleotides + 1]int
	*/
	Dangle5 [][]int
	/**
	Same as above for bases on the 3' side of a helix.
	```
				       5'- A-3'
				       3'-AU-5'
	```
	corresponds to entry Dangle3[5][1] (AU=5, A=1).
	size: [NBDistinguishablePairs + 1][NBDistinguishableNucleotides + 1]int
	*/
	Dangle3 [][]int
	/**
	Free energies for symmetric size 2 interior loops. `7*7*5*5` entries.
	Example:
	```
							5'-CUU-3'
							3'-GCA-5'
	```
	corresponds to entry Interior1x1Loop[1][5][4][2] (CG=1, AU=5, U=4, C=2),
	which should be identical to Interior1x1Loop[5][1][2][4] (AU=5, CG=1, C=2, U=4).
	size: [NBDistinguishablePairs + 1][NBDistinguishablePairs + 1][NBDistinguishableNucleotides + 1][NBDistinguishableNucleotides + 1]int
	*/
	Interior1x1Loop [][][][]int
	/**
	Free energies for 2x1 interior loops, where 2 is the number of unpaired
	nucelobases on the larger 'side' of the interior loop and 1 is the number of
	unpaired nucelobases on the smaller 'side' of the interior loop.
	`7*7*5*5*5` entries.
	Example:
	```
							5'-CUUU-3'
							3'-GC A-5'
	```
	corresponds to entry Interior2x1Loop[1][5][4][4][2] (CG=1, AU=5, U=4, U=4, C=2).
	Note that this matrix is always accessed in the 5' to 3' direction with the
	larger number of unpaired nucleobases first.
	size: [NBDistinguishablePairs + 1][NBDistinguishablePairs + 1][NBDistinguishableNucleotides + 1][NBDistinguishableNucleotides + 1][NBDistinguishableNucleotides + 1]int
	*/
	Interior2x1Loop [][][][][]int
	/**
	Free energies for symmetric size 4 interior loops. `7*7*5*5*5*5` entries.
	Example:
	```
							5'-CUAU-3'
							3'-GCCA-5'
	```
	corresponds to entry Interior2x2Loop[1][5][4][1][2][2] (CG=1, AU=5, U=4, A=1, C=2, C=2),
	which should be identical to Interior2x2Loop[5][1][2][2][1][4] (AU=5, CG=1, C=2, C=2, A=1, U=4).
	size: [NBDistinguishablePairs + 1][NBDistinguishablePairs + 1][NBDistinguishableNucleotides + 1][NBDistinguishableNucleotides + 1][NBDistinguishableNucleotides + 1][NBDistinguishableNucleotides + 1]int
	*/
	Interior2x2Loop [][][][][][]int
	/**
	Parameter for logarithmic loop energy extrapolation. Used to scale energy
	parameters for hairpin, bulge and interior loops when the length of the loop
	is greather than `MaxLenLoop`.
	*/
	LogExtrapolationConstant                                                            float64
	MultiLoopUnpairedNucelotideBonus, MultiLoopClosingPenalty, TerminalAUPenalty, Ninio int
	// size: [NBDistinguishablePairs + 1]int
	MultiLoopIntern []int
	/**
	Some tetraloops particularly stable tetraloops are assigned an energy
	bonus. For example:
	```
		GAAA    -200
	```
	assigns a bonus energy of -2 kcal/mol to TetraLoop containing
	the sequence GAAA.
	*/
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
	file, err := embededEnergyParamsDirectory.Open(energyParamsDirectory + "/" + fileName)
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
				rawEnergyParams.stackingPairEnergy37C = parseParamValues(scanner, NBDistinguishableBasePairs, NBDistinguishableBasePairs).([][]int)
				rawEnergyParams.stackingPairEnergy37C = addOffset(rawEnergyParams.stackingPairEnergy37C, preOffset, 1, 1).([][]int)

			case "stack_enthalpies":
				rawEnergyParams.stackingPairEnthalpy = parseParamValues(scanner, NBDistinguishableBasePairs, NBDistinguishableBasePairs).([][]int)
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
				rawEnergyParams.mismatchExteriorLoopEnergy37C = parseParamValues(scanner, NBDistinguishableBasePairs, NBDistinguishableNucleotides+1, NBDistinguishableNucleotides+1).([][][]int)
				rawEnergyParams.mismatchExteriorLoopEnergy37C = addOffset(rawEnergyParams.mismatchExteriorLoopEnergy37C, preOffset, 1, 0, 0).([][][]int)

			case "mismatch_exterior_enthalpies":
				rawEnergyParams.mismatchExteriorLoopEnthalpy = parseParamValues(scanner, NBDistinguishableBasePairs, NBDistinguishableNucleotides+1, NBDistinguishableNucleotides+1).([][][]int)
				rawEnergyParams.mismatchExteriorLoopEnthalpy = addOffset(rawEnergyParams.mismatchExteriorLoopEnthalpy, preOffset, 1, 0, 0).([][][]int)

			case "mismatch_hairpin":
				rawEnergyParams.mismatchHairpinLoopEnergy37C = parseParamValues(scanner, NBDistinguishableBasePairs, NBDistinguishableNucleotides+1, NBDistinguishableNucleotides+1).([][][]int)
				rawEnergyParams.mismatchHairpinLoopEnergy37C = addOffset(rawEnergyParams.mismatchHairpinLoopEnergy37C, preOffset, 1, 0, 0).([][][]int)

			case "mismatch_hairpin_enthalpies":
				rawEnergyParams.mismatchHairpinLoopEnthalpy = parseParamValues(scanner, NBDistinguishableBasePairs, NBDistinguishableNucleotides+1, NBDistinguishableNucleotides+1).([][][]int)
				rawEnergyParams.mismatchHairpinLoopEnthalpy = addOffset(rawEnergyParams.mismatchHairpinLoopEnthalpy, preOffset, 1, 0, 0).([][][]int)

			case "mismatch_interior":
				rawEnergyParams.mismatchInteriorLoopEnergy37C = parseParamValues(scanner, NBDistinguishableBasePairs, NBDistinguishableNucleotides+1, NBDistinguishableNucleotides+1).([][][]int)
				rawEnergyParams.mismatchInteriorLoopEnergy37C = addOffset(rawEnergyParams.mismatchInteriorLoopEnergy37C, preOffset, 1, 0, 0).([][][]int)

			case "mismatch_interior_enthalpies":
				rawEnergyParams.mismatchInteriorLoopEnthalpy = parseParamValues(scanner, NBDistinguishableBasePairs, NBDistinguishableNucleotides+1, NBDistinguishableNucleotides+1).([][][]int)
				rawEnergyParams.mismatchInteriorLoopEnthalpy = addOffset(rawEnergyParams.mismatchInteriorLoopEnthalpy, preOffset, 1, 0, 0).([][][]int)

			case "mismatch_interior_1n":
				rawEnergyParams.mismatch1xnInteriorLoopEnergy37C = parseParamValues(scanner, NBDistinguishableBasePairs, NBDistinguishableNucleotides+1, NBDistinguishableNucleotides+1).([][][]int)
				rawEnergyParams.mismatch1xnInteriorLoopEnergy37C = addOffset(rawEnergyParams.mismatch1xnInteriorLoopEnergy37C, preOffset, 1, 0, 0).([][][]int)

			case "mismatch_interior_1n_enthalpies":
				rawEnergyParams.mismatch1xnInteriorLoopEnthalpy = parseParamValues(scanner, NBDistinguishableBasePairs, NBDistinguishableNucleotides+1, NBDistinguishableNucleotides+1).([][][]int)
				rawEnergyParams.mismatch1xnInteriorLoopEnthalpy = addOffset(rawEnergyParams.mismatch1xnInteriorLoopEnthalpy, preOffset, 1, 0, 0).([][][]int)

			case "mismatch_interior_23":
				rawEnergyParams.mismatch2x3InteriorLoopEnergy37C = parseParamValues(scanner, NBDistinguishableBasePairs, NBDistinguishableNucleotides+1, NBDistinguishableNucleotides+1).([][][]int)
				rawEnergyParams.mismatch2x3InteriorLoopEnergy37C = addOffset(rawEnergyParams.mismatch2x3InteriorLoopEnergy37C, preOffset, 1, 0, 0).([][][]int)

			case "mismatch_interior_23_enthalpies":
				rawEnergyParams.mismatch2x3InteriorLoopEnthalpy = parseParamValues(scanner, NBDistinguishableBasePairs, NBDistinguishableNucleotides+1, NBDistinguishableNucleotides+1).([][][]int)
				rawEnergyParams.mismatch2x3InteriorLoopEnthalpy = addOffset(rawEnergyParams.mismatch2x3InteriorLoopEnthalpy, preOffset, 1, 0, 0).([][][]int)

			case "mismatch_multi":
				rawEnergyParams.mismatchMultiLoopEnergy37C = parseParamValues(scanner, NBDistinguishableBasePairs, NBDistinguishableNucleotides+1, NBDistinguishableNucleotides+1).([][][]int)
				rawEnergyParams.mismatchMultiLoopEnergy37C = addOffset(rawEnergyParams.mismatchMultiLoopEnergy37C, preOffset, 1, 0, 0).([][][]int)

			case "mismatch_multi_enthalpies":
				rawEnergyParams.mismatchMultiLoopEnthalpy = parseParamValues(scanner, NBDistinguishableBasePairs, NBDistinguishableNucleotides+1, NBDistinguishableNucleotides+1).([][][]int)
				rawEnergyParams.mismatchMultiLoopEnthalpy = addOffset(rawEnergyParams.mismatchMultiLoopEnthalpy, preOffset, 1, 0, 0).([][][]int)

			case "int11":
				rawEnergyParams.interior1x1LoopEnergy37C = parseParamValues(scanner, NBDistinguishableBasePairs, NBDistinguishableBasePairs, NBDistinguishableNucleotides+1, NBDistinguishableNucleotides+1).([][][][]int)
				rawEnergyParams.interior1x1LoopEnergy37C = addOffset(rawEnergyParams.interior1x1LoopEnergy37C, preOffset, 1, 1, 0, 0).([][][][]int)

			case "int11_enthalpies":
				rawEnergyParams.interior1x1LoopEnthalpy = parseParamValues(scanner, NBDistinguishableBasePairs, NBDistinguishableBasePairs, NBDistinguishableNucleotides+1, NBDistinguishableNucleotides+1).([][][][]int)
				rawEnergyParams.interior1x1LoopEnthalpy = addOffset(rawEnergyParams.interior1x1LoopEnthalpy, preOffset, 1, 1, 0, 0).([][][][]int)

			case "int21":
				rawEnergyParams.interior2x1LoopEnergy37C = parseParamValues(scanner, NBDistinguishableBasePairs, NBDistinguishableBasePairs, NBDistinguishableNucleotides+1, NBDistinguishableNucleotides+1, NBDistinguishableNucleotides+1).([][][][][]int)
				rawEnergyParams.interior2x1LoopEnergy37C = addOffset(rawEnergyParams.interior2x1LoopEnergy37C, preOffset, 1, 1, 0, 0, 0).([][][][][]int)

			case "int21_enthalpies":
				rawEnergyParams.interior2x1LoopEnthalpy = parseParamValues(scanner, NBDistinguishableBasePairs, NBDistinguishableBasePairs, NBDistinguishableNucleotides+1, NBDistinguishableNucleotides+1, NBDistinguishableNucleotides+1).([][][][][]int)
				rawEnergyParams.interior2x1LoopEnthalpy = addOffset(rawEnergyParams.interior2x1LoopEnthalpy, preOffset, 1, 1, 0, 0, 0).([][][][][]int)

			case "int22":
				rawEnergyParams.interior2x2LoopEnergy37C = parseParamValues(scanner, NBDistinguishableBasePairs-1, NBDistinguishableBasePairs-1, NBDistinguishableNucleotides, NBDistinguishableNucleotides, NBDistinguishableNucleotides, NBDistinguishableNucleotides).([][][][][][]int)
				rawEnergyParams.interior2x2LoopEnergy37C = addOffset(rawEnergyParams.interior2x2LoopEnergy37C, preOffset, 1, 1, 1, 1, 1, 1).([][][][][][]int)
				rawEnergyParams.interior2x2LoopEnergy37C = addOffset(rawEnergyParams.interior2x2LoopEnergy37C, postOffset, 1, 1, 0, 0, 0, 0).([][][][][][]int)
				// rawEnergyParams.interior2x2LoopEnergy37C = updateNonStandardBasePairEnergyContributions(rawEnergyParams.interior2x2LoopEnergy37C)

			case "int22_enthalpies":
				rawEnergyParams.interior2x2LoopEnthalpy = parseParamValues(scanner, NBDistinguishableBasePairs-1, NBDistinguishableBasePairs-1, NBDistinguishableNucleotides, NBDistinguishableNucleotides, NBDistinguishableNucleotides, NBDistinguishableNucleotides).([][][][][][]int)
				rawEnergyParams.interior2x2LoopEnthalpy = addOffset(rawEnergyParams.interior2x2LoopEnthalpy, preOffset, 1, 1, 1, 1, 1, 1).([][][][][][]int)
				rawEnergyParams.interior2x2LoopEnthalpy = addOffset(rawEnergyParams.interior2x2LoopEnthalpy, postOffset, 1, 1, 0, 0, 0, 0).([][][][][][]int)
				// rawEnergyParams.interior2x2LoopEnthalpy = updateNonStandardBasePairEnergyContributions(rawEnergyParams.interior2x2LoopEnthalpy)

			case "dangle5":
				rawEnergyParams.dangle5Energy37C = parseParamValues(scanner, NBDistinguishableBasePairs, NBDistinguishableNucleotides+1).([][]int)
				rawEnergyParams.dangle5Energy37C = addOffset(rawEnergyParams.dangle5Energy37C, preOffset, 1, 0).([][]int)

			case "dangle5_enthalpies":
				rawEnergyParams.dangle5Enthalpy = parseParamValues(scanner, NBDistinguishableBasePairs, NBDistinguishableNucleotides+1).([][]int)
				rawEnergyParams.dangle5Enthalpy = addOffset(rawEnergyParams.dangle5Enthalpy, preOffset, 1, 0).([][]int)

			case "dangle3":
				rawEnergyParams.dangle3Energy37C = parseParamValues(scanner, NBDistinguishableBasePairs, NBDistinguishableNucleotides+1).([][]int)
				rawEnergyParams.dangle3Energy37C = addOffset(rawEnergyParams.dangle3Energy37C, preOffset, 1, 0).([][]int)

			case "dangle3_enthalpies":
				rawEnergyParams.dangle3Enthalpy = parseParamValues(scanner, NBDistinguishableBasePairs, NBDistinguishableNucleotides+1).([][]int)
				rawEnergyParams.dangle3Enthalpy = addOffset(rawEnergyParams.dangle3Enthalpy, preOffset, 1, 0).([][]int)

			case "ML_params":
				mlParams := parseParamValues(scanner, 6).([]int)
				rawEnergyParams.setMultiLoopParams(mlParams)

			case "NINIO":
				ninioParams := parseParamValues(scanner, 3).([]int)
				rawEnergyParams.setNinioParams(ninioParams)

			case "Triloops":
				rawEnergyParams.triLoops,
					rawEnergyParams.triLoopEnergy37C,
					rawEnergyParams.triLoopEnthalpy = parseTriTetraHexaLoopParams(scanner)

			case "Tetraloops":
				rawEnergyParams.tetraLoops,
					rawEnergyParams.tetraLoopEnergy37C,
					rawEnergyParams.tetraLoopEnthalpy = parseTriTetraHexaLoopParams(scanner)

			case "Hexaloops":
				rawEnergyParams.hexaLoops,
					rawEnergyParams.hexaLoopEnergy37C,
					rawEnergyParams.hexaLoopEnthalpy = parseTriTetraHexaLoopParams(scanner)

			case "Misc":
				miscParams := parseFloatParamLineIntoSlice(scanner)
				rawEnergyParams.terminalAU37C = int(miscParams[2])
				rawEnergyParams.terminalAUEnthalpy = int(miscParams[3])
				if len(miscParams) > 4 {
					rawEnergyParams.logExtrapolationConstant = miscParams[5]
				} else {
					rawEnergyParams.logExtrapolationConstant = logExtrapolationConstantAt37C
				}
			}
		}

		line, lineAvailable = readLine(scanner)
	}

	return
}

func parseTriTetraHexaLoopParams(scanner *bufio.Scanner) (loops map[string]bool, energies map[string]int, enthalpies map[string]int) {
	loops, energies, enthalpies = make(map[string]bool), make(map[string]int), make(map[string]int)
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
				loops[loop] = true
				energies[loop] = energy
				enthalpies[loop] = enthalpy
			} // else line only contain comments so continue reading

			line, lineAvailable = readLine(scanner)
		}
	}
	return
}

// returns a slice of length `length` with all values set to `inf`
func infSlice(length int) (ret []int) {
	ret = make([]int, length)
	for i := 0; i < length; i++ {
		ret[i] = inf
	}
	return
}

// adds `length` `inf`s to the front of a slice
func prependInfsToSlice(slice []int, length int) []int {
	return append(infSlice(length), slice...)
}

func appendInfsToSlice(slice []int, length int) []int {
	return append(slice, infSlice(length)...)
}

func addOffset(values interface{}, offsetType offsetType, dims ...int) interface{} {
	switch len(dims) {
	case 0:
		panic("invalid number of dims passed to parseParamValues")
	case 1:
		return addOffset1Dim(values.([]int), dims[0], offsetType)
	case 2:
		return addOffset2Dim(values.([][]int), dims[0], dims[1], offsetType)
	case 3:
		return addOffset3Dim(values.([][][]int), dims[0], dims[1], dims[2], offsetType)
	case 4:
		return addOffset4Dim(values.([][][][]int), dims[0], dims[1], dims[2], dims[3], offsetType)
	case 5:
		return addOffset5Dim(values.([][][][][]int), dims[0], dims[1], dims[2], dims[3], dims[4], offsetType)
	case 6:
		return addOffset6Dim(values.([][][][][][]int), dims[0], dims[1], dims[2], dims[3], dims[4], dims[5], offsetType)
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
	currLenDim1 := len(values)
	newLenDim1 := currLenDim1 + dim1Offset
	ret = make([][]int, 0, newLenDim1)

	var newLenDim2 int = dim2Offset
	if currLenDim1 > 0 {
		// add check to ensure we don't index an empty slice
		newLenDim2 += len(values[0])
	}

	switch offsetType {
	case preOffset:
		for i := 0; i < dim1Offset; i++ {
			infSlice := infSlice(newLenDim2)
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
			infSlice := infSlice(newLenDim2)
			ret = append(ret, infSlice)
		}
	}
	return
}

func addOffset3Dim(values [][][]int, dim1Offset, dim2Offset, dim3Offset int, offsetType offsetType) (ret [][][]int) {
	currLenDim1 := len(values)
	newLenDim1 := currLenDim1 + dim1Offset
	ret = make([][][]int, 0, newLenDim1)

	var newLenDim2 int = dim2Offset
	var currLenDim2 int = 0
	if currLenDim1 > 0 {
		// add check to ensure we don't index an empty slice
		currLenDim2 = len(values[0])
		newLenDim2 += currLenDim2
	}
	var newLenDim3 int = dim3Offset
	if currLenDim2 > 0 {
		// add check to ensure we don't index an empty slice
		newLenDim3 += len(values[0][0])
	}

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

type offsetType int

const (
	preOffset offsetType = iota
	postOffset
)

func addOffset4Dim(values [][][][]int, dim1Offset, dim2Offset, dim3Offset, dim4Offset int, offsetType offsetType) (ret [][][][]int) {
	currLenDim1 := len(values)
	newLenDim1 := currLenDim1 + dim1Offset
	ret = make([][][][]int, 0, newLenDim1)

	var newLenDim2 int = dim2Offset
	var currLenDim2 int = 0
	if currLenDim1 > 0 {
		// add check to ensure we don't index an empty slice
		currLenDim2 = len(values[0])
		newLenDim2 += currLenDim2
	}
	var newLenDim3 int = dim3Offset
	var currLenDim3 int = 0
	if currLenDim2 > 0 {
		// add check to ensure we don't index an empty slice
		currLenDim3 = len(values[0][0])
		newLenDim3 += currLenDim3
	}
	var newLenDim4 int = dim4Offset
	if currLenDim3 > 0 {
		// add check to ensure we don't index an empty slice
		newLenDim4 += len(values[0][0][0])
	}

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
	currLenDim1 := len(values)
	newLenDim1 := currLenDim1 + dim1Offset
	ret = make([][][][][]int, 0, newLenDim1)

	var newLenDim2 int = dim2Offset
	var currLenDim2 int = 0
	if currLenDim1 > 0 {
		// add check to ensure we don't index an empty slice
		currLenDim2 = len(values[0])
		newLenDim2 += currLenDim2
	}
	var newLenDim3 int = dim3Offset
	var currLenDim3 int = 0
	if currLenDim2 > 0 {
		// add check to ensure we don't index an empty slice
		currLenDim3 = len(values[0][0])
		newLenDim3 += currLenDim3
	}
	var newLenDim4 int = dim4Offset
	var currLenDim4 int = 0
	if currLenDim3 > 0 {
		// add check to ensure we don't index an empty slice
		currLenDim4 = len(values[0][0][0])
		newLenDim4 += currLenDim4
	}
	var newLenDim5 int = dim5Offset
	if currLenDim4 > 0 {
		// add check to ensure we don't index an empty slice
		newLenDim5 += len(values[0][0][0][0])
	}

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
	currLenDim1 := len(values)
	newLenDim1 := currLenDim1 + dim1Offset
	ret = make([][][][][][]int, 0, newLenDim1)

	var newLenDim2 int = dim2Offset
	var currLenDim2 int = 0
	if currLenDim1 > 0 {
		// add check to ensure we don't index an empty slice
		currLenDim2 = len(values[0])
		newLenDim2 += currLenDim2
	}
	var newLenDim3 int = dim3Offset
	var currLenDim3 int = 0
	if currLenDim2 > 0 {
		// add check to ensure we don't index an empty slice
		currLenDim3 = len(values[0][0])
		newLenDim3 += currLenDim3
	}
	var newLenDim4 int = dim4Offset
	var currLenDim4 int = 0
	if currLenDim3 > 0 {
		// add check to ensure we don't index an empty slice
		currLenDim4 = len(values[0][0][0])
		newLenDim4 += currLenDim4
	}
	var newLenDim5 int = dim5Offset
	var currLenDim5 int = 0
	if currLenDim4 > 0 {
		// add check to ensure we don't index an empty slice
		currLenDim5 = len(values[0][0][0][0])
		newLenDim5 += currLenDim5
	}
	var newLenDim6 int = dim6Offset
	if currLenDim5 > 0 {
		// add check to ensure we don't index an empty slice
		newLenDim6 += len(values[0][0][0][0][0])
	}

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

func parseLineIntoSlice(scanner *bufio.Scanner) (paramArray []int) {
	line := readUncommentedLine(scanner)
	values := strings.Fields(line)

	for _, value := range values {
		if value == "INF" {
			paramArray = append(paramArray, inf)
		} else {
			paramArray = append(paramArray, parseInt(value))
		}
	}

	return
}

func parseUntilEnoughValuesIntoSlice(scanner *bufio.Scanner, numValuesToParse int) (ret []int) {
	totalValuesParsed := 0
	ret = make([]int, 0, numValuesToParse)
	for totalValuesParsed < numValuesToParse {
		parsedParamsSlice := parseLineIntoSlice(scanner)
		ret = append(ret, parsedParamsSlice...)
		totalValuesParsed += len(parsedParamsSlice)
	}
	if totalValuesParsed > numValuesToParse {
		panic("parsed too many values")
	}
	return
}

func parseInt(token string) int {
	valueInt64, err := strconv.ParseInt(token, 10, 0)
	if err != nil {
		panic(err)
	}
	return int(valueInt64)
}

func parseFloat64(token string) (ret float64) {
	ret, err := strconv.ParseFloat(token, 64)
	if err != nil {
		panic(err)
	}
	return
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

/* update nonstandard nucleotide/basepair involved contributions for int22 */
// func updateNonStandardBasePairEnergyContributions(interior2x2Loop [][][][][][]int) [][][][][][]int {
// 	var max, max2, max3, max4, max5, max6 int

// 	/* get maxima for one nonstandard nucleotide */
// 	for i := 1; i < NBDistinguishableBasePairs; i++ {
// 		for j := 1; j < NBDistinguishableBasePairs; j++ {
// 			for k := 1; k <= NBDistinguishableNucleotides; k++ {
// 				for l := 1; l <= NBDistinguishableNucleotides; l++ {
// 					for m := 1; m <= NBDistinguishableNucleotides; m++ {
// 						/* max of {CGAU} */
// 						interior2x2Loop[i][j][k][l][m][0] = maxIntSlice(interior2x2Loop[i][j][k][l][m][:][1 : NBDistinguishableNucleotides+1])
// 						interior2x2Loop[i][j][k][l][0][m] = maxIntSlice(interior2x2Loop[i][j][k][l][:][m][1 : NBDistinguishableNucleotides+1])
// 						interior2x2Loop[i][j][k][0][l][m] = maxIntSlice(interior2x2Loop[i][j][k][:][l][m][1 : NBDistinguishableNucleotides+1])
// 						interior2x2Loop[i][j][0][k][l][m] = maxIntSlice(interior2x2Loop[i][j][:][k][l][m][1 : NBDistinguishableNucleotides+1])
// 					}
// 				}
// 			}
// 		}
// 	}
// 	/* get maxima for two nonstandard nucleotides */
// 	for i := 1; i < NBDistinguishableBasePairs; i++ {
// 		for j := 1; j < NBDistinguishableBasePairs; j++ {
// 			for k := 1; k <= NBDistinguishableNucleotides; k++ {
// 				for l := 1; l <= NBDistinguishableNucleotides; l++ {
// 					max, max2, max3, max4, max5, max6 = -inf, -inf, -inf, -inf, -inf, -inf /* max of {CGAU} */
// 					for m := 1; m <= NBDistinguishableNucleotides; m++ {
// 						max = maxInt(max, interior2x2Loop[i][j][k][l][m][0])
// 						max2 = maxInt(max2, interior2x2Loop[i][j][k][m][0][l])
// 						max3 = maxInt(max3, interior2x2Loop[i][j][m][0][k][l])
// 						max4 = maxInt(max4, interior2x2Loop[i][j][0][k][l][m])
// 						max5 = maxInt(max5, interior2x2Loop[i][j][0][k][m][l])
// 						max6 = maxInt(max6, interior2x2Loop[i][j][k][0][l][m])
// 					}
// 					interior2x2Loop[i][j][k][l][0][0] = max
// 					interior2x2Loop[i][j][k][0][0][l] = max2
// 					interior2x2Loop[i][j][0][0][k][l] = max3
// 					interior2x2Loop[i][j][k][0][l][0] = max6
// 					interior2x2Loop[i][j][0][k][0][l] = max5
// 					interior2x2Loop[i][j][0][k][l][0] = max4
// 				}
// 			}
// 		}
// 	}
// 	/* get maxima for three nonstandard nucleotides */
// 	for i := 1; i < NBDistinguishableBasePairs; i++ {
// 		for j := 1; j < NBDistinguishableBasePairs; j++ {
// 			for k := 1; k <= NBDistinguishableNucleotides; k++ {
// 				max, max2, max3, max4 = -inf, -inf, -inf, -inf /* max of {CGAU} */
// 				for l := 1; l <= NBDistinguishableNucleotides; l++ {
// 					/* should be arbitrary where index l resides in last 3 possible locations */
// 					max = maxInt(max, interior2x2Loop[i][j][k][l][0][0])
// 					max2 = maxInt(max2, interior2x2Loop[i][j][0][k][l][0])
// 					max3 = maxInt(max3, interior2x2Loop[i][j][0][0][k][l])
// 					max4 = maxInt(max4, interior2x2Loop[i][j][0][0][l][k])
// 				}
// 				interior2x2Loop[i][j][k][0][0][0] = max
// 				interior2x2Loop[i][j][0][k][0][0] = max2
// 				interior2x2Loop[i][j][0][0][k][0] = max3
// 				interior2x2Loop[i][j][0][0][0][k] = max4
// 			}
// 		}
// 	}
// 	/* get maxima for 4 nonstandard nucleotides */
// 	for i := 1; i < NBDistinguishableBasePairs; i++ {
// 		for j := 1; j < NBDistinguishableBasePairs; j++ {
// 			max = -inf /* max of {CGAU} */
// 			for k := 1; k <= NBDistinguishableNucleotides; k++ {
// 				max = maxInt(max, interior2x2Loop[i][j][k][0][0][0])
// 			}
// 			interior2x2Loop[i][j][0][0][0][0] = max
// 		}
// 	}

// 	/*
// 	 * now compute contributions for nonstandard base pairs ...
// 	 * first, 1 nonstandard bp
// 	 */
// 	for i := 1; i < NBDistinguishableBasePairs; i++ {
// 		for k := 0; k <= NBDistinguishableNucleotides; k++ {
// 			for l := 0; l <= NBDistinguishableNucleotides; l++ {
// 				for m := 0; m <= NBDistinguishableNucleotides; m++ {
// 					for n := 0; n <= NBDistinguishableNucleotides; n++ {
// 						max, max2 = -inf, -inf
// 						for j := 1; j < NBDistinguishableBasePairs; j++ {
// 							max = maxInt(max, interior2x2Loop[i][j][k][l][m][n])
// 							max2 = maxInt(max2, interior2x2Loop[j][i][k][l][m][n])
// 						}
// 						interior2x2Loop[i][NBDistinguishableBasePairs][k][l][m][n] = max
// 						interior2x2Loop[NBDistinguishableBasePairs][i][k][l][m][n] = max2
// 					}
// 				}
// 			}
// 		}
// 	}

// 	/* now 2 nst base pairs */
// 	for k := 0; k <= NBDistinguishableNucleotides; k++ {
// 		for l := 0; l <= NBDistinguishableNucleotides; l++ {
// 			for m := 0; m <= NBDistinguishableNucleotides; m++ {
// 				for n := 0; n <= NBDistinguishableNucleotides; n++ {
// 					max = -inf
// 					for j := 1; j < NBDistinguishableBasePairs; j++ {
// 						max = maxInt(max, interior2x2Loop[NBDistinguishableBasePairs][j][k][l][m][n])
// 					}
// 					interior2x2Loop[NBDistinguishableBasePairs][NBDistinguishableBasePairs][k][l][m][n] = max
// 				}
// 			}
// 		}
// 	}
// 	return interior2x2Loop
// }

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

// type alias
type intFunc = func(int) int

// idInt is the identify func for int values
func idInt(x int) int {
	return x
}

func onlyNegative(x int) int {
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

	// set the non-matix energy parameters
	var params *EnergyParams = &EnergyParams{
		LogExtrapolationConstant:         rescaleDgFloat64(rawEnergyParams.logExtrapolationConstant, 0, temperature),
		TerminalAUPenalty:                rescaleDg(rawEnergyParams.terminalAU37C, rawEnergyParams.terminalAUEnthalpy, temperature),
		MultiLoopUnpairedNucelotideBonus: rescaleDg(rawEnergyParams.multiLoopBase37C, rawEnergyParams.multiLoopBaseEnthalpy, temperature),
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
	for k := range rawEnergyParams.tetraLoops {
		params.TetraLoop[k] = rescaleDg(rawEnergyParams.tetraLoopEnergy37C[k], rawEnergyParams.tetraLoopEnthalpy[k], temperature)
	}

	params.TriLoop = make(map[string]int)
	for k := range rawEnergyParams.tetraLoops {
		params.TriLoop[k] = rescaleDg(rawEnergyParams.triLoopEnergy37C[k], rawEnergyParams.triLoopEnthalpy[k], temperature)
	}

	params.HexaLoop = make(map[string]int)
	for k := range rawEnergyParams.hexaLoops {
		params.HexaLoop[k] = rescaleDg(rawEnergyParams.hexaLoopEnergy37C[k], rawEnergyParams.hexaLoopEnthalpy[k], temperature)
	}

	/* stacks    G(T) = H - [H - G(T0)]*T/T0 */
	params.StackingPair = rescaleDgSlice(rawEnergyParams.stackingPairEnergy37C, rawEnergyParams.stackingPairEnthalpy, temperature, idInt).([][]int)

	/* mismatches */
	params.MismatchInteriorLoop = rescaleDgSlice(rawEnergyParams.mismatchInteriorLoopEnergy37C, rawEnergyParams.mismatchInteriorLoopEnthalpy, temperature, idInt).([][][]int)
	params.MismatchHairpinLoop = rescaleDgSlice(rawEnergyParams.mismatchHairpinLoopEnergy37C, rawEnergyParams.mismatchHairpinLoopEnthalpy, temperature, idInt).([][][]int)
	params.Mismatch1xnInteriorLoop = rescaleDgSlice(rawEnergyParams.mismatch1xnInteriorLoopEnergy37C, rawEnergyParams.mismatch1xnInteriorLoopEnthalpy, temperature, idInt).([][][]int)
	params.Mismatch2x3InteriorLoop = rescaleDgSlice(rawEnergyParams.mismatch2x3InteriorLoopEnergy37C, rawEnergyParams.mismatch2x3InteriorLoopEnthalpy, temperature, idInt).([][][]int)

	params.MismatchMultiLoop = rescaleDgSlice(rawEnergyParams.mismatchMultiLoopEnergy37C, rawEnergyParams.mismatchMultiLoopEnthalpy, temperature, onlyNegative).([][][]int)
	params.MismatchExteriorLoop = rescaleDgSlice(rawEnergyParams.mismatchExteriorLoopEnergy37C, rawEnergyParams.mismatchExteriorLoopEnthalpy, temperature, onlyNegative).([][][]int)

	/* dangles */
	params.Dangle5 = rescaleDgSlice(rawEnergyParams.dangle5Energy37C, rawEnergyParams.dangle5Enthalpy, temperature, onlyNegative).([][]int)
	params.Dangle3 = rescaleDgSlice(rawEnergyParams.dangle3Energy37C, rawEnergyParams.dangle3Enthalpy, temperature, onlyNegative).([][]int)

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
	// if temperate == energyParamsTemperature then below calculation will always
	// return dG. So we save some computation with this check.
	if temperature == measurementTemperatureInCelcius {
		return dG
	}

	defaultEnergyParamsTemperatureKelvin := measurementTemperatureInCelcius + ZeroCKelvin
	temperatureKelvin := temperature + ZeroCKelvin
	var T float64 = float64(temperatureKelvin / defaultEnergyParamsTemperatureKelvin)

	dGFloat64 := float64(dG)
	dHFloat64 := float64(dH)

	dSFloat64 := dHFloat64 - dGFloat64

	return int(dHFloat64 - dSFloat64*T)
}

// same as rescaleDg, but for floats
func rescaleDgFloat64(dG, dH, temperature float64) float64 {
	// if temperate == energyParamsTemperature then below calculation will always
	// return dG. So we save some computation with this check.
	if temperature == measurementTemperatureInCelcius {
		return dG
	}

	defaultEnergyParamsTemperatureKelvin := measurementTemperatureInCelcius + ZeroCKelvin
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
