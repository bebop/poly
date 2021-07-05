package energy_params

import (
	"bufio"
	"embed"
	"fmt"
	"log"
	"strconv"
	"strings"
)

/**
* These energy parameters are nearest neighbor parameters for RNA folding
* compiled by the Turner group. The names of variables from the original repo
* (ViennaRNA) have been changed to be more descriptive and the original name
* of the variable has been commented out in the line above the variable
* declaration (to allow for easy comparison and updating of the variables).
* `mfe.go` includes an explaination of these variables in the declaration of the
* `energyParams` struct.
 */

const (
	// The number of distinguishable base pairs:
	// CG, GC, GU, UG, AU, UA, & non-standard (see `energy_params.md`)
	NBPairs int = 7
	// The number of distinguishable nucleotides: A, C, G, U
	NBBases int = 4
	// The maximum loop length
	MaxLenLoop int = 30

	// ZeroCKelvin is 0 deg Celsius in Kelvin
	ZeroCKelvin float64 = 273.15
	// lxc37
	defaultLXC float64 = 107.856

	// The temperature at which the energy parameters have been measured at
	measurementTemperature float64 = 37.0

	// INF is infinity as used in minimization routines (INT_MAX/10)
	INF int = 10000000
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
	lxc                      float64

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
	size: [NBPairs + 1][NBPairs + 1]int
	*/
	StackingPair [NBPairs + 1][NBPairs + 1]int
	/**
	Free energies of hairpin loops as a function of size. The list should
	contain 31 (MaxLenLoop + 1) entries. Since the minimum size of a hairpin loop
	is 3 and we start counting with 0, the first three values should be INF to
	indicate a forbidden value.
	*/
	HairpinLoop [MaxLenLoop + 1]int
	/**
	Free energies of Bulge loops. Should contain 31 (MaxLenLoop + 1) entries, the
	first one being INF.
	*/
	Bulge [MaxLenLoop + 1]int
	/**
	Free energies of interior loops. Should contain 31 (MaxLenLoop + 1) entries,
	the first 4 being INF (since smaller loops are tabulated).

	This field was previous called internal_loop, but has been renamed to
	interior loop to remain consistent with the names of other interior loops
	*/
	InteriorLoop [MaxLenLoop + 1]int
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
	*/
	MismatchInteriorLoop    [NBPairs + 1][NBBases + 1][NBBases + 1]int
	Mismatch1xnInteriorLoop [NBPairs + 1][NBBases + 1][NBBases + 1]int
	Mismatch2x3InteriorLoop [NBPairs + 1][NBBases + 1][NBBases + 1]int
	MismatchExteriorLoop    [NBPairs + 1][NBBases + 1][NBBases + 1]int
	MismatchHairpinLoop     [NBPairs + 1][NBBases + 1][NBBases + 1]int
	MismatchMultiLoop       [NBPairs + 1][NBBases + 1][NBBases + 1]int
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
	*/
	Dangle5 [NBPairs + 1][NBBases + 1]int
	/**
	Same as above for bases on the 3' side of a helix.
	```
				       5'- A-3'
				       3'-AU-5'
	```
	corresponds to entry Dangle3[5][1] (AU=5, A=1).
	*/
	Dangle3 [NBPairs + 1][NBBases + 1]int
	/**
	Free energies for symmetric size 2 interior loops. `7*7*5*5` entries.
	Example:
	```
							5'-CUU-3'
							3'-GCA-5'
	```
	corresponds to entry Interior1x1Loop[1][5][4][2] (CG=1, AU=5, U=4, C=2),
	which should be identical to Interior1x1Loop[5][1][2][4] (AU=5, CG=1, C=2, U=4).
	*/
	Interior1x1Loop [NBPairs + 1][NBPairs + 1][NBBases + 1][NBBases + 1]int
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
	*/
	Interior2x1Loop [NBPairs + 1][NBPairs + 1][NBBases + 1][NBBases + 1][NBBases + 1]int
	/**
	Free energies for symmetric size 4 interior loops. `7*7*5*5*5*5` entries.
	Example:
	```
							5'-CUAU-3'
							3'-GCCA-5'
	```
	corresponds to entry Interior2x2Loop[1][5][4][1][2][2] (CG=1, AU=5, U=4, A=1, C=2, C=2),
	which should be identical to Interior2x2Loop[5][1][2][2][1][4] (AU=5, CG=1, C=2, C=2, A=1, U=4).
	*/
	Interior2x2Loop [NBPairs + 1][NBPairs + 1][NBBases + 1][NBBases + 1][NBBases + 1][NBBases + 1]int
	/**
	Parameter for logarithmic loop energy extrapolation. Used to scale energy
	parameters for hairpin, bulge and interior loops when the length of the loop
	is greather than `MaxLenLoop`.
	*/
	LXC                                                                                 float64
	MultiLoopUnpairedNucelotideBonus, MultiLoopClosingPenalty, TerminalAUPenalty, Ninio int
	MultiLoopIntern                                                                     [NBPairs + 1]int
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

func NewEnergyParams(energyParamsSet EnergyParamsSet, temperatureInCelcius float64) *EnergyParams {
	return newRawEnergyParams(energyParamsSet).scaleByTemperature(temperatureInCelcius)
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
				rawEnergyParams.stackingPairEnergy37C = parseParamValues(scanner, NBPairs, NBPairs).([][]int)
				rawEnergyParams.stackingPairEnergy37C = addPreOffset2Dim(rawEnergyParams.stackingPairEnergy37C, 1, 1)

			case "stack_enthalpies":
				rawEnergyParams.stackingPairEnthalpy = parseParamValues(scanner, NBPairs, NBPairs).([][]int)
				rawEnergyParams.stackingPairEnthalpy = addPreOffset(rawEnergyParams.stackingPairEnthalpy, 1, 1).([][]int)

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
				rawEnergyParams.mismatchExteriorLoopEnergy37C = parseParamValues(scanner, NBPairs, NBBases+1, NBBases+1).([][][]int)
				rawEnergyParams.mismatchExteriorLoopEnergy37C = addPreOffset(rawEnergyParams.mismatchExteriorLoopEnergy37C, 1, 0, 0).([][][]int)

			case "mismatch_exterior_enthalpies":
				rawEnergyParams.mismatchExteriorLoopEnthalpy = parseParamValues(scanner, NBPairs, NBBases+1, NBBases+1).([][][]int)
				rawEnergyParams.mismatchExteriorLoopEnthalpy = addPreOffset(rawEnergyParams.mismatchExteriorLoopEnthalpy, 1, 0, 0).([][][]int)

			case "mismatch_hairpin":
				rawEnergyParams.mismatchHairpinLoopEnergy37C = parseParamValues(scanner, NBPairs, NBBases+1, NBBases+1).([][][]int)
				rawEnergyParams.mismatchHairpinLoopEnergy37C = addPreOffset(rawEnergyParams.mismatchHairpinLoopEnergy37C, 1, 0, 0).([][][]int)

			case "mismatch_hairpin_enthalpies":
				rawEnergyParams.mismatchHairpinLoopEnthalpy = parseParamValues(scanner, NBPairs, NBBases+1, NBBases+1).([][][]int)
				rawEnergyParams.mismatchHairpinLoopEnthalpy = addPreOffset(rawEnergyParams.mismatchHairpinLoopEnthalpy, 1, 0, 0).([][][]int)

			case "mismatch_interior":
				rawEnergyParams.mismatchInteriorLoopEnergy37C = parseParamValues(scanner, NBPairs, NBBases+1, NBBases+1).([][][]int)
				rawEnergyParams.mismatchInteriorLoopEnergy37C = addPreOffset(rawEnergyParams.mismatchInteriorLoopEnergy37C, 1, 0, 0).([][][]int)

			case "mismatch_interior_enthalpies":
				rawEnergyParams.mismatchInteriorLoopEnthalpy = parseParamValues(scanner, NBPairs, NBBases+1, NBBases+1).([][][]int)
				rawEnergyParams.mismatchInteriorLoopEnthalpy = addPreOffset(rawEnergyParams.mismatchInteriorLoopEnthalpy, 1, 0, 0).([][][]int)

			case "mismatch_interior_1n":
				rawEnergyParams.mismatch1xnInteriorLoopEnergy37C = parseParamValues(scanner, NBPairs, NBBases+1, NBBases+1).([][][]int)
				rawEnergyParams.mismatch1xnInteriorLoopEnergy37C = addPreOffset(rawEnergyParams.mismatch1xnInteriorLoopEnergy37C, 1, 0, 0).([][][]int)

			case "mismatch_interior_1n_enthalpies":
				rawEnergyParams.mismatch1xnInteriorLoopEnthalpy = parseParamValues(scanner, NBPairs, NBBases+1, NBBases+1).([][][]int)
				rawEnergyParams.mismatch1xnInteriorLoopEnthalpy = addPreOffset(rawEnergyParams.mismatch1xnInteriorLoopEnthalpy, 1, 0, 0).([][][]int)

			case "mismatch_interior_23":
				rawEnergyParams.mismatch2x3InteriorLoopEnergy37C = parseParamValues(scanner, NBPairs, NBBases+1, NBBases+1).([][][]int)
				rawEnergyParams.mismatch2x3InteriorLoopEnergy37C = addPreOffset(rawEnergyParams.mismatch2x3InteriorLoopEnergy37C, 1, 0, 0).([][][]int)

			case "mismatch_interior_23_enthalpies":
				rawEnergyParams.mismatch2x3InteriorLoopEnthalpy = parseParamValues(scanner, NBPairs, NBBases+1, NBBases+1).([][][]int)
				rawEnergyParams.mismatch2x3InteriorLoopEnthalpy = addPreOffset(rawEnergyParams.mismatch2x3InteriorLoopEnthalpy, 1, 0, 0).([][][]int)

			case "mismatch_multi":
				rawEnergyParams.mismatchMultiLoopEnergy37C = parseParamValues(scanner, NBPairs, NBBases+1, NBBases+1).([][][]int)
				rawEnergyParams.mismatchMultiLoopEnergy37C = addPreOffset(rawEnergyParams.mismatchMultiLoopEnergy37C, 1, 0, 0).([][][]int)

			case "mismatch_multi_enthalpies":
				rawEnergyParams.mismatchMultiLoopEnthalpy = parseParamValues(scanner, NBPairs, NBBases+1, NBBases+1).([][][]int)
				rawEnergyParams.mismatchMultiLoopEnthalpy = addPreOffset(rawEnergyParams.mismatchMultiLoopEnthalpy, 1, 0, 0).([][][]int)

			case "int11":
				rawEnergyParams.interior1x1LoopEnergy37C = parseParamValues(scanner, NBPairs, NBPairs, NBBases+1, NBBases+1).([][][][]int)
				rawEnergyParams.interior1x1LoopEnergy37C = addPreOffset(rawEnergyParams.interior1x1LoopEnergy37C, 1, 1, 0, 0).([][][][]int)

			case "int11_enthalpies":
				rawEnergyParams.interior1x1LoopEnthalpy = parseParamValues(scanner, NBPairs, NBPairs, NBBases+1, NBBases+1).([][][][]int)
				rawEnergyParams.interior1x1LoopEnthalpy = addPreOffset(rawEnergyParams.interior1x1LoopEnthalpy, 1, 1, 0, 0).([][][][]int)

			case "int21":
				rawEnergyParams.interior2x1LoopEnergy37C = parseParamValues(scanner, NBPairs, NBPairs, NBBases+1, NBBases+1, NBBases+1).([][][][][]int)
				rawEnergyParams.interior2x1LoopEnergy37C = addPreOffset(rawEnergyParams.interior2x1LoopEnergy37C, 1, 1, 0, 0, 0).([][][][][]int)

			case "int21_enthalpies":
				rawEnergyParams.interior2x1LoopEnthalpy = parseParamValues(scanner, NBPairs, NBPairs, NBBases+1, NBBases+1, NBBases+1).([][][][][]int)
				rawEnergyParams.interior2x1LoopEnthalpy = addPreOffset(rawEnergyParams.interior2x1LoopEnthalpy, 1, 1, 0, 0, 0).([][][][][]int)

			case "int22":
				rawEnergyParams.interior2x2LoopEnergy37C = parseParamValues(scanner, NBPairs-1, NBPairs-1, NBBases, NBBases, NBBases, NBBases).([][][][][][]int)
				rawEnergyParams.interior2x2LoopEnergy37C = addPreOffset(rawEnergyParams.interior2x2LoopEnergy37C, 1, 1, 1, 1, 1, 1).([][][][][][]int)
				rawEnergyParams.interior2x2LoopEnergy37C = addPostOffset(rawEnergyParams.interior2x2LoopEnergy37C, 1, 1, 0, 0, 0, 0).([][][][][][]int)
				rawEnergyParams.interior2x2LoopEnergy37C = updateNST(rawEnergyParams.interior2x2LoopEnergy37C)

			case "int22_enthalpies":
				rawEnergyParams.interior2x2LoopEnthalpy = parseParamValues(scanner, NBPairs-1, NBPairs-1, NBBases, NBBases, NBBases, NBBases).([][][][][][]int)
				rawEnergyParams.interior2x2LoopEnthalpy = addPreOffset(rawEnergyParams.interior2x2LoopEnthalpy, 1, 1, 1, 1, 1, 1).([][][][][][]int)
				rawEnergyParams.interior2x2LoopEnthalpy = addPostOffset(rawEnergyParams.interior2x2LoopEnthalpy, 1, 1, 0, 0, 0, 0).([][][][][][]int)
				rawEnergyParams.interior2x2LoopEnthalpy = updateNST(rawEnergyParams.interior2x2LoopEnthalpy)

			case "dangle5":
				rawEnergyParams.dangle5Energy37C = parseParamValues(scanner, NBPairs, NBBases+1).([][]int)
				rawEnergyParams.dangle5Energy37C = addPreOffset(rawEnergyParams.dangle5Energy37C, 1, 0).([][]int)

			case "dangle5_enthalpies":
				rawEnergyParams.dangle5Enthalpy = parseParamValues(scanner, NBPairs, NBBases+1).([][]int)
				rawEnergyParams.dangle5Enthalpy = addPreOffset(rawEnergyParams.dangle5Enthalpy, 1, 0).([][]int)

			case "dangle3":
				rawEnergyParams.dangle3Energy37C = parseParamValues(scanner, NBPairs, NBBases+1).([][]int)
				rawEnergyParams.dangle3Energy37C = addPreOffset(rawEnergyParams.dangle3Energy37C, 1, 0).([][]int)

			case "dangle3_enthalpies":
				rawEnergyParams.dangle3Enthalpy = parseParamValues(scanner, NBPairs, NBBases+1).([][]int)
				rawEnergyParams.dangle3Enthalpy = addPreOffset(rawEnergyParams.dangle3Enthalpy, 1, 0).([][]int)

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
					rawEnergyParams.lxc = miscParams[5]
				} else {
					rawEnergyParams.lxc = defaultLXC
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

// returns a slice of length numINF with all values set to INF
func infSlice(numINF int) (ret []int) {
	ret = make([]int, numINF)
	for i := 0; i < numINF; i++ {
		ret[i] = INF
	}
	return
}

// adds numINF `INF`s to the front of a slice
func prependInfsToSlice(slice []int, numINF int) []int {
	return append(infSlice(numINF), slice...)
}

func appendInfsToSlice(slice []int, numINF int) []int {
	return append(slice, infSlice(numINF)...)
}

func addPreOffset(values interface{}, dims ...int) interface{} {
	switch len(dims) {
	case 0:
		panic("invalid number of dims passed to parseParamValues")
	case 1:
		return prependInfsToSlice(values.([]int), dims[0])
	case 2:
		return addPreOffset2Dim(values.([][]int), dims[0], dims[1])
	case 3:
		return addPreOffset3Dim(values.([][][]int), dims[0], dims[1], dims[2])
	case 4:
		return addPreOffset4Dim(values.([][][][]int), dims[0], dims[1], dims[2], dims[3])
	case 5:
		return addPreOffset5Dim(values.([][][][][]int), dims[0], dims[1], dims[2], dims[3], dims[4])
	case 6:
		return addPreOffset6Dim(values.([][][][][][]int), dims[0], dims[1], dims[2], dims[3], dims[4], dims[5])
	}
	return nil
}

func addPreOffset2Dim(values [][]int, dim1Offset, dim2Offset int) (ret [][]int) {
	currLenDim1 := len(values)
	newLenDim1 := currLenDim1 + dim1Offset
	ret = make([][]int, 0, newLenDim1)

	var newLenDim2 int = dim2Offset
	if currLenDim1 > 0 {
		// add check to ensure we don't index an empty slice
		newLenDim2 += len(values[0])
	}
	for i := 0; i < dim1Offset; i++ {
		infSlice := infSlice(newLenDim2)
		ret = append(ret, infSlice)
	}

	for i := 0; i < currLenDim1; i++ {
		ret = append(ret, prependInfsToSlice(values[i], dim2Offset))
	}
	return
}

func addPreOffset3Dim(values [][][]int, dim1Offset, dim2Offset, dim3Offset int) (ret [][][]int) {
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
	for i := 0; i < dim1Offset; i++ {
		var infMatrix [][]int
		infMatrix = addPreOffset2Dim(infMatrix, newLenDim2, newLenDim3)
		ret = append(ret, infMatrix)
	}

	for i := 0; i < currLenDim1; i++ {
		ret = append(ret, addPreOffset2Dim(values[i], dim2Offset, dim3Offset))
	}
	return
}

func addPreOffset4Dim(values [][][][]int, dim1Offset, dim2Offset, dim3Offset, dim4Offset int) (ret [][][][]int) {
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
	for i := 0; i < dim1Offset; i++ {
		var infMatrix [][][]int
		infMatrix = addPreOffset3Dim(infMatrix, newLenDim2, newLenDim3, newLenDim4)
		ret = append(ret, infMatrix)
	}

	for i := 0; i < len(values); i++ {
		ret = append(ret, addPreOffset3Dim(values[i], dim2Offset, dim3Offset, dim4Offset))
	}
	return
}

func addPreOffset5Dim(values [][][][][]int, dim1Offset, dim2Offset, dim3Offset, dim4Offset, dim5Offset int) (ret [][][][][]int) {
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
	for i := 0; i < dim1Offset; i++ {
		var infMatrix [][][][]int
		infMatrix = addPreOffset4Dim(infMatrix, newLenDim2, newLenDim3, newLenDim4, newLenDim5)
		ret = append(ret, infMatrix)
	}

	for i := 0; i < len(values); i++ {
		ret = append(ret, addPreOffset4Dim(values[i], dim2Offset, dim3Offset, dim4Offset, dim5Offset))
	}
	return
}

func addPreOffset6Dim(values [][][][][][]int, dim1Offset, dim2Offset, dim3Offset, dim4Offset, dim5Offset, dim6Offset int) (ret [][][][][][]int) {
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
	for i := 0; i < dim1Offset; i++ {
		var infMatrix [][][][][]int
		infMatrix = addPreOffset5Dim(infMatrix, newLenDim2, newLenDim3, newLenDim4, newLenDim5, newLenDim6)
		ret = append(ret, infMatrix)
	}

	for i := 0; i < len(values); i++ {
		ret = append(ret, addPreOffset5Dim(values[i], dim2Offset, dim3Offset, dim4Offset, dim5Offset, dim6Offset))
	}
	return
}

// func readLineIntoIntArray(line string, arraySize int) ([]int, error) {
// 	lineSlice := readLineIntoSlice(line)
// 	var slice [4]int
// 	return slice
// }

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
			ret = append(ret, float64(INF))
		} else {
			ret = append(ret, parseFloat64(value))
		}
	}
	return
}

func parseLineIntoSlice(scanner *bufio.Scanner) (paramArray []int) {
	line := readUncommentedLine(scanner)
	values := strings.Fields(line)
	// if len(values) < lenArray {
	// 	panic("too few values to parse into array")
	// }

	for _, value := range values {
		if value == "INF" {
			paramArray = append(paramArray, INF)
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

func addPostOffset(values interface{}, dims ...int) interface{} {
	switch len(dims) {
	case 0:
		panic("invalid number of dims passed to parseParamValues")
	case 1:
		return appendInfsToSlice(values.([]int), dims[0])
	case 2:
		return addPostOffset2Dim(values.([][]int), dims[0], dims[1])
	case 3:
		return addPostOffset3Dim(values.([][][]int), dims[0], dims[1], dims[2])
	case 4:
		return addPostOffset4Dim(values.([][][][]int), dims[0], dims[1], dims[2], dims[3])
	case 5:
		return addPostOffset5Dim(values.([][][][][]int), dims[0], dims[1], dims[2], dims[3], dims[4])
	case 6:
		return addPostOffset6Dim(values.([][][][][][]int), dims[0], dims[1], dims[2], dims[3], dims[4], dims[5])
	}
	return nil
}

func addPostOffset2Dim(values [][]int, dim1Offset, dim2Offset int) (ret [][]int) {
	currLenDim1 := len(values)
	newLenDim1 := currLenDim1 + dim1Offset
	ret = make([][]int, 0, newLenDim1)

	var newLenDim2 int = dim2Offset
	if currLenDim1 > 0 {
		// add check to ensure we don't index an empty slice
		newLenDim2 += len(values[0])
	}

	for i := 0; i < currLenDim1; i++ {
		ret = append(ret, appendInfsToSlice(values[i], dim2Offset))
	}

	for i := 0; i < dim1Offset; i++ {
		infSlice := infSlice(newLenDim2)
		ret = append(ret, infSlice)
	}

	return
}

func addPostOffset3Dim(values [][][]int, dim1Offset, dim2Offset, dim3Offset int) (ret [][][]int) {
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

	for i := 0; i < currLenDim1; i++ {
		ret = append(ret, addPostOffset2Dim(values[i], dim2Offset, dim3Offset))
	}

	for i := 0; i < dim1Offset; i++ {
		var infMatrix [][]int
		infMatrix = addPostOffset2Dim(infMatrix, newLenDim2, newLenDim3)
		ret = append(ret, infMatrix)
	}
	return
}

func addPostOffset4Dim(values [][][][]int, dim1Offset, dim2Offset, dim3Offset, dim4Offset int) (ret [][][][]int) {
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

	for i := 0; i < currLenDim1; i++ {
		ret = append(ret, addPostOffset3Dim(values[i], dim2Offset, dim3Offset, dim4Offset))
	}

	for i := 0; i < dim1Offset; i++ {
		var infMatrix [][][]int
		infMatrix = addPostOffset3Dim(infMatrix, newLenDim2, newLenDim3, newLenDim4)
		ret = append(ret, infMatrix)
	}

	return
}

func addPostOffset5Dim(values [][][][][]int, dim1Offset, dim2Offset, dim3Offset, dim4Offset, dim5Offset int) (ret [][][][][]int) {
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

	for i := 0; i < currLenDim1; i++ {
		ret = append(ret, addPostOffset4Dim(values[i], dim2Offset, dim3Offset, dim4Offset, dim5Offset))
	}

	for i := 0; i < dim1Offset; i++ {
		var infMatrix [][][][]int
		infMatrix = addPostOffset4Dim(infMatrix, newLenDim2, newLenDim3, newLenDim4, newLenDim5)
		ret = append(ret, infMatrix)
	}

	return
}

func addPostOffset6Dim(values [][][][][][]int, dim1Offset, dim2Offset, dim3Offset, dim4Offset, dim5Offset, dim6Offset int) (ret [][][][][][]int) {
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

	for i := 0; i < currLenDim1; i++ {
		ret = append(ret, addPostOffset5Dim(values[i], dim2Offset, dim3Offset, dim4Offset, dim5Offset, dim6Offset))
	}

	for i := 0; i < dim1Offset; i++ {
		var infMatrix [][][][][]int
		infMatrix = addPostOffset5Dim(infMatrix, newLenDim2, newLenDim3, newLenDim4, newLenDim5, newLenDim6)
		ret = append(ret, infMatrix)
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
func updateNST(array [][][][][][]int) [][][][][][]int {
	var max, max2, max3, max4, max5, max6 int
	/* get maxima for one nonstandard nucleotide */
	for i := 1; i < NBPairs; i++ {
		for j := 1; j < NBPairs; j++ {
			for k := 1; k < 5; k++ {
				for l := 1; l < 5; l++ {
					for m := 1; m < 5; m++ {
						max, max2, max3, max4 = -INF, -INF, -INF, -INF /* max of {CGAU} */
						for n := 1; n < 5; n++ {
							max = maxInt(max, array[i][j][k][l][m][n])
							max2 = maxInt(max2, array[i][j][k][l][n][m])
							max3 = maxInt(max3, array[i][j][k][n][l][m])
							max4 = maxInt(max4, array[i][j][n][k][l][m])
						}
						array[i][j][k][l][m][0] = max
						array[i][j][k][l][0][m] = max2
						array[i][j][k][0][l][m] = max3
						array[i][j][0][k][l][m] = max4
					}
				}
			}
		}
	}
	/* get maxima for two nonstandard nucleotides */
	for i := 1; i < NBPairs; i++ {
		for j := 1; j < NBPairs; j++ {
			for k := 1; k < 5; k++ {
				for l := 1; l < 5; l++ {
					max, max2, max3, max4, max5, max6 = -INF, -INF, -INF, -INF, -INF, -INF /* max of {CGAU} */
					for m := 1; m < 5; m++ {
						max = maxInt(max, array[i][j][k][l][m][0])
						max2 = maxInt(max2, array[i][j][k][m][0][l])
						max3 = maxInt(max3, array[i][j][m][0][k][l])
						max4 = maxInt(max4, array[i][j][0][k][l][m])
						max5 = maxInt(max5, array[i][j][0][k][m][l])
						max6 = maxInt(max6, array[i][j][k][0][l][m])
					}
					array[i][j][k][l][0][0] = max
					array[i][j][k][0][0][l] = max2
					array[i][j][0][0][k][l] = max3
					array[i][j][k][0][l][0] = max6
					array[i][j][0][k][0][l] = max5
					array[i][j][0][k][l][0] = max4
				}
			}
		}
	}
	/* get maxima for three nonstandard nucleotides */
	for i := 1; i < NBPairs; i++ {
		for j := 1; j < NBPairs; j++ {
			for k := 1; k < 5; k++ {
				max, max2, max3, max4 = -INF, -INF, -INF, -INF /* max of {CGAU} */
				for l := 1; l < 5; l++ {
					/* should be arbitrary where index l resides in last 3 possible locations */
					max = maxInt(max, array[i][j][k][l][0][0])
					max2 = maxInt(max2, array[i][j][0][k][l][0])
					max3 = maxInt(max3, array[i][j][0][0][k][l])
					max4 = maxInt(max4, array[i][j][0][0][l][k])
				}
				array[i][j][k][0][0][0] = max
				array[i][j][0][k][0][0] = max2
				array[i][j][0][0][k][0] = max3
				array[i][j][0][0][0][k] = max4
			}
		}
	}
	/* get maxima for 4 nonstandard nucleotides */
	for i := 1; i < NBPairs; i++ {
		for j := 1; j < NBPairs; j++ {
			max = -INF /* max of {CGAU} */
			for k := 1; k < 5; k++ {
				max = maxInt(max, array[i][j][k][0][0][0])
			}
			array[i][j][0][0][0][0] = max
		}
	}

	/*
	 * now compute contributions for nonstandard base pairs ...
	 * first, 1 nonstandard bp
	 */
	for i := 1; i < NBPairs; i++ {
		for k := 0; k < 5; k++ {
			for l := 0; l < 5; l++ {
				for m := 0; m < 5; m++ {
					for n := 0; n < 5; n++ {
						max, max2 = -INF, -INF
						for j := 1; j < NBPairs; j++ {
							max = maxInt(max, array[i][j][k][l][m][n])
							max2 = maxInt(max2, array[j][i][k][l][m][n])
						}
						array[i][NBPairs][k][l][m][n] = max
						array[NBPairs][i][k][l][m][n] = max2
					}
				}
			}
		}
	}

	/* now 2 nst base pairs */
	for k := 0; k < 5; k++ {
		for l := 0; l < 5; l++ {
			for m := 0; m < 5; m++ {
				for n := 0; n < 5; n++ {
					max = -INF
					for j := 1; j < NBPairs; j++ {
						max = maxInt(max, array[NBPairs][j][k][l][m][n])
					}
					array[NBPairs][NBPairs][k][l][m][n] = max
				}
			}
		}
	}
	return array
}

func maxInt(a, b int) int {
	if a > b {
		return a
	}
	return b
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

// scales energy paramaters according to the specificed temperatue
func (rawEnergyParams rawEnergyParams) scaleByTemperature(temperature float64) *EnergyParams {

	// set the non-matix energy parameters
	var params *EnergyParams = &EnergyParams{
		LXC:                              rescaleDgFloat64(rawEnergyParams.lxc, 0, temperature),
		TerminalAUPenalty:                rescaleDg(rawEnergyParams.terminalAU37C, rawEnergyParams.terminalAUEnthalpy, temperature),
		MultiLoopUnpairedNucelotideBonus: rescaleDg(rawEnergyParams.multiLoopBase37C, rawEnergyParams.multiLoopBaseEnthalpy, temperature),
		MultiLoopClosingPenalty:          rescaleDg(rawEnergyParams.multiLoopClosing37C, rawEnergyParams.multiLoopClosingEnthalpy, temperature),
		Ninio:                            rescaleDg(rawEnergyParams.ninio37C, rawEnergyParams.ninioEnthalpy, temperature),
		MaxNinio:                         rawEnergyParams.maxNinio,
	}

	// scale and set hairpin, bulge and interior energy params
	for i := 0; i <= MaxLenLoop; i++ {
		params.HairpinLoop[i] = rescaleDg(rawEnergyParams.hairpinLoopEnergy37C[i], rawEnergyParams.hairpinLoopEnthalpy[i], temperature)
		params.Bulge[i] = rescaleDg(rawEnergyParams.bulgeEnergy37C[i], rawEnergyParams.bulgeEnthalpy[i], temperature)
		params.InteriorLoop[i] = rescaleDg(rawEnergyParams.interiorLoopEnergy37C[i], rawEnergyParams.interiorLoopEnthalpy[i], temperature)
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

	for i := 0; i <= NBPairs; i++ {
		params.MultiLoopIntern[i] = rescaleDg(rawEnergyParams.multiLoopIntern37C, rawEnergyParams.multiLoopInternEnthalpy, temperature)
	}

	/* stacks    G(T) = H - [H - G(T0)]*T/T0 */
	for i := 0; i <= NBPairs; i++ {
		for j := 0; j <= NBPairs; j++ {
			params.StackingPair[i][j] = rescaleDg(rawEnergyParams.stackingPairEnergy37C[i][j],
				rawEnergyParams.stackingPairEnthalpy[i][j],
				temperature)
		}
	}

	/* mismatches */
	for i := 0; i <= NBPairs; i++ {
		for j := 0; j <= NBBases; j++ {
			for k := 0; k <= NBBases; k++ {
				var mismatch int
				params.MismatchInteriorLoop[i][j][k] = rescaleDg(rawEnergyParams.mismatchInteriorLoopEnergy37C[i][j][k],
					rawEnergyParams.mismatchInteriorLoopEnthalpy[i][j][k],
					temperature)
				params.MismatchHairpinLoop[i][j][k] = rescaleDg(rawEnergyParams.mismatchHairpinLoopEnergy37C[i][j][k],
					rawEnergyParams.mismatchHairpinLoopEnthalpy[i][j][k],
					temperature)
				params.Mismatch1xnInteriorLoop[i][j][k] = rescaleDg(rawEnergyParams.mismatch1xnInteriorLoopEnergy37C[i][j][k],
					rawEnergyParams.mismatch1xnInteriorLoopEnthalpy[i][j][k],
					temperature)
				params.Mismatch2x3InteriorLoop[i][j][k] = rescaleDg(rawEnergyParams.mismatch2x3InteriorLoopEnergy37C[i][j][k],
					rawEnergyParams.mismatch2x3InteriorLoopEnthalpy[i][j][k],
					temperature)

				mismatch = rescaleDg(rawEnergyParams.mismatchMultiLoopEnergy37C[i][j][k],
					rawEnergyParams.mismatchMultiLoopEnthalpy[i][j][k],
					temperature)
				params.MismatchMultiLoop[i][j][k] = min(0, mismatch)

				mismatch = rescaleDg(rawEnergyParams.mismatchExteriorLoopEnergy37C[i][j][k],
					rawEnergyParams.mismatchExteriorLoopEnthalpy[i][j][k],
					temperature)
				params.MismatchExteriorLoop[i][j][k] = min(0, mismatch)
			}
		}
	}

	/* dangles */
	for i := 0; i <= NBPairs; i++ {
		for j := 0; j <= NBBases; j++ {
			var dd int
			dd = rescaleDg(rawEnergyParams.dangle5Energy37C[i][j],
				rawEnergyParams.dangle5Enthalpy[i][j],
				temperature)
			params.Dangle5[i][j] = min(0, dd)

			dd = rescaleDg(rawEnergyParams.dangle3Energy37C[i][j],
				rawEnergyParams.dangle3Enthalpy[i][j],
				temperature)
			params.Dangle3[i][j] = min(0, dd)
		}
	}

	/* interior 1x1 loops */
	for i := 0; i <= NBPairs; i++ {
		for j := 0; j <= NBPairs; j++ {
			for k := 0; k <= NBBases; k++ {
				for l := 0; l <= NBBases; l++ {
					params.Interior1x1Loop[i][j][k][l] = rescaleDg(rawEnergyParams.interior1x1LoopEnergy37C[i][j][k][l],
						rawEnergyParams.interior1x1LoopEnthalpy[i][j][k][l],
						temperature)
				}
			}
		}
	}

	/* interior 2x1 loops */
	for i := 0; i < NBPairs; i++ {
		for j := 0; j <= NBPairs; j++ {
			for k := 0; k <= NBBases; k++ {
				for l := 0; l <= NBBases; l++ {
					for m := 0; m <= NBBases; m++ {
						params.Interior2x1Loop[i][j][k][l][m] = rescaleDg(rawEnergyParams.interior2x1LoopEnergy37C[i][j][k][l][m],
							rawEnergyParams.interior2x1LoopEnthalpy[i][j][k][l][m],
							temperature)
					}
				}
			}
		}
	}

	/* interior 2x2 loops */
	for i := 0; i < NBPairs; i++ {
		for j := 0; j <= NBPairs; j++ {
			for k := 0; k <= NBBases; k++ {
				for l := 0; l <= NBBases; l++ {
					for m := 0; m <= NBBases; m++ {
						for n := 0; n <= NBBases; n++ {
							params.Interior2x2Loop[i][j][k][l][m][n] = rescaleDg(rawEnergyParams.interior2x2LoopEnergy37C[i][j][k][l][m][n],
								rawEnergyParams.interior2x2LoopEnthalpy[i][j][k][l][m][n],
								temperature)
						}
					}
				}
			}
		}
	}

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
	if temperature == measurementTemperature {
		return dG
	}

	defaultEnergyParamsTemperatureKelvin := measurementTemperature + ZeroCKelvin
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
	if temperature == measurementTemperature {
		return dG
	}

	defaultEnergyParamsTemperatureKelvin := measurementTemperature + ZeroCKelvin
	temperatureKelvin := temperature + ZeroCKelvin
	var T float64 = temperatureKelvin / defaultEnergyParamsTemperatureKelvin

	dS := dH - dG
	return dH - dS*T
}

// Returns the minimum of two ints
func min(a, b int) int {
	if a < b {
		return a
	}
	return b
}
