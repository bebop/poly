package mfe

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

var (
	// zeroCKelvin is 0 deg Celsius in Kelvin
	zeroCKelvin float64 = 273.15
	// lxc37
	defaultLXC float64 = 107.856
	// INF is infinity as used in minimization routines (INT_MAX/10)
	INF int = 10000000
)

const (
	Langdon2018 = iota
	Andronescu2007
	Turner2004
	Turner1999
)

//go:embed energy_params/*
var embededEnergyParamsDirectory embed.FS
var energyParamsDirectory = "energy_params"

var energyParamFileNames map[int]string = map[int]string{
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

func readLine(scanner *bufio.Scanner) (string, bool) {
	// scanner.Split(bufio.ScanLines)
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

func rawEnergyParamsFromFile(fileName string) (rawEnergyParams rawEnergyParams) {
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

	// size_t      line_no;
	// char        *line, ident[256];
	// enum parset type;
	// int         r;

	// line_no = 0

	// if ((!file_content) ||
	//     (!file_content[line_no])) {
	//   return 0
	// }

	// if (strncmp(file_content[line_no++], "## RNAfold parameter file v2.0", 30) != 0) {
	//   vrna_message_warning("Missing header line in file.\n"
	//                        "May be this file has not v2.0 format.\n"
	//                        "Use INTERRUPT-key to stop.");
	// }
	line, lineAvailable := readLine(scanner)
	for lineAvailable {
		var energyParamSet string
		n, _ := fmt.Sscanf(line, "# %255s", &energyParamSet)

		if n == 1 {
			switch energyParamSet {
			case "END":
				break
			case "stack":
				rawEnergyParams.stackingPairEnergy37C = parseParamValues(scanner, nbPairs, nbPairs).([][]int)
				rawEnergyParams.stackingPairEnergy37C = addPreOffset2Dim(rawEnergyParams.stackingPairEnergy37C, 1, 1)
				// fmt.Println(testEq(rawEnergyParams.stackingPairEnergy37C, stackingPairEnergy37C, 2))
			case "stack_enthalpies":
				rawEnergyParams.stackingPairEnthalpy = parseParamValues(scanner, nbPairs, nbPairs).([][]int)
				rawEnergyParams.stackingPairEnthalpy = addPreOffset(rawEnergyParams.stackingPairEnthalpy, 1, 1).([][]int)
				// fmt.Println(testEq(rawEnergyParams.stackingPairEnthalpy, stackingPairEnthalpy, 2))
			case "hairpin":
				rawEnergyParams.hairpinLoopEnergy37C = parseParamValues(scanner, 31).([]int)
				// fmt.Println(testEq(rawEnergyParams.hairpinLoopEnergy37C, hairpinLoopEnergy37C[:], 1))
			case "hairpin_enthalpies":
				rawEnergyParams.hairpinLoopEnthalpy = parseParamValues(scanner, 31).([]int)
				// fmt.Println(testEq(rawEnergyParams.hairpinLoopEnthalpy, hairpinLoopEnthalpy[:], 1))
			case "bulge":
				rawEnergyParams.bulgeEnergy37C = parseParamValues(scanner, 31).([]int)
				// fmt.Println(testEq(rawEnergyParams.bulgeEnergy37C, bulgeEnergy37C[:], 1))
			case "bulge_enthalpies":
				rawEnergyParams.bulgeEnthalpy = parseParamValues(scanner, 31).([]int)
				// fmt.Println(testEq(rawEnergyParams.bulgeEnthalpy, bulgeEnthalpy[:], 1))
			case "interior":
				rawEnergyParams.interiorLoopEnergy37C = parseParamValues(scanner, 31).([]int)
				// fmt.Println(testEq(rawEnergyParams.interiorLoopEnergy37C, interiorLoopEnergy37C[:], 1))
			case "interior_enthalpies":
				rawEnergyParams.interiorLoopEnthalpy = parseParamValues(scanner, 31).([]int)
				// fmt.Println(testEq(rawEnergyParams.interiorLoopEnthalpy, interiorLoopEnthalpy[:], 1))
			case "mismatch_exterior":
				rawEnergyParams.mismatchExteriorLoopEnergy37C = parseParamValues(scanner, nbPairs, nbNucleobase+1, nbNucleobase+1).([][][]int)
				rawEnergyParams.mismatchExteriorLoopEnergy37C = addPreOffset(rawEnergyParams.mismatchExteriorLoopEnergy37C, 1, 0, 0).([][][]int)
				// fmt.Println(testEq(rawEnergyParams.mismatchExteriorLoopEnergy37C, mismatchExteriorLoopEnergy37C, 3))
			case "mismatch_exterior_enthalpies":
				rawEnergyParams.mismatchExteriorLoopEnthalpy = parseParamValues(scanner, nbPairs, nbNucleobase+1, nbNucleobase+1).([][][]int)
				rawEnergyParams.mismatchExteriorLoopEnthalpy = addPreOffset(rawEnergyParams.mismatchExteriorLoopEnthalpy, 1, 0, 0).([][][]int)
				// fmt.Println(testEq(rawEnergyParams.mismatchExteriorLoopEnthalpy, mismatchExteriorLoopEnthalpy, 3))
			case "mismatch_hairpin":
				rawEnergyParams.mismatchHairpinLoopEnergy37C = parseParamValues(scanner, nbPairs, nbNucleobase+1, nbNucleobase+1).([][][]int)
				rawEnergyParams.mismatchHairpinLoopEnergy37C = addPreOffset(rawEnergyParams.mismatchHairpinLoopEnergy37C, 1, 0, 0).([][][]int)
				// fmt.Println(testEq(rawEnergyParams.mismatchHairpinLoopEnergy37C, mismatchHairpinLoopEnergy37C, 3))
			case "mismatch_hairpin_enthalpies":
				rawEnergyParams.mismatchHairpinLoopEnthalpy = parseParamValues(scanner, nbPairs, nbNucleobase+1, nbNucleobase+1).([][][]int)
				rawEnergyParams.mismatchHairpinLoopEnthalpy = addPreOffset(rawEnergyParams.mismatchHairpinLoopEnthalpy, 1, 0, 0).([][][]int)
				// fmt.Println(testEq(rawEnergyParams.mismatchHairpinLoopEnthalpy, mismatchHairpinLoopEnthalpy, 3))
			case "mismatch_interior":
				rawEnergyParams.mismatchInteriorLoopEnergy37C = parseParamValues(scanner, nbPairs, nbNucleobase+1, nbNucleobase+1).([][][]int)
				rawEnergyParams.mismatchInteriorLoopEnergy37C = addPreOffset(rawEnergyParams.mismatchInteriorLoopEnergy37C, 1, 0, 0).([][][]int)
				// fmt.Println(testEq(rawEnergyParams.mismatchInteriorLoopEnergy37C, mismatchInteriorLoopEnergy37C, 3))
			case "mismatch_interior_enthalpies":
				rawEnergyParams.mismatchInteriorLoopEnthalpy = parseParamValues(scanner, nbPairs, nbNucleobase+1, nbNucleobase+1).([][][]int)
				rawEnergyParams.mismatchInteriorLoopEnthalpy = addPreOffset(rawEnergyParams.mismatchInteriorLoopEnthalpy, 1, 0, 0).([][][]int)
				// fmt.Println(testEq(rawEnergyParams.mismatchInteriorLoopEnthalpy, mismatchInteriorLoopEnthalpy, 3))
			case "mismatch_interior_1n":
				rawEnergyParams.mismatch1xnInteriorLoopEnergy37C = parseParamValues(scanner, nbPairs, nbNucleobase+1, nbNucleobase+1).([][][]int)
				rawEnergyParams.mismatch1xnInteriorLoopEnergy37C = addPreOffset(rawEnergyParams.mismatch1xnInteriorLoopEnergy37C, 1, 0, 0).([][][]int)
				// fmt.Println(testEq(rawEnergyParams.mismatch1xnInteriorLoopEnergy37C, mismatch1xnInteriorLoopEnergy37C, 3))
			case "mismatch_interior_1n_enthalpies":
				rawEnergyParams.mismatch1xnInteriorLoopEnthalpy = parseParamValues(scanner, nbPairs, nbNucleobase+1, nbNucleobase+1).([][][]int)
				rawEnergyParams.mismatch1xnInteriorLoopEnthalpy = addPreOffset(rawEnergyParams.mismatch1xnInteriorLoopEnthalpy, 1, 0, 0).([][][]int)
				// fmt.Println(testEq(rawEnergyParams.mismatch1xnInteriorLoopEnthalpy, mismatch1xnInteriorLoopEnthalpy, 3))
			case "mismatch_interior_23":
				rawEnergyParams.mismatch2x3InteriorLoopEnergy37C = parseParamValues(scanner, nbPairs, nbNucleobase+1, nbNucleobase+1).([][][]int)
				rawEnergyParams.mismatch2x3InteriorLoopEnergy37C = addPreOffset(rawEnergyParams.mismatch2x3InteriorLoopEnergy37C, 1, 0, 0).([][][]int)
				// fmt.Println(testEq(rawEnergyParams.mismatch2x3InteriorLoopEnergy37C, mismatch2x3InteriorLoopEnergy37C, 3))
			case "mismatch_interior_23_enthalpies":
				rawEnergyParams.mismatch2x3InteriorLoopEnthalpy = parseParamValues(scanner, nbPairs, nbNucleobase+1, nbNucleobase+1).([][][]int)
				rawEnergyParams.mismatch2x3InteriorLoopEnthalpy = addPreOffset(rawEnergyParams.mismatch2x3InteriorLoopEnthalpy, 1, 0, 0).([][][]int)
				// fmt.Println(testEq(rawEnergyParams.mismatch2x3InteriorLoopEnthalpy, mismatch2x3InteriorLoopEnthalpy, 3))

			case "mismatch_multi":
				rawEnergyParams.mismatchMultiLoopEnergy37C = parseParamValues(scanner, nbPairs, nbNucleobase+1, nbNucleobase+1).([][][]int)
				rawEnergyParams.mismatchMultiLoopEnergy37C = addPreOffset(rawEnergyParams.mismatchMultiLoopEnergy37C, 1, 0, 0).([][][]int)
				// fmt.Println(testEq(rawEnergyParams.mismatchMultiLoopEnergy37C, mismatchMultiLoopEnergy37C, 3))
			case "mismatch_multi_enthalpies":
				rawEnergyParams.mismatchMultiLoopEnthalpy = parseParamValues(scanner, nbPairs, nbNucleobase+1, nbNucleobase+1).([][][]int)
				rawEnergyParams.mismatchMultiLoopEnthalpy = addPreOffset(rawEnergyParams.mismatchMultiLoopEnthalpy, 1, 0, 0).([][][]int)
				// fmt.Println(testEq(rawEnergyParams.mismatchMultiLoopEnthalpy, mismatchMultiLoopEnthalpy, 3))
			case "int11":
				rawEnergyParams.interior1x1LoopEnergy37C = parseParamValues(scanner, nbPairs, nbPairs, nbNucleobase+1, nbNucleobase+1).([][][][]int)
				rawEnergyParams.interior1x1LoopEnergy37C = addPreOffset(rawEnergyParams.interior1x1LoopEnergy37C, 1, 1, 0, 0).([][][][]int)
				// fmt.Println(testEq(rawEnergyParams.interior1x1LoopEnergy37C, interior1x1LoopEnergy37C, 4))
			case "int11_enthalpies":
				rawEnergyParams.interior1x1LoopEnthalpy = parseParamValues(scanner, nbPairs, nbPairs, nbNucleobase+1, nbNucleobase+1).([][][][]int)
				rawEnergyParams.interior1x1LoopEnthalpy = addPreOffset(rawEnergyParams.interior1x1LoopEnthalpy, 1, 1, 0, 0).([][][][]int)
				// fmt.Println(testEq(rawEnergyParams.interior1x1LoopEnthalpy, interior1x1LoopEnthalpy, 4))

			case "int21":
				rawEnergyParams.interior2x1LoopEnergy37C = parseParamValues(scanner, nbPairs, nbPairs, nbNucleobase+1, nbNucleobase+1, nbNucleobase+1).([][][][][]int)
				rawEnergyParams.interior2x1LoopEnergy37C = addPreOffset(rawEnergyParams.interior2x1LoopEnergy37C, 1, 1, 0, 0, 0).([][][][][]int)
				// fmt.Println(testEq(rawEnergyParams.interior2x1LoopEnergy37C, interior2x1LoopEnergy37, 5))
			case "int21_enthalpies":
				rawEnergyParams.interior2x1LoopEnthalpy = parseParamValues(scanner, nbPairs, nbPairs, nbNucleobase+1, nbNucleobase+1, nbNucleobase+1).([][][][][]int)
				rawEnergyParams.interior2x1LoopEnthalpy = addPreOffset(rawEnergyParams.interior2x1LoopEnthalpy, 1, 1, 0, 0, 0).([][][][][]int)
				// fmt.Println(testEq(rawEnergyParams.interior2x1LoopEnthalpy, interior2x1LoopEnthalpy, 5))
			case "int22":
				rawEnergyParams.interior2x2LoopEnergy37C = parseParamValues(scanner, nbPairs-1, nbPairs-1, nbNucleobase, nbNucleobase, nbNucleobase, nbNucleobase).([][][][][][]int)
				rawEnergyParams.interior2x2LoopEnergy37C = addPreOffset(rawEnergyParams.interior2x2LoopEnergy37C, 1, 1, 1, 1, 1, 1).([][][][][][]int)
				rawEnergyParams.interior2x2LoopEnergy37C = addPostOffset(rawEnergyParams.interior2x2LoopEnergy37C, 1, 1, 0, 0, 0, 0).([][][][][][]int)
				rawEnergyParams.interior2x2LoopEnergy37C = updateNST(rawEnergyParams.interior2x2LoopEnergy37C)
				// fmt.Println(testEq(rawEnergyParams.interior2x2LoopEnergy37C, interior2x2LoopEnergy37, 6))
			case "int22_enthalpies":
				rawEnergyParams.interior2x2LoopEnthalpy = parseParamValues(scanner, nbPairs-1, nbPairs-1, nbNucleobase, nbNucleobase, nbNucleobase, nbNucleobase).([][][][][][]int)
				rawEnergyParams.interior2x2LoopEnthalpy = addPreOffset(rawEnergyParams.interior2x2LoopEnthalpy, 1, 1, 1, 1, 1, 1).([][][][][][]int)
				rawEnergyParams.interior2x2LoopEnthalpy = addPostOffset(rawEnergyParams.interior2x2LoopEnthalpy, 1, 1, 0, 0, 0, 0).([][][][][][]int)
				rawEnergyParams.interior2x2LoopEnthalpy = updateNST(rawEnergyParams.interior2x2LoopEnthalpy)
				// fmt.Println(testEq(rawEnergyParams.interior2x2LoopEnthalpy, interior2x2LoopEnthalpy, 6))
			case "dangle5":
				rawEnergyParams.dangle5Energy37C = parseParamValues(scanner, nbPairs, nbNucleobase+1).([][]int)
				rawEnergyParams.dangle5Energy37C = addPreOffset(rawEnergyParams.dangle5Energy37C, 1, 0).([][]int)
				// fmt.Println(testEq(rawEnergyParams.dangle5Energy37C, dangle5Energy37C, 2))
			case "dangle5_enthalpies":
				rawEnergyParams.dangle5Enthalpy = parseParamValues(scanner, nbPairs, nbNucleobase+1).([][]int)
				rawEnergyParams.dangle5Enthalpy = addPreOffset(rawEnergyParams.dangle5Enthalpy, 1, 0).([][]int)
				// fmt.Println(testEq(rawEnergyParams.dangle5Enthalpy, dangle5Enthalpy, 2))
			case "dangle3":
				rawEnergyParams.dangle3Energy37C = parseParamValues(scanner, nbPairs, nbNucleobase+1).([][]int)
				rawEnergyParams.dangle3Energy37C = addPreOffset(rawEnergyParams.dangle3Energy37C, 1, 0).([][]int)
				// fmt.Println(testEq(rawEnergyParams.dangle3Energy37C, dangle3Energy37C, 2))
			case "dangle3_enthalpies":
				rawEnergyParams.dangle3Enthalpy = parseParamValues(scanner, nbPairs, nbNucleobase+1).([][]int)
				rawEnergyParams.dangle3Enthalpy = addPreOffset(rawEnergyParams.dangle3Enthalpy, 1, 0).([][]int)
				// fmt.Println(testEq(rawEnergyParams.dangle3Enthalpy, dangle3Enthalpy, 2))

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
		// else { // skip all other lines
		// 	if len(strings.TrimSpace(removeComments(line))) == 0 {
		// 		fmt.Println(line)
		// 	}
		// }
		line, lineAvailable = readLine(scanner)
	}

	// fmt.Println("and the answer is...")
	// fmt.Println(reflect.DeepEqual(turner2004Params(), rawEnergyParams))
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

func printDims(arr interface{}, dims int) {
	switch dims {
	case 0:
		fmt.Println()
	case 1:
		array := arr.([]int)
		fmt.Printf("%v ", len(array))
		printDims(array[0], dims-1)
	case 2:
		array := arr.([][]int)
		fmt.Printf("%v ", len(array))
		printDims(array[0], dims-1)
	case 3:
		array := arr.([][][]int)
		fmt.Printf("%v ", len(array))
		printDims(array[0], dims-1)
	case 4:
		array := arr.([][][][]int)
		fmt.Printf("%v ", len(array))
		printDims(array[0], dims-1)
	case 5:
		array := arr.([][][][][]int)
		fmt.Printf("%v ", len(array))
		printDims(array[0], dims-1)
	case 6:
		array := arr.([][][][][][]int)
		fmt.Printf("%v ", len(array))
		printDims(array[0], dims-1)
	}
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

// func parseUntilEnoughValuesIntoSlice(scanner *bufio.Scanner, numValuesToParse int) (ret []int) {
// 	scanner.Split(bufio.ScanWords)
// 	totalValuesParsed := 0
// 	ret = make([]int, 0, numValuesToParse)
// 	for totalValuesParsed < numValuesToParse {
// 		scanner.Scan()
// 		token := scanner.Text()
// 		if token == "/*" {
// 			readUntilCommentEnd(scanner)
// 		} else {
// 			ret = append(ret, parseToken(token))
// 			totalValuesParsed++
// 		}
// 		// parsedParamsSlice := parseLineIntoSlice(scanner)
// 		// ret = append(ret, parsedParamsSlice...)
// 		// totalValuesParsed += len(parsedParamsSlice)
// 	}
// 	if totalValuesParsed > numValuesToParse {
// 		panic("parsed too many values")
// 	}
// 	return
// }

// func readUntilCommentEnd(scanner *bufio.Scanner) {
// 	scanner.Split(bufio.ScanWords)
// 	token := ""
// 	for token != "*/" {
// 		scanner.Scan()
// 		token = scanner.Text()
// 	}
// }

// func parseToken(token string) int {
// 	if token == "INF" {
// 		return INF
// 	} else {
// 		valueInt64, err := strconv.ParseInt(token, 10, 0)
// 		if err != nil {
// 			panic(err)
// 		}
// 		return int(valueInt64)
// 	}
// }

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
// array[nbPairs + 1][nbPairs + 1][5][5][5][5]
func updateNST(arrayPointer [][][][][][]int) (array [][][][][][]int) {
	array = arrayPointer
	var max, max2, max3, max4, max5, max6 int
	/* get maxima for one nonstandard nucleotide */
	for i := 1; i < nbPairs; i++ {
		for j := 1; j < nbPairs; j++ {
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
	for i := 1; i < nbPairs; i++ {
		for j := 1; j < nbPairs; j++ {
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
	for i := 1; i < nbPairs; i++ {
		for j := 1; j < nbPairs; j++ {
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
	for i := 1; i < nbPairs; i++ {
		for j := 1; j < nbPairs; j++ {
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
	for i := 1; i < nbPairs; i++ {
		for k := 0; k < 5; k++ {
			for l := 0; l < 5; l++ {
				for m := 0; m < 5; m++ {
					for n := 0; n < 5; n++ {
						max, max2 = -INF, -INF
						for j := 1; j < nbPairs; j++ {
							max = maxInt(max, array[i][j][k][l][m][n])
							max2 = maxInt(max2, array[j][i][k][l][m][n])
						}
						array[i][nbPairs][k][l][m][n] = max
						array[nbPairs][i][k][l][m][n] = max2
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
					for j := 1; j < nbPairs; j++ {
						max = maxInt(max, array[nbPairs][j][k][l][m][n])
					}
					array[nbPairs][nbPairs][k][l][m][n] = max
				}
			}
		}
	}
	return
}

func maxInt(a, b int) int {
	if a > b {
		return a
	}
	return b
}

// func rd_2dim_slice(content []string, line_no int, array []int, dim [2]int, shift [2]int) {
// 	var post [2]int = int{0, 0}
// 	delta_pre := shift[0] + shift[1]
// 	delta_post := post[0] + post[1]

// 	if delta_pre + delta_post == 0 {
// 		rd_1dim(content, line_no, array, dim[0]*dim[1], 0)
// 		return
// 	}

// 	for i := shift[0]; i < dim[0]-post[0]; i++ {
// 		rd_1dim_slice(content, line_no, array+(i*dim[1]), dim[1], shift[1], post[1])
// 	}
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

func testEq1DimSlice(a, b []int) bool {
	if len(a) != len(b) {
		return false
	}
	for i := range a {
		if a[i] != b[i] {
			return false
		}
	}
	return true
}

func testEq2DimSlice(a, b [][]int) bool {
	if len(a) != len(b) {
		return false
	}
	for i := range a {
		if !testEq1DimSlice(a[i], b[i]) {
			return false
		}
	}
	return true
}

func testEq3DimSlice(a, b [][][]int) bool {
	if len(a) != len(b) {
		return false
	}
	for i := range a {
		if !testEq2DimSlice(a[i], b[i]) {
			return false
		}
	}
	return true
}

func testEq4DimSlice(a, b [][][][]int) bool {
	if len(a) != len(b) {
		return false
	}
	for i := range a {
		if !testEq3DimSlice(a[i], b[i]) {
			return false
		}
	}
	return true
}

func testEq5DimSlice(a, b [][][][][]int) bool {
	if len(a) != len(b) {
		return false
	}
	for i := range a {
		if !testEq4DimSlice(a[i], b[i]) {
			return false
		}
	}
	return true
}

func testEq6DimSlice(a, b [][][][][][]int) bool {
	if len(a) != len(b) {
		return false
	}
	for i := range a {
		if !testEq5DimSlice(a[i], b[i]) {
			return false
		}
	}
	return true
}

func testEq(a, b interface{}, dims int) bool {
	switch dims {
	case 1:
		return testEq1DimSlice(a.([]int), b.([]int))
	case 2:
		return testEq2DimSlice(a.([][]int), b.([][]int))
	case 3:
		return testEq3DimSlice(a.([][][]int), b.([][][]int))
	case 4:
		return testEq4DimSlice(a.([][][][]int), b.([][][][]int))
	case 5:
		return testEq5DimSlice(a.([][][][][]int), b.([][][][][]int))
	case 6:
		return testEq6DimSlice(a.([][][][][][]int), b.([][][][][][]int))
	}
	return false
}

func getUncommentedLines(scanner *bufio.Scanner, numLines int) (ret []string) {
	line, lineAvailable := readLine(scanner)
	numUncommentedLines := 0

	for lineAvailable {
		line = removeComments(line)
		if len(strings.TrimSpace(line)) != 0 {
			ret = append(ret, line)
			numUncommentedLines++
		} // skip white space lines and full line comments
	}

	return
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
