package energy_params

import (
	"bufio"
	"fmt"
	"log"
	"strconv"
	"strings"
)

/******************************************************************************

This file defines the structs and functions needed to parse RNA energy
parameters as specified by the `RNAfold parameter file v2.0` file format from
[ViennaRNA](https://github.com/ViennaRNA/ViennaRNA#energy-parameters).

******************************************************************************/

const (
	// The constant used to extrapolate energy values when length of loop > `MaxLenLoop`
	defaultLogExtrapolationConstantAt37C float64 = 107.856

	// inf is infinity as used in minimization routines (INT_MAX/10)
	inf int = 10000000
)

// rawEnergyParams is an un-exported intermediate struct used to store parsed
// energy param values. It is converted into usable energy params by the
// function `scaleByTemperature()`. For more information on the fields of
// this struct, read the documentation of the fields of the `EnergyParams`
// struct.
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

	danglingEndsFivePrimeEnergy37C  [][]int
	danglingEndsFivePrimeEnthalpy   [][]int
	danglingEndsThreePrimeEnergy37C [][]int
	danglingEndsThreePrimeEnthalpy  [][]int
	stackingPairEnergy37C           [][]int
	stackingPairEnthalpy            [][]int

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

	tetraLoopEnergy37C map[string]int
	tetraLoopEnthalpy  map[string]int

	triLoopEnergy37C map[string]int
	triLoopEnthalpy  map[string]int

	hexaLoopEnergy37C map[string]int
	hexaLoopEnthalpy  map[string]int
}

// The main function where the parsing of energy param files starts from.
//
// Note that the `RNAfold parameter file v2.0` file format has a flaw.
// Since the nucleotide encoded int map (`NucleotideEncodedIntMap`) starts at
// 1 (instead of 0), all parameters specified for nucleotides in the energy
// params file have an extra first dimension even though the values can never
// be accessed.
// For example, for `mismatch_exterior`, we parse parameters of the dimesions
// `NbDistinguishableBasePairs` x `NbDistinguishableNucleotides+1` x
// `NbDistinguishableNucleotides+1`.
// The extra `+ 1` comes from the fact that `NucleotideEncodedIntMap` starts
// at 1 instead of 0. The first dimension that we parse doesn't have any
// useful information since it can't be accessed, but we still have to parse
// it since the parameter file includes the dimension. Thus, in the future,
// if we'd like to make `NucleotideEncodedIntMap` start at 0, we have to
// include a function that removes the first dimension whenever a
// `NbDistinguishableNucleotides+1` dimension is parsed.
//
// Note that for some reason int22 and int22_enthalpies parameters have the
// right dimension sizes for the nucleotide dimensions, so we have to add
// a pre offset to those nucleotide dimensions.
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
		panic("Missing header line in file.\nThis file is not of the RNAfold parameter file v2.0 format.\n")
	}

	var typeOfEnergyParamsToParse string
	for line, lineAvailable := readLine(scanner); lineAvailable; line, lineAvailable = readLine(scanner) {
		// scan the required string into energyParamsToParse
		numItemsParsed, _ := fmt.Sscanf(line, "# %255s", &typeOfEnergyParamsToParse)

		if numItemsParsed == 0 {
			// we have parsed a blank line so continue
			continue
		}

		if typeOfEnergyParamsToParse == "END" {
			// we have finished parsing the file
			break
		}

		switch typeOfEnergyParamsToParse {
		case "stack":
			rawEnergyParams.stackingPairEnergy37C = parseItemsIntoIntMatrix(scanner, NbDistinguishableBasePairs, NbDistinguishableBasePairs).([][]int)

		case "stack_enthalpies":
			rawEnergyParams.stackingPairEnthalpy = parseItemsIntoIntMatrix(scanner, NbDistinguishableBasePairs, NbDistinguishableBasePairs).([][]int)

		case "hairpin":
			rawEnergyParams.hairpinLoopEnergy37C = parseItemsIntoIntMatrix(scanner, MaxLenLoop+1).([]int)

		case "hairpin_enthalpies":
			rawEnergyParams.hairpinLoopEnthalpy = parseItemsIntoIntMatrix(scanner, MaxLenLoop+1).([]int)

		case "bulge":
			rawEnergyParams.bulgeEnergy37C = parseItemsIntoIntMatrix(scanner, MaxLenLoop+1).([]int)

		case "bulge_enthalpies":
			rawEnergyParams.bulgeEnthalpy = parseItemsIntoIntMatrix(scanner, MaxLenLoop+1).([]int)

		case "interior":
			rawEnergyParams.interiorLoopEnergy37C = parseItemsIntoIntMatrix(scanner, MaxLenLoop+1).([]int)

		case "interior_enthalpies":
			rawEnergyParams.interiorLoopEnthalpy = parseItemsIntoIntMatrix(scanner, MaxLenLoop+1).([]int)

		case "mismatch_exterior":
			rawEnergyParams.mismatchExteriorLoopEnergy37C = parseItemsIntoIntMatrix(scanner, NbDistinguishableBasePairs, NbDistinguishableNucleotides+1, NbDistinguishableNucleotides+1).([][][]int)

		case "mismatch_exterior_enthalpies":
			rawEnergyParams.mismatchExteriorLoopEnthalpy = parseItemsIntoIntMatrix(scanner, NbDistinguishableBasePairs, NbDistinguishableNucleotides+1, NbDistinguishableNucleotides+1).([][][]int)

		case "mismatch_hairpin":
			rawEnergyParams.mismatchHairpinLoopEnergy37C = parseItemsIntoIntMatrix(scanner, NbDistinguishableBasePairs, NbDistinguishableNucleotides+1, NbDistinguishableNucleotides+1).([][][]int)

		case "mismatch_hairpin_enthalpies":
			rawEnergyParams.mismatchHairpinLoopEnthalpy = parseItemsIntoIntMatrix(scanner, NbDistinguishableBasePairs, NbDistinguishableNucleotides+1, NbDistinguishableNucleotides+1).([][][]int)

		case "mismatch_interior":
			rawEnergyParams.mismatchInteriorLoopEnergy37C = parseItemsIntoIntMatrix(scanner, NbDistinguishableBasePairs, NbDistinguishableNucleotides+1, NbDistinguishableNucleotides+1).([][][]int)

		case "mismatch_interior_enthalpies":
			rawEnergyParams.mismatchInteriorLoopEnthalpy = parseItemsIntoIntMatrix(scanner, NbDistinguishableBasePairs, NbDistinguishableNucleotides+1, NbDistinguishableNucleotides+1).([][][]int)

		case "mismatch_interior_1n":
			rawEnergyParams.mismatch1xnInteriorLoopEnergy37C = parseItemsIntoIntMatrix(scanner, NbDistinguishableBasePairs, NbDistinguishableNucleotides+1, NbDistinguishableNucleotides+1).([][][]int)

		case "mismatch_interior_1n_enthalpies":
			rawEnergyParams.mismatch1xnInteriorLoopEnthalpy = parseItemsIntoIntMatrix(scanner, NbDistinguishableBasePairs, NbDistinguishableNucleotides+1, NbDistinguishableNucleotides+1).([][][]int)

		case "mismatch_interior_23":
			rawEnergyParams.mismatch2x3InteriorLoopEnergy37C = parseItemsIntoIntMatrix(scanner, NbDistinguishableBasePairs, NbDistinguishableNucleotides+1, NbDistinguishableNucleotides+1).([][][]int)

		case "mismatch_interior_23_enthalpies":
			rawEnergyParams.mismatch2x3InteriorLoopEnthalpy = parseItemsIntoIntMatrix(scanner, NbDistinguishableBasePairs, NbDistinguishableNucleotides+1, NbDistinguishableNucleotides+1).([][][]int)

		case "mismatch_multi":
			rawEnergyParams.mismatchMultiLoopEnergy37C = parseItemsIntoIntMatrix(scanner, NbDistinguishableBasePairs, NbDistinguishableNucleotides+1, NbDistinguishableNucleotides+1).([][][]int)

		case "mismatch_multi_enthalpies":
			rawEnergyParams.mismatchMultiLoopEnthalpy = parseItemsIntoIntMatrix(scanner, NbDistinguishableBasePairs, NbDistinguishableNucleotides+1, NbDistinguishableNucleotides+1).([][][]int)

		case "int11":
			rawEnergyParams.interior1x1LoopEnergy37C = parseItemsIntoIntMatrix(scanner, NbDistinguishableBasePairs, NbDistinguishableBasePairs, NbDistinguishableNucleotides+1, NbDistinguishableNucleotides+1).([][][][]int)

		case "int11_enthalpies":
			rawEnergyParams.interior1x1LoopEnthalpy = parseItemsIntoIntMatrix(scanner, NbDistinguishableBasePairs, NbDistinguishableBasePairs, NbDistinguishableNucleotides+1, NbDistinguishableNucleotides+1).([][][][]int)

		case "int21":
			rawEnergyParams.interior2x1LoopEnergy37C = parseItemsIntoIntMatrix(scanner, NbDistinguishableBasePairs, NbDistinguishableBasePairs, NbDistinguishableNucleotides+1, NbDistinguishableNucleotides+1, NbDistinguishableNucleotides+1).([][][][][]int)

		case "int21_enthalpies":
			rawEnergyParams.interior2x1LoopEnthalpy = parseItemsIntoIntMatrix(scanner, NbDistinguishableBasePairs, NbDistinguishableBasePairs, NbDistinguishableNucleotides+1, NbDistinguishableNucleotides+1, NbDistinguishableNucleotides+1).([][][][][]int)

		case "int22":
			rawEnergyParams.interior2x2LoopEnergy37C = parseItemsIntoIntMatrix(scanner, NbDistinguishableBasePairs-1, NbDistinguishableBasePairs-1, NbDistinguishableNucleotides, NbDistinguishableNucleotides, NbDistinguishableNucleotides, NbDistinguishableNucleotides).([][][][][][]int)
			rawEnergyParams.interior2x2LoopEnergy37C = addOffset(rawEnergyParams.interior2x2LoopEnergy37C, preOffset, 0, 0, 1, 1, 1, 1).([][][][][][]int)
			rawEnergyParams.interior2x2LoopEnergy37C = addOffset(rawEnergyParams.interior2x2LoopEnergy37C, postOffset, 1, 1, 0, 0, 0, 0).([][][][][][]int)

		case "int22_enthalpies":
			rawEnergyParams.interior2x2LoopEnthalpy = parseItemsIntoIntMatrix(scanner, NbDistinguishableBasePairs-1, NbDistinguishableBasePairs-1, NbDistinguishableNucleotides, NbDistinguishableNucleotides, NbDistinguishableNucleotides, NbDistinguishableNucleotides).([][][][][][]int)
			rawEnergyParams.interior2x2LoopEnthalpy = addOffset(rawEnergyParams.interior2x2LoopEnthalpy, preOffset, 0, 0, 1, 1, 1, 1).([][][][][][]int)
			rawEnergyParams.interior2x2LoopEnthalpy = addOffset(rawEnergyParams.interior2x2LoopEnthalpy, postOffset, 1, 1, 0, 0, 0, 0).([][][][][][]int)

		case "dangle5":
			rawEnergyParams.danglingEndsFivePrimeEnergy37C = parseItemsIntoIntMatrix(scanner, NbDistinguishableBasePairs, NbDistinguishableNucleotides+1).([][]int)

		case "dangle5_enthalpies":
			rawEnergyParams.danglingEndsFivePrimeEnthalpy = parseItemsIntoIntMatrix(scanner, NbDistinguishableBasePairs, NbDistinguishableNucleotides+1).([][]int)

		case "dangle3":
			rawEnergyParams.danglingEndsThreePrimeEnergy37C = parseItemsIntoIntMatrix(scanner, NbDistinguishableBasePairs, NbDistinguishableNucleotides+1).([][]int)

		case "dangle3_enthalpies":
			rawEnergyParams.danglingEndsThreePrimeEnthalpy = parseItemsIntoIntMatrix(scanner, NbDistinguishableBasePairs, NbDistinguishableNucleotides+1).([][]int)

		case "ML_params":
			mlParams := parseItemsIntoIntMatrix(scanner, 6).([]int)
			rawEnergyParams.setMultiLoopParams(mlParams)

		case "NINIO":
			ninioParams := parseItemsIntoIntMatrix(scanner, 3).([]int)
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

	return
}

/*****************************************************************************
Section: Parsing Into Matrices

The following section contains funcs needed to parse parameter values
into `int` matrices (up to 6 dimensions).
*****************************************************************************/

func parseItemsIntoIntMatrix(scanner *bufio.Scanner, dims ...int) interface{} {
	switch len(dims) {
	case 0:
		panic("invalid number of dims passed to parseParamValues")
	case 1:
		return parseUntilEnoughItemsIntoIntSlice(scanner, dims[0])
	case 2:
		return parseItemsInto2DimIntMatrix(scanner, dims[0], dims[1])
	case 3:
		return parseItemsInto3DimIntMatrix(scanner, dims[0], dims[1], dims[2])
	case 4:
		return parseItemsInto4DimIntMatrix(scanner, dims[0], dims[1], dims[2], dims[3])
	case 5:
		return parseItemsInto5DimIntMatrix(scanner, dims[0], dims[1], dims[2], dims[3], dims[4])
	case 6:
		return parseItemsInto6DimIntMatrix(scanner, dims[0], dims[1], dims[2], dims[3], dims[4], dims[5])
	}
	return nil
}

func parseUntilEnoughItemsIntoIntSlice(scanner *bufio.Scanner, numValuesToParse int) (ret []int) {
	totalValuesParsed := 0
	ret = make([]int, 0, numValuesToParse)
	for totalValuesParsed < numValuesToParse {
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

func parseLineIntoSlice(scanner *bufio.Scanner) (values []string) {
	line := readUncommentedLine(scanner)
	values = strings.Fields(line)
	return
}

// readUncommentedLine continues reading until a line without comments is read
func readUncommentedLine(scanner *bufio.Scanner) (ret string) {
	for line, lineAvailable := readLine(scanner); lineAvailable; line, lineAvailable = readLine(scanner) {
		line = removeComments(line)
		if len(strings.TrimSpace(line)) == 0 {
			// line only contain comments so continue reading
			continue
		} else {
			ret = line
			return
		}
	}
	panic("no more lines to read")
}

// readLine reads a line and returns the line and whether a line could be read
func readLine(scanner *bufio.Scanner) (line string, lineAvailable bool) {
	lineAvailable = scanner.Scan()

	if err := scanner.Err(); err != nil {
		log.Fatal(err)
	}

	if lineAvailable {
		line = scanner.Text()
	}
	return
}

// removeComments removes C-style inline comments from a string
func removeComments(source string) string {

	for commentStartIdx := strings.Index(source, "/*"); commentStartIdx != -1; commentStartIdx = strings.Index(source, "/*") {
		commentEndIdx := strings.Index(source, "*/")
		if commentEndIdx == -1 {
			panic("unclosed comment in parameter file")
		}

		// remove current comment from `source`
		source = source[:commentStartIdx] + source[commentEndIdx+2:]
	}
	return source
}

func convertToInt(values []string) (ret []int) {
	for _, value := range values {
		ret = append(ret, parseInt(value))
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

func parseItemsInto2DimIntMatrix(scanner *bufio.Scanner, lenDim1, lenDim2 int) (ret [][]int) {
	ret = make([][]int, lenDim1)
	for i := 0; i < lenDim1; i++ {
		ret[i] = parseUntilEnoughItemsIntoIntSlice(scanner, lenDim2)
	}

	return ret
}

func parseItemsInto3DimIntMatrix(scanner *bufio.Scanner, lenDim1, lenDim2, lenDim3 int) (ret [][][]int) {
	ret = make([][][]int, lenDim1)
	for i := 0; i < lenDim1; i++ {
		ret[i] = parseItemsInto2DimIntMatrix(scanner, lenDim2, lenDim3)
	}
	return
}

func parseItemsInto4DimIntMatrix(scanner *bufio.Scanner, lenDim1, lenDim2, lenDim3, lenDim4 int) (ret [][][][]int) {
	ret = make([][][][]int, lenDim1)
	for i := 0; i < lenDim1; i++ {
		ret[i] = parseItemsInto3DimIntMatrix(scanner, lenDim2, lenDim3, lenDim4)
	}
	return
}

func parseItemsInto5DimIntMatrix(scanner *bufio.Scanner, lenDim1, lenDim2, lenDim3, lenDim4, lenDim5 int) (ret [][][][][]int) {
	ret = make([][][][][]int, lenDim1)
	for i := 0; i < lenDim1; i++ {
		ret[i] = parseItemsInto4DimIntMatrix(scanner, lenDim2, lenDim3, lenDim4, lenDim5)
	}
	return
}

func parseItemsInto6DimIntMatrix(scanner *bufio.Scanner, lenDim1, lenDim2, lenDim3, lenDim4, lenDim5, lenDim6 int) (ret [][][][][][]int) {
	ret = make([][][][][][]int, lenDim1)
	for i := 0; i < lenDim1; i++ {
		ret[i] = parseItemsInto5DimIntMatrix(scanner, lenDim2, lenDim3, lenDim4, lenDim5, lenDim6)
	}
	return
}

/*****************************************************************************
End Section: Parsing Into Matrices
*****************************************************************************/

/*****************************************************************************
Section: Miscellaneous Parsing Functions

Parsing functions that don't fit in other sections.
*****************************************************************************/

func parseTriTetraHexaLoopParams(scanner *bufio.Scanner) (energies map[string]int, enthalpies map[string]int) {
	energies, enthalpies = make(map[string]int), make(map[string]int)

	for line, lineAvailable := readLine(scanner); lineAvailable; line, lineAvailable = readLine(scanner) {
		if len(strings.TrimSpace(line)) == 0 {
			// blank line encountered so we're done parsing tri-, tetra-, hexa- loop values
			return
		} else {
			line = removeComments(line)

			if len(strings.TrimSpace(line)) == 0 {
				// encountered a line with only comments so continue parsing
				continue
			}

			values := strings.Fields(line)
			if len(values) != 3 {
				panic(fmt.Sprintf("encountered incorrect number of values. expected 3, got %v", len(values)))
			}
			loop, energy, enthalpy := values[0], parseInt(values[1]), parseInt(values[2])
			energies[loop] = energy
			enthalpies[loop] = enthalpy
		}
	}
	return
}

func convertToFloat64(values []string) (ret []float64) {
	for _, value := range values {
		ret = append(ret, parseFloat64(value))
	}
	return
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

/*****************************************************************************
End Section: Miscellaneous Parsing Functions
*****************************************************************************/

/*****************************************************************************
Section: Adding Offsets to Matrices

The following section contains funcs needed to add pre and post offsets to
`int` matrices. Please read the documentation for `newRawEnergyParams` to
understand why pre and post offsets are needed.
*****************************************************************************/
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

// adds `length` `inf`s to the front of a slice
func prependInfsToSlice(slice []int, length int) []int {
	return append(newIntSlice(inf, length), slice...)
}

func appendInfsToSlice(slice []int, length int) []int {
	return append(slice, newIntSlice(inf, length)...)
}

// returns a slice of length `length` with all values set to `value`
func newIntSlice(value, length int) (ret []int) {
	ret = make([]int, length)
	for i := 0; i < length; i++ {
		ret[i] = value
	}
	return
}

func addOffset2Dim(values [][]int, dim1Offset, dim2Offset int, offsetType offsetType) (ret [][]int) {
	valuesDims := getIntMatrixDims(values)
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

// getIntMatrixDims returns the dimensions of a int matrix (up to 6 dimensions)
func getIntMatrixDims(values interface{}) (ret []int) {
	switch values := values.(type) {
	case []int:
		return []int{len(values)}

	case [][]int:
		lenCurrDim := len(values)
		ret = []int{lenCurrDim}
		if lenCurrDim > 0 {
			ret = append(ret, getIntMatrixDims(values[0])...)
		} else {
			// The rest of the dimensions are of length 0
			ret = append(ret, newIntSlice(0, 1)...)
		}

	case [][][]int:
		lenCurrDim := len(values)
		ret = []int{lenCurrDim}
		if lenCurrDim > 0 {
			ret = append(ret, getIntMatrixDims(values[0])...)
		} else {
			ret = append(ret, newIntSlice(0, 2)...)
		}

	case [][][][]int:
		lenCurrDim := len(values)
		ret = []int{lenCurrDim}
		if lenCurrDim > 0 {
			ret = append(ret, getIntMatrixDims(values[0])...)
		} else {
			ret = append(ret, newIntSlice(0, 3)...)
		}

	case [][][][][]int:
		lenCurrDim := len(values)
		ret = []int{lenCurrDim}
		if lenCurrDim > 0 {
			ret = append(ret, getIntMatrixDims(values[0])...)
		} else {
			ret = append(ret, newIntSlice(0, 4)...)
		}

	case [][][][][][]int:
		lenCurrDim := len(values)
		ret = []int{lenCurrDim}
		if lenCurrDim > 0 {
			ret = append(ret, getIntMatrixDims(values[0])...)
		} else {
			ret = append(ret, newIntSlice(0, 5)...)
		}

	}

	return
}

func addOffset3Dim(values [][][]int, dim1Offset, dim2Offset, dim3Offset int, offsetType offsetType) (ret [][][]int) {
	valuesDims := getIntMatrixDims(values)
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
	valuesDims := getIntMatrixDims(values)
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
	valuesDims := getIntMatrixDims(values)
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
	valuesDims := getIntMatrixDims(values)
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

/*****************************************************************************
End Section: Adding Offsets to Matrices
*****************************************************************************/

/*****************************************************************************
Section: Miscellaneous Funcs to Set Energy Param Values

Miscellaneous functions to set energy param values for the `rawEnergyParams`
struct
*****************************************************************************/

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

/*****************************************************************************
End Section: Miscellaneous Funcs to Set Energy Param Values
*****************************************************************************/
