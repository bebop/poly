package energy_params

/******************************************************************************

This file defines the functions needed to scale free energy params by a
specified temperature.

******************************************************************************/

const (
	// The temperature at which the energy parameters have been measured at
	measurementTemperatureInCelsius float64 = 37.0
)

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
