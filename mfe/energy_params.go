package mfe

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
	lxc37 float64 = 107.856
	// INF is infinity as used in minimization routines (INT_MAX/10)
	INF int = 10000000
)