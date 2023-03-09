package seqfold

import (
	"math"

	"golang.org/x/exp/constraints"
)

func max[T constraints.Ordered](a, b T) T {
	if a > b {
		return a
	}
	return b
}

func min[T constraints.Ordered](a, b T) T {
	if a < b {
		return a
	}
	return b
}

func abs(x int) int {
	if x < 0 {
		return -x
	}
	return x
}

// RoundFloat runds at the given decimal place
func RoundFloat(number float64, decimalPlace int) float64 {
	// Calculate the 10 to the power of decimal place
	temp := math.Pow(10, float64(decimalPlace))
	// Multiply floating-point number with 10**decimalPlace and round it
	// Divide the rounded number with 10**decimalPlace to get decimal place rounding
	return math.Round(number*temp) / temp
}
