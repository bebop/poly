/*
Package variants contains a function for generating all variants of a sequence.

Sometimes sequencers will only give you an *estimate* of what the basepair at
a given position is. This package provides a function for generating all
possible deterministic variants of a sequence given a sequence
with ambiguous bases.
*/
package variants

import (
	"errors"
	"strings"
)

// AllVariantsIUPAC takes a string as input
// and returns all iupac variants as output
func AllVariantsIUPAC(seq string) ([]string, error) {
	seqVariantList := [][]rune{}
	seqVariants := []string{}

	iupac := map[rune][]rune{ // rune map of all iupac nucleotide variants
		'G': []rune{'G'},
		'A': []rune{'A'},
		'T': []rune{'T'},
		'C': []rune{'C'},
		'R': []rune{'G', 'A'},
		'Y': []rune{'T', 'C'},
		'M': []rune{'A', 'C'},
		'K': []rune{'G', 'T'},
		'S': []rune{'G', 'C'},
		'W': []rune{'A', 'T'},
		'H': []rune{'A', 'C', 'T'},
		'B': []rune{'G', 'T', 'C'},
		'V': []rune{'G', 'C', 'A'},
		'D': []rune{'G', 'A', 'T'},
		'N': []rune{'G', 'A', 'T', 'C'},
	}

	for _, s := range strings.ToUpper(seq) {
		variantsIUPAC, ok := iupac[s]
		if ok {
			seqVariantList = append(seqVariantList, variantsIUPAC)
		} else {
			return seqVariants, errors.New("Error:" + string(s) + " is not a supported IUPAC character")
		}

	}

	cartesianProducts := cartRune(seqVariantList...)
	for _, product := range cartesianProducts {
		seqVariants = append(seqVariants, string(product))
	}
	return seqVariants, nil
}
func cartRune(inList ...[]rune) [][]rune {
	// An iteratitive approach to calculate Cartesian product of two or more lists
	// Adapted from https://rosettacode.org/wiki/Cartesian_product_of_two_or_more_lists
	// supposedly "minimizes allocations and computes and fills the result sequentially"

	var possibleVariants int = 1 // a counter used to determine the possible number of variants
	for _, inList := range inList {
		possibleVariants *= len(inList)
	}
	if possibleVariants == 0 {
		return nil // in the future this could be part of an error return?
	}
	allVariants := make([][]rune, possibleVariants)              // this is the 2D slice where all variants will be stored
	variantHolders := make([]rune, possibleVariants*len(inList)) // this is an empty slice with a length totaling the size of all input characters
	variantChoices := make([]int, len(inList))                   // these will be all the possible variants
	start := 0
	for variant := range allVariants {
		end := start + len(inList) // define end point
		variantHolder := variantHolders[start:end]

		allVariants[variant] = variantHolder

		start = end // start at end point

		for variantChoicesIndex, variantChoice := range variantChoices {
			variantHolder[variantChoicesIndex] = inList[variantChoicesIndex][variantChoice]
		}
		for variantChoicesIndex := len(variantChoices) - 1; variantChoicesIndex >= 0; variantChoicesIndex-- {
			variantChoices[variantChoicesIndex]++
			if variantChoices[variantChoicesIndex] < len(inList[variantChoicesIndex]) {
				break
			}
			variantChoices[variantChoicesIndex] = 0
		}
	}
	return allVariants
}
