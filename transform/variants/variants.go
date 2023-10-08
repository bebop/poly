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

var iupacToBases = map[rune][]rune{ // rune map of all iupac nucleotide variants
	'G': {'G'},
	'A': {'A'},
	'T': {'T'},
	'C': {'C'},
	'R': {'G', 'A'},
	'Y': {'T', 'C'},
	'M': {'A', 'C'},
	'K': {'G', 'T'},
	'S': {'G', 'C'},
	'W': {'A', 'T'},
	'H': {'A', 'C', 'T'},
	'B': {'G', 'T', 'C'},
	'V': {'G', 'C', 'A'},
	'D': {'G', 'A', 'T'},
	'N': {'G', 'A', 'T', 'C'},
}

var iupacToAAs = map[rune][]rune{
	'X': []rune("ACDEFGHIJLMNPQRSTVWY"), // Unknown
	'B': {'D', 'N'},
	'Z': {'E', 'Q'},
	'J': {'I', 'L'},
	'+': {'K', 'R', 'H'}, // Positively charged
	'-': {'D', 'E'},      // Negatively charged
}

func IUPAC2Bases() map[rune][]rune {
	return iupacToBases
}

func IUPAC2AAs() map[rune][]rune {
	return iupacToAAs
}

// AllVariantsIUPAC takes a string as input
// and returns all iupac variants as output
func AllVariantsIUPAC(seq string) ([]string, error) {
	seqVariantList := [][]rune{}
	seqVariants := []string{}

	for _, s := range strings.ToUpper(seq) {
		variantsIUPAC, ok := iupacToBases[s]
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

func cartRune(inLists ...[]rune) [][]rune {
	// An iterative approach to calculate Cartesian product of two or more lists
	// Adapted from https://rosettacode.org/wiki/Cartesian_product_of_two_or_more_lists
	// supposedly "minimizes allocations and computes and fills the result sequentially"

	possibleVariants := 1 // a counter used to determine the possible number of variants
	for _, inList := range inLists {
		possibleVariants *= len(inList)
	}
	if possibleVariants == 0 {
		return nil // in the future this could be part of an error return?
	}
	allVariants := make([][]rune, possibleVariants)               // this is the 2D slice where all variants will be stored
	variantHolders := make([]rune, possibleVariants*len(inLists)) // this is an empty slice with a length totaling the size of all input characters
	variantChoices := make([]int, len(inLists))                   // these will be all the possible variants
	start := 0
	for variant := range allVariants {
		end := start + len(inLists) // define end point
		variantHolder := variantHolders[start:end]

		allVariants[variant] = variantHolder

		start = end // start at end point

		for variantChoicesIndex, variantChoice := range variantChoices {
			variantHolder[variantChoicesIndex] = inLists[variantChoicesIndex][variantChoice]
		}
		for variantChoicesIndex := len(variantChoices) - 1; variantChoicesIndex >= 0; variantChoicesIndex-- {
			variantChoices[variantChoicesIndex]++
			if variantChoices[variantChoicesIndex] < len(inLists[variantChoicesIndex]) {
				break
			}
			variantChoices[variantChoicesIndex] = 0
		}
	}
	return allVariants
}
