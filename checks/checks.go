package checks

import (
	"fmt"
	"regexp"

	"github.com/Open-Science-Global/poly/transform"
)

// IsPalindromic accepts a sequence of even length and returns if it is
// palindromic. More here - https://en.wikipedia.org/wiki/Palindromic_sequence
func IsPalindromic(sequence string) bool {
	return sequence == transform.ReverseComplement(sequence)
}

// GC content calculates the percentage of G and C base pairs content in a given sequence
func GcContent(sequence string) float64 {

	if len(sequence) < 1 {
	  return 0.0
	}
  
	const c = 67
	const g = 71
	count := 0.0
  
	for _, character := range(sequence) {
	  if character == g || character == c {
		count += 1.0 
	  }
	}
  
	return count / float64(len(sequence))
}

// IsValidDotBracketStructure accepts a string and checks if uses valid
// dot-bracket notation. See the `secondary_structure` package for more info
// on dot-bracket notation.
func IsValidDotBracketStructure(structure string) (bool, error) {
	dotBracketRegex := "^[().]+"
	return checkRegexpMatchesFullString(structure, dotBracketRegex, "found invalid characters in structure. Only dot-bracket notation allowed")
}

// IsValidRNA accepts a string and checks if it is a valid RNA sequence.
func IsValidRNA(sequence string) (bool, error) {
	rnaRegex := "^[ACGU]+"
	return checkRegexpMatchesFullString(sequence, rnaRegex, "expected valid RNA sequence (A, C, G, or U only).")
}

func checkRegexpMatchesFullString(str, regex, errMsg string) (bool, error) {
	regexp, err := regexp.Compile(regex)
	if err != nil {
		return false, fmt.Errorf("regexp.Compile() error: %s", err)
	}

	if !doCheckRegexpMatchesFullString(str, regexp) {
		return false, fmt.Errorf(errMsg)
	}
	return true, nil
}

func doCheckRegexpMatchesFullString(str string, regexp *regexp.Regexp) bool {
	strIdxsThatMatchRegexp := regexp.FindStringIndex(str)

	return len(strIdxsThatMatchRegexp) != 0 && strIdxsThatMatchRegexp[0] == 0 && strIdxsThatMatchRegexp[1] == len(str)
}
