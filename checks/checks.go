package checks

import (
	"fmt"
	"regexp"

	"github.com/TimothyStiles/poly/transform"
)

// IsPalindromic accepts a sequence of even length and returns if it is
// palindromic. More here - https://en.wikipedia.org/wiki/Palindromic_sequence
func IsPalindromic(sequence string) bool {
	return sequence == transform.ReverseComplement(sequence)
}

// IsValidDotBracketStructure accepts a string and checks if it uses valid
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
