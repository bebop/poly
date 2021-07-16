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

// IsValidDotBracketStructure accepts a string and checks if uses valid
// dot-bracket notation. See the `secondary_structure` package for more info
// on dot-bracket notation.
func IsValidDotBracketStructure(structure string) (bool, error) {
	dotBracketStructureRegex, _ := regexp.Compile("^[().]+")

	if !checkRegexpMatchesFullString(structure, dotBracketStructureRegex) {
		return false, fmt.Errorf("found invalid characters in structure. Only dot-bracket notation allowed")
	}
	return true, nil
}

func checkRegexpMatchesFullString(str string, regexp *regexp.Regexp) bool {
	strIdxsThatMatchRegexp := regexp.FindStringIndex(str)

	return len(strIdxsThatMatchRegexp) != 0 && strIdxsThatMatchRegexp[0] == 0 && strIdxsThatMatchRegexp[1] == len(str)
}
