/*
Package checks provides utilities to check for certain properties of a sequence.
*/
package checks

import (
	"strings"

	"github.com/TimothyStiles/poly/transform"
)

// IsPalindromic accepts a sequence of even length and returns if it is
// palindromic. More here - https://en.wikipedia.org/wiki/Palindromic_sequence
func IsPalindromic(sequence string) bool {
	return sequence == transform.ReverseComplement(sequence)
}

// GcContent checks the GcContent of a given sequence.
func GcContent(sequence string) float64 {
	sequence = strings.ToUpper(sequence)
	GuanineCount := strings.Count(sequence, "G")
	CytosineCount := strings.Count(sequence, "C")
	GuanineAndCytosinePercentage := float64(GuanineCount+CytosineCount) / float64(len(sequence))
	return GuanineAndCytosinePercentage
}

func IsDNA(seq string) bool {
	for _, base := range seq {
		switch base {
		case 'A', 'C', 'T', 'G':
			continue
		default:
			return false
		}
	}
	return true
}

func IsRNA(seq string) bool {
	for _, base := range seq {
		switch base {
		case 'A', 'C', 'U', 'G':
			continue
		default:
			return false
		}
	}
	return true
}
