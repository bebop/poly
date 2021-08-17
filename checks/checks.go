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

func GcContent(sequence string) float64 {
	sequence = strings.ToUpper(sequence)
	return float64(strings.Count(sequence, "G")+strings.Count(sequence, "C")) / float64(len(sequence))
}
