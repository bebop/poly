package checks

import "github.com/TimothyStiles/poly/transformations"

// IsPalindromic accepts a sequence of even length and returns if it is
// palindromic. More here - https://en.wikipedia.org/wiki/Palindromic_sequence
func IsPalindromic(sequence string) bool {
	return sequence == transformations.ReverseComplement(sequence)
}
