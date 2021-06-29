package checks

import "github.com/TimothyStiles/poly/transform/reverse"

// IsPalindromic accepts a sequence of even length and returns if it is
// palindromic. More here - https://en.wikipedia.org/wiki/Palindromic_sequence
func IsPalindromic(sequence string) bool {
	return sequence == reverse.ReverseComplement(sequence)
}
