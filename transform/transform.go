/*
Package transform provides functions for transforming sequences.

Complement takes the complement of a sequence.
(returns a sequence string where each nucleotide has been swapped with its complement A<->T, C<->G)

Reverse takes the reverse of a sequence.
(literally just reverses a string. Exists in stdlib but hey why not have it here too?)

ReverseComplement takes the reverse complement of a sequence.
(Reverses the sequence string and returns the complement of the reversed sequence.)
*/
package transform

import "unsafe"

// ReverseComplement takes the reverse complement of a sequence.
// ReverseComplement expects ASCII input.
func ReverseComplement(sequence string) string {
	n := len(sequence)
	newSeq := make([]byte, n)
	for i := 0; i < n; i++ {
		newSeq[i] = complementTable[sequence[n-i-1]]
	}
	// This is how strings.Builder works with the String() method. If Mr. Go says it's safe...
	return *(*string)(unsafe.Pointer(&newSeq))
}

// Complement takes the complement of a sequence.
func Complement(sequence string) string {
	n := len(sequence)
	newSeq := make([]byte, n)
	for i := 0; i < n; i++ {
		newSeq[i] = complementTable[sequence[i]]
	}
	// This is how strings.Builder works with the String() method. If Mr. Go says it's safe...
	return *(*string)(unsafe.Pointer(&newSeq))
}

// Reverse takes the reverse of a sequence.
func Reverse(sequence string) string {
	n := len(sequence)
	newSeq := make([]byte, n)
	for i := 0; i < n; i++ {
		newSeq[i] = sequence[n-i-1]
	}
	// This is how strings.Builder works with the String() method. If Mr. Go says it's safe...
	return *(*string)(unsafe.Pointer(&newSeq))
}

// ComplementBase accepts a base pair and returns its complement base pair
func ComplementBase(basePair rune) rune {
	got := rune(complementTable[basePair])
	if got == 0 {
		return ' ' // invalid sequence returns empty space.
	}
	return got
}

// complementTable provides 1:1 mapping between bases and their complements
var complementTable = [256]byte{
	'A': 'T',
	'B': 'V',
	'C': 'G',
	'D': 'H',
	'G': 'C',
	'H': 'D',
	'K': 'M',
	'M': 'K',
	'N': 'N',
	'R': 'Y',
	'S': 'S',
	'T': 'A',
	'U': 'A',
	'V': 'B',
	'W': 'W',
	'Y': 'R',
	'a': 't',
	'b': 'v',
	'c': 'g',
	'd': 'h',
	'g': 'c',
	'h': 'd',
	'k': 'm',
	'm': 'k',
	'n': 'n',
	'r': 'y',
	's': 's',
	't': 'a',
	'u': 'a',
	'v': 'b',
	'w': 'w',
	'y': 'r',
}
