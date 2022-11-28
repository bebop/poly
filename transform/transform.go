/*
Package transform provides functions for transforming sequences.
*/
package transform

import "unsafe"

// ReverseComplement returns the reversed complement of sequence.
// It is the equivalent of calling
//
//	revComplement := Reverse(Complement(sequence))
//
// This function expects byte characters in the range a-z and A-Z and
// will not check for non-byte characters, i.e. utf-8 encoding.
func ReverseComplement(sequence string) string {
	n := len(sequence)
	newSeq := make([]byte, n)
	for i := 0; i < n; i++ {
		newSeq[i] = complementTable[sequence[n-i-1]]
	}
	// This is how strings.Builder works with the String() method. If Mr. Go says it's safe...
	return *(*string)(unsafe.Pointer(&newSeq))
}

// Complement returns the complement of sequence. In [DNA] each nucleotide
// (A, T, C or G) is has a deterministic pair. A is paired with T and
// C is paired with G. The complement of a sequence is formed by exchanging
// these letters for their pair:
//
//	'A' becomes 'T'
//	'T' becomes 'A'
//	'C' becomes 'G'
//	'G' becomes 'C'
//
// This function expects byte characters in the range a-z and A-Z and
// will not check for non-byte characters, i.e. utf-8 encoding.
//
// [DNA]: https://en.wikipedia.org/wiki/DNA
func Complement(sequence string) string {
	n := len(sequence)
	newSeq := make([]byte, n)
	for i := 0; i < n; i++ {
		newSeq[i] = complementTable[sequence[i]]
	}
	// This is how strings.Builder works with the String() method. If Mr. Go says it's safe...
	return *(*string)(unsafe.Pointer(&newSeq))
}

// Reverse returns the reverse of sequence. It performs a basic
// string reversal by working on the bytes.
//
// This function expects byte characters in the range a-z and A-Z and
// will not check for non-byte character, i.e. utf-8 encoding.
func Reverse(sequence string) string {
	n := len(sequence)
	newSeq := make([]byte, n)
	for i := 0; i < n; i++ {
		newSeq[i] = sequence[n-i-1]
	}
	// This is how strings.Builder works with the String() method. If Mr. Go says it's safe...
	return *(*string)(unsafe.Pointer(&newSeq))
}

// ComplementBase accepts a base pair and returns its complement base pair. See Complement.
//
// This function expects byte characters in the range a-z and A-Z and
// will return a space ' ' (U+0020) for characters that are not matched
// to any known base. This is subject to change.
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
