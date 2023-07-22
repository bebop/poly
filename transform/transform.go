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
	sequenceLength := len(sequence)
	newSequence := make([]byte, sequenceLength)
	for index := 0; index < sequenceLength; index++ {
		newSequence[index] = complementTable[sequence[sequenceLength-index-1]]
	}
	// This is how strings.Builder works with the String() method. If Mr. Go says it's safe...
	return *(*string)(unsafe.Pointer(&newSequence))
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
	sequenceLength := len(sequence)
	newSequence := make([]byte, sequenceLength)
	for index := 0; index < sequenceLength; index++ {
		newSequence[index] = complementTable[sequence[index]]
	}
	// This is how strings.Builder works with the String() method. If Mr. Go says it's safe...
	return *(*string)(unsafe.Pointer(&newSequence))
}

// Reverse returns the reverse of sequence. It performs a basic
// string reversal by working on the bytes.
//
// This function expects byte characters in the range a-z and A-Z and
// will not check for non-byte character, i.e. utf-8 encoding.
func Reverse(sequence string) string {
	sequenceLength := len(sequence)
	newSequence := make([]byte, sequenceLength)
	for index := 0; index < sequenceLength; index++ {
		newSequence[index] = sequence[sequenceLength-index-1]
	}
	// This is how strings.Builder works with the String() method. If Mr. Go says it's safe...
	return *(*string)(unsafe.Pointer(&newSequence))
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
	'v': 'b',
	'w': 'w',
	'y': 'r',
}

// ReverseComplementRNA returns the reversed complement of sequence.
// It is the equivalent of calling
//
//	revComplement := Reverse(ComplementRNA(sequence))
//
// This function expects byte characters in the range a-z and A-Z and
// will not check for non-byte characters, i.e. utf-8 encoding.
func ReverseComplementRNA(sequence string) string {
	sequenceLength := len(sequence)
	newSequence := make([]byte, sequenceLength)
	for index := 0; index < sequenceLength; index++ {
		newSequence[index] = complementTableRNA[sequence[sequenceLength-index-1]]
	}
	// This is how strings.Builder works with the String() method. If Mr. Go says it's safe...
	return *(*string)(unsafe.Pointer(&newSequence))
}

// ComplementRNA returns the complement of sequence. In [RNA] each nucleotide
// (A, U, C or G) is has a deterministic pair. A is paired with U and
// C is paired with G. The complement of a sequence is formed by exchanging
// these letters for their pair:
//
//	'A' becomes 'U'
//	'U' becomes 'A'
//	'C' becomes 'G'
//	'G' becomes 'C'
//
// This function expects byte characters in the range a-z and A-Z and
// will not check for non-byte characters, i.e. utf-8 encoding.
//
// [RNA]: https://en.wikipedia.org/wiki/RNA
func ComplementRNA(sequence string) string {
	sequenceLength := len(sequence)
	newSequence := make([]byte, sequenceLength)
	for index := 0; index < sequenceLength; index++ {
		newSequence[index] = complementTableRNA[sequence[index]]
	}
	// This is how strings.Builder works with the String() method. If Mr. Go says it's safe...
	return *(*string)(unsafe.Pointer(&newSequence))
}

// ComplementBaseRNA accepts a RNA base pair and returns its complement base
// pair. See Complement.
//
// This function expects byte characters in the range a-z and A-Z and
// will return a space ' ' (U+0020) for characters that are not matched
// to any known base. This is subject to change.
func ComplementBaseRNA(basePair rune) rune {
	got := rune(complementTableRNA[basePair])
	if got == 0 {
		return ' ' // invalid sequence returns empty space.
	}
	return got
}

// complementTable provides 1:1 mapping between bases and their complements
// see https://www.dnabaser.com/articles/IUPAC%20ambiguity%20codes.html
var complementTableRNA = [256]byte{
	'A': 'U',
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
	'U': 'A',
	'V': 'B',
	'W': 'W',
	'Y': 'R',
	'X': 'X',
	'a': 'u',
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
	'u': 'a',
	'v': 'b',
	'w': 'w',
	'y': 'r',
	'x': 'x',
}
