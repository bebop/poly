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

import "strings"

// complementBaseRuneMap provides 1:1 mapping between bases and their complements
var complementBaseRuneMap = map[rune]rune{
	65:  84,  // A -> T
	66:  86,  // B -> V
	67:  71,  // C -> G
	68:  72,  // D -> H
	71:  67,  // G -> C
	72:  68,  // H -> D
	75:  77,  // K -> M
	77:  75,  // M -> K
	78:  78,  // N -> N
	82:  89,  // R -> Y
	83:  83,  // S -> S
	84:  65,  // T -> A
	85:  65,  // U -> A
	86:  66,  // V -> B
	87:  87,  // W -> W
	89:  82,  // Y -> R
	97:  116, // a -> t
	98:  118, // b -> v
	99:  103, // c -> g
	100: 104, // d -> h
	103: 99,  // g -> a mistake?
	104: 100, // h -> d
	107: 109, // k -> m
	109: 107, // m -> k
	110: 110, // n -> n
	114: 121, // r -> y
	115: 115, // s -> s
	116: 97,  // t -> a
	117: 97,  // u -> a
	118: 98,  // v -> b
	119: 119, // w -> w
	121: 114, // y -> r
}

const aa = 'a'

var complementBaseRuneMap2 = map[rune]rune{
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
	'g': 'a',
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

// ReverseComplement takes the reverse complement of a sequence.
func ReverseComplement(sequence string) string {
	complementString := strings.Map(ComplementBase, sequence)
	length := len(complementString)
	newString := make([]rune, length)
	for _, base := range complementString {
		length--
		newString[length] = base
	}
	return string(newString)
}

// Complement takes the complement of a sequence.
func Complement(sequence string) string {
	complementString := strings.Map(ComplementBase, sequence)
	return complementString
}

// Reverse takes the reverse of a sequence.
func Reverse(sequence string) string {
	length := len(sequence)
	newString := make([]rune, length)
	for _, base := range sequence {
		length--
		newString[length] = base
	}
	return string(newString)
}

// ComplementBase accepts a base pair and returns its complement base pair
func ComplementBase(basePair rune) rune {
	return complementBaseRuneMap[basePair]
}
