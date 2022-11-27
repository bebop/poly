package transform_test

import (
	"fmt"
	"testing"

	"github.com/TimothyStiles/poly/transform"
)

func ExampleReverseComplement() {
	sequence := "GATTACA"
	reverseComplement := transform.ReverseComplement(sequence)
	fmt.Println(reverseComplement)

	// Output: TGTAATC
}

func ExampleComplement() {
	sequence := "GATTACA"
	complement := transform.Complement(sequence)
	fmt.Println(complement)

	// Output: CTAATGT
}

func ExampleReverse() {
	sequence := "GATTACA"
	reverse := transform.Reverse(sequence)
	fmt.Println(reverse)

	// Output: ACATTAG
}
func TestA(t *testing.T) {
	for k, v := range complementBaseRuneMap2 {
		got := transform.ComplementBase(k)
		if v != got {
			t.Errorf("%q: %q %q", k, v, got)
		}
		gotInverse := transform.ComplementBase(got)
		if gotInverse != k {
			t.Errorf("%q: %q %q", got, k, gotInverse)
		}
	}
}

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
