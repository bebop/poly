package transform

import (
	"fmt"
)

func ExampleReverseComplement() {
	sequence := "GATTACA"
	reverseComplement := ReverseComplement(sequence)
	fmt.Println(reverseComplement)

	// Output: TGTAATC
}

func ExampleComplement() {
	sequence := "GATTACA"
	complement := Complement(sequence)
	fmt.Println(complement)

	// Output: CTAATGT
}

func ExampleReverse() {
	sequence := "GATTACA"
	reverse := Reverse(sequence)
	fmt.Println(reverse)

	// Output: ACATTAG
}

func ExampleGetKmerTable() {
	kmerTable := GetKmerTable(2, "ATCG")

	fmt.Println(kmerTable)
	//Output: map[AT:true CG:true TC:true]
}
