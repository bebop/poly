package transform

import (
	"fmt"
	"strconv"
	"testing"
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

func exampleKmerTable(t *testing.T) {
	var kmers []string
	kmerTable := GetKmerTable(2, "ATCG")

	for k := range kmerTable {
		kmers = append(kmers, k)
	}

	if len(kmers) == 3 && kmers[0] == "AT" && kmers[1] == "TC" && kmers[2] == "CG" {
		t.Errorf("Expected kmers table to have length 3, got: %s", strconv.Itoa(len(kmers)))
		t.Errorf("Expected kmers table to have first Element 1 AT, got: %s", kmers[0])
		t.Errorf("Expected kmers table to have first Element 2 TC, got: %s", kmers[1])
		t.Errorf("Expected kmers table to have first Element 3 CG, got: %s", kmers[2])
	}
}
