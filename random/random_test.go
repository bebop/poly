package random

import (
	"fmt"
	"strconv"
	"testing"
)

func ExampleProteinSequence() {
	// RandomProteinSequence builds a Protein Sequence by only passing through arguments a length and a seed that will be use to generate a randomly the sequence. The length needs to be greater than two because every sequence already have a start and stop codon. Seed makes this test deterministic.
	randomProtein, _ := ProteinSequence(15, 2)
	fmt.Println(randomProtein)

	// Output: MHHPAFRMFNTMYG*
}

func TestRandomProteinSequence(t *testing.T) {
	const length = 10
	const seed = 2
	sequence, _ := ProteinSequence(length, seed)

	if sequence[0] != 'M' {
		t.Errorf("Random sequence doesn't have the correct initial aminoacid in sequence 'RandomSequence(10, 2)'. Got this: \n%s instead of \n%s", string(sequence[0]), "M")
	}

	if sequence[len(sequence)-1] != '*' {
		t.Errorf("Random sequence doesn't have correct last aminoacid in sequence 'RandomSequence(10, 2)'. Got this: \n%s instead of \n%s", string(sequence[len(sequence)-1]), "*")
	}

	if len(sequence) != length {
		t.Errorf("Random sequence doesn't have the sequence size equal parameter passed through n 'RandomSequence(10, 2)'. Got this: \n%s instead of \n%s", strconv.Itoa(len(sequence)), strconv.Itoa(length))
	}
}

// Write a new case of test when you have a n inferior than 3
func TestRandomProteinSequenceError(t *testing.T) {
	const length = 2
	const seed = 4
	sequence, _ := ProteinSequence(length, seed)

	if len(sequence) != 0 {
		t.Errorf("Random sequence must have sequence size equals 0 'RandomSequence(2, 4)'. Got this: \n%s instead of \n%s", strconv.Itoa(len(sequence)), strconv.Itoa(length))
	}
}
func ExampleDNASequence() {
	// RandomDNASequence builds a DNA Sequence by only passing through arguments a length and a seed that will be use to generate a randomly the sequence. The length needs to be greater than two because every sequence already have a start and stop codon. Seed makes this test deterministic.
	randomDNA, _ := DNASequence(15, 2)
	fmt.Println(randomDNA)

	// Output: TTAAATTAGATGCAA
}
