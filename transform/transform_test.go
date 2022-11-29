package transform

import (
	"math/rand"
	"testing"

	"github.com/TimothyStiles/poly/random"
)

func TestReverse(t *testing.T) {
	seed := rand.New(rand.NewSource(1)).Int63()
	sequence, err := random.DNASequence(20, seed)
	if err != nil {
		t.Errorf("failed to generate random DNA sequence")
	}
	// Test even and odd lengthed string.
	evenAndOddLengthSequences := []string{sequence, sequence[:len(sequence)-1]} // indexing just removes last character of sequence for even/odd check
	for _, testSequence := range evenAndOddLengthSequences {
		reversedSequence := Reverse(testSequence)
		for index := range reversedSequence[:len(reversedSequence)/2+1] {
			gotbase := reversedSequence[index]
			expect := testSequence[len(testSequence)-index-1]
			if gotbase != expect {
				t.Errorf("mismatch at pos %d, got %q, expect %q", index, gotbase, expect)
			}
		}
	}
}

func TestComplement(t *testing.T) {
	seed := rand.New(rand.NewSource(1)).Int63()
	sequence, err := random.DNASequence(20, seed)
	if err != nil {
		t.Errorf("failed to generate random DNA sequence")
	}
	for _, testSequence := range []string{sequence} { // loop is unneeded but makes it easier to add more test cases in the future.
		complementSequence := Complement(testSequence)
		for index, complement := range complementSequence {
			expectedBase := ComplementBase(rune(testSequence[index]))
			if complement != expectedBase {
				t.Errorf("bad %q complement: got %q, expect %q", testSequence[index], complement, expectedBase)
			}
		}
	}
}

func TestReverseComplement(t *testing.T) {
	seed := rand.New(rand.NewSource(1)).Int63()
	sequence, err := random.DNASequence(20, seed)
	if err != nil {
		t.Errorf("failed to generate random DNA sequence")
	}
	evenAndOddLengthSequences := []string{sequence, sequence[:len(sequence)-1]}
	for _, testSequence := range evenAndOddLengthSequences {
		got := ReverseComplement(testSequence)
		expect := Reverse(Complement(testSequence))
		if got != expect {
			t.Errorf("mismatch with individual Reverse and Complement call:\n%q\n%q", got, expect)
		}
	}
}

func TestComplementBaseError(t *testing.T) {
	complementBase := ComplementBase('!')
	if complementBase != ' ' {
		t.Errorf("expected space for invalid base")
	}
}

func BenchmarkReverseComplement(b *testing.B) {
	seed := rand.New(rand.NewSource(1)).Int63()
	sequence, _ := random.DNASequence(8, seed)
	for index := 0; index < b.N; index++ {
		sequence = ReverseComplement(sequence)
	}
	_ = sequence
}

func BenchmarkComplement(b *testing.B) {
	seed := rand.New(rand.NewSource(1)).Int63()
	sequence, _ := random.DNASequence(8, seed)
	for index := 0; index < b.N; index++ {
		sequence = Complement(sequence)
	}
	_ = sequence
}

func BenchmarkReverse(b *testing.B) {
	seed := rand.New(rand.NewSource(1)).Int63()
	sequence, _ := random.DNASequence(8, seed)
	for index := 0; index < b.N; index++ {
		sequence = Reverse(sequence)
	}
	_ = sequence
}
