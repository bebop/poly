package transform

import (
	"math/rand"
	"testing"

	"github.com/bebop/poly/random"
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

func TestComplementBase(t *testing.T) {
	var letters = [...]byte{'A', 'B', 'C', 'D', 'G', 'H', 'K', 'M', 'N', 'R', 'S', 'T', 'V', 'W', 'Y', 'a', 'b', 'c', 'd', 'g', 'h', 'k', 'm', 'n', 'r', 's', 't', 'v', 'w', 'y'}
	for _, c := range letters {
		got := ComplementBase(rune(c))
		gotI := ComplementBase(got)
		gotII := ComplementBase(gotI)
		if rune(c) != gotI || gotII != got {
			t.Errorf("complement transform mismatch: %q->%q->%q->%q", c, got, gotI, gotII)
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

func TestComplementBaseRNA(t *testing.T) {
	bases := [...]rune{'A', 'B', 'C', 'D', 'G', 'H', 'K', 'M', 'N', 'R', 'S', 'U', 'V', 'W', 'Y', 'a', 'b', 'c', 'd', 'g', 'h', 'k', 'm', 'n', 'r', 's', 'u', 'v', 'w', 'y'}
	for _, c := range bases {
		got := ComplementBaseRNA(c)
		gotI := ComplementBaseRNA(got)
		gotII := ComplementBaseRNA(gotI)
		if c != gotI || gotII != got {
			t.Errorf("complement transform mismatch: %q->%q->%q->%q", c, got, gotI, gotII)
		}
	}
}

func TestComplementRNANoRandom(t *testing.T) {
	t.Skip()
	// keys of complementTableRNA plus '!'
	sequence := string([]rune{'A', 'B', 'C', 'D', 'G', 'H', 'K', 'M', 'N', 'R', 'S', 'U', 'V', 'W', 'Y', 'a', 'b', 'c', 'd', 'g', 'h', 'k', 'm', 'n', 'r', 's', 'u', 'v', 'w', 'y', '!'})

	// values of complementTableRNA plus ' '
	complement := string([]rune{'U', 'V', 'G', 'H', 'C', 'D', 'M', 'K', 'N', 'Y', 'S', 'A', 'B', 'W', 'R', 'u', 'v', 'g', 'h', 'c', 'd', 'm', 'k', 'n', 'y', 's', 'a', 'b', 'w', 'r', ' '})

	result := ComplementRNA(sequence)
	if result != complement {
		t.Errorf("bad complement: got %q, expect %q", result, complement)
	}
}

func TestComplementRNA(t *testing.T) {
	seed := rand.New(rand.NewSource(1)).Int63()
	sequence, err := random.RNASequence(20, seed)
	if err != nil {
		t.Errorf("failed to generate random DNA sequence")
	}
	for _, testSequence := range []string{sequence} { // loop is unneeded but makes it easier to add more test cases in the future.
		complementSequence := ComplementRNA(testSequence)
		for index, complement := range complementSequence {
			expectedBase := ComplementBaseRNA(rune(testSequence[index]))
			if complement != expectedBase {
				t.Errorf("bad %q complement: got %q, expect %q", testSequence[index], complement, expectedBase)
			}
		}
	}
}

func TestReverseComplementRNA(t *testing.T) {
	seed := rand.New(rand.NewSource(1)).Int63()
	sequence, err := random.RNASequence(20, seed)
	if err != nil {
		t.Errorf("failed to generate random RNA sequence")
	}
	evenAndOddLengthSequences := []string{sequence, sequence[:len(sequence)-1]}
	for _, testSequence := range evenAndOddLengthSequences {
		got := ReverseComplementRNA(testSequence)
		expect := Reverse(ComplementRNA(testSequence))
		if got != expect {
			t.Errorf("mismatch with individual Reverse and Complement call:\n%q\n%q", got, expect)
		}
	}
}
