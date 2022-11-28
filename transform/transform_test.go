package transform

import (
	"math/rand"
	"testing"
)

func TestReverse(t *testing.T) {
	rng := rand.New(rand.NewSource(1))
	sequence := randomSequence(rng)
	// Test even and odd lengthed string.
	for _, test := range []string{sequence, sequence[:len(sequence)-1]} {
		got := Reverse(test)
		for i := range got[:len(got)/2+1] {
			gotbase := got[i]
			expect := test[len(test)-i-1]
			if gotbase != expect {
				t.Errorf("mismatch at pos %d, got %q, expect %q", i, gotbase, expect)
			}
		}
	}
}

func TestComplement(t *testing.T) {
	rng := rand.New(rand.NewSource(1))
	original := randomSequence(rng)
	for _, test := range []string{original} {
		gotSequence := Complement(test)
		for i, got := range gotSequence {
			expect := ComplementBase(rune(test[i]))
			if got != expect {
				t.Errorf("bad %q complement: got %q, expect %q", test[i], got, expect)
			}
		}
	}
}

func TestReverseComplement(t *testing.T) {
	rng := rand.New(rand.NewSource(1))
	original := randomSequence(rng)
	for _, test := range []string{original, original[:len(original)-1]} {
		got := ReverseComplement(test)
		expect := Reverse(Complement(test))
		if got != expect {
			t.Errorf("mismatch with individual Reverse and Complement call:\n%q\n%q", got, expect)
		}
	}
}

func TestComplementBase(t *testing.T) {
	for _, c := range letters {
		got := ComplementBase(rune(c))
		gotI := ComplementBase(got)
		gotII := ComplementBase(gotI)
		if (c == 'U' && got == 'A') || (c == 'u' && got == 'a') {
			continue // Edge case: RNA Uracil base-pairs with Adenine.
		}
		if rune(c) != gotI || gotII != got {
			t.Errorf("complement transform mismatch: %q->%q->%q->%q", c, got, gotI, gotII)
		}
	}
}

func randomSequence(rnd *rand.Rand) string {
	// TODO(soypat): make a more professional random sequencer
	// Probably could have a length argument.
	lettersCopy := letters
	rnd.Shuffle(len(lettersCopy), func(i, j int) { lettersCopy[i], lettersCopy[j] = lettersCopy[j], lettersCopy[i] })
	return string(lettersCopy[:])
}

var letters = [...]byte{'A', 'B', 'C', 'D', 'G', 'H', 'K', 'M', 'N', 'R', 'S', 'T', 'U', 'V', 'W', 'Y', 'a', 'b', 'c', 'd', 'g', 'h', 'k', 'm', 'n', 'r', 's', 't', 'u', 'v', 'w', 'y'}

func BenchmarkReverseComplement(b *testing.B) {
	rng := rand.New(rand.NewSource(1))
	seq := randomSequence(rng)
	for i := 0; i < b.N; i++ {
		seq = ReverseComplement(seq)
	}
	_ = seq
}

func BenchmarkComplement(b *testing.B) {
	rng := rand.New(rand.NewSource(1))
	seq := randomSequence(rng)
	for i := 0; i < b.N; i++ {
		seq = Complement(seq)
	}
	_ = seq
}

func BenchmarkReverse(b *testing.B) {
	rng := rand.New(rand.NewSource(1))
	seq := randomSequence(rng)
	for i := 0; i < b.N; i++ {
		seq = Reverse(seq)
	}
	_ = seq
}
