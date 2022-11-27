package transform

import (
	"math/rand"
	"testing"
)

// goos: linux
// goarch: amd64
// pkg: github.com/TimothyStiles/poly/transform
// cpu: Intel(R) Core(TM) i5-8265U CPU @ 1.60GHz
// BenchmarkReverseComplement-8   	 1804766	       669.8 ns/op	     224 B/op	       3 allocs/op
// PASS
func BenchmarkReverseComplement(b *testing.B) {
	rng := rand.New(rand.NewSource(1))
	seq := randomSequence(rng)
	for i := 0; i < b.N; i++ {
		seq = ReverseComplement(seq)
	}
	_ = seq
}

func TestComplementBase(t *testing.T) {
	for _, c := range letters {
		got := ComplementBase(rune(c))
		gotI := ComplementBase(got)
		gotII := ComplementBase(gotI)
		if rune(c) != gotI {
			t.Errorf("double complement should yield start: %q->%q->%q", c, got, gotI)
		}
		if gotII != got {
			t.Errorf("double complement should yield start: %q->%q->%q", got, gotI, gotII)
		}
	}
}

func randomSequence(rnd *rand.Rand) string {
	cp := letters
	rnd.Shuffle(len(cp), func(i, j int) { cp[i], cp[j] = cp[j], cp[i] })
	return string(cp[:])
}

var letters = [...]byte{'A', 'B', 'C', 'D', 'G', 'H', 'K', 'M', 'N', 'R', 'S', 'T', 'U', 'V', 'W', 'Y', 'a', 'b', 'c', 'd', 'g', 'h', 'k', 'm', 'n', 'r', 's', 't', 'u', 'v', 'w', 'y'}
