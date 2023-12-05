package bwt

import (
	"testing"
)

type BWTCountTestCase struct {
	seq      string
	expected int
}

func TestBWT_Count(t *testing.T) {
	bwt := New("BANANA")

	testTable := []BWTCountTestCase{
		{"NANA", 1},
		{"ANA", 2},
		{"NA", 2},
		{"B", 1},
		{"N", 2},
		{"BA", 1},
		{"ANANA", 1},
		{"QWERTY", 0},
		{"ANANANA", 0},
		{"ABCD", 0},
		{"ABA", 0},
	}

	for _, v := range testTable {
		count := bwt.Count(v.seq)
		if count != v.expected {
			t.Fatalf("seq=%s expectedCount=%v actualCount=%v", v.seq, v.expected, count)
		}
	}
}

func BenchmarkBWTBuildPower12(b *testing.B) {
	base := "!BANANA!"
	BaseBenchmarkBWTBuild(base, 12, b)
}

//go:noinline
func BaseBenchmarkBWTBuild(base string, power int, b *testing.B) {
	for n := 0; n < b.N; n++ {
		buildBWTForBench(base, power)
	}
}

func buildBWTForBench(base string, power int) BWT {
	test := base
	for i := 0; i < power; i++ {
		test += test
	}

	return New(test)
}

func BenchmarkBWTQueryPower12(b *testing.B) {
	base := "!BANANA!"
	bwt := buildBWTForBench(base, 12)
	BaseBenchmarkBWTQuery(bwt, "ANANABANANA", b)
}

//go:noinline
func BaseBenchmarkBWTQuery(bwt BWT, seq string, b *testing.B) {
	for n := 0; n < b.N; n++ {
		bwt.Count(seq)
	}
}
