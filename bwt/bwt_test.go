package bwt

import (
	"testing"
)

type QueryTest struct {
	seq      string
	expected bool
}

func TestQueryBWT(t *testing.T) {
	bwt := New("BANANA")

	testTable := []QueryTest{
		{"NANA", true},
		{"ANA", true},
		{"NA", true},
		{"B", true},
		{"N", true},
		{"BA", true},
		{"ANANA", true},
		{"QWERTY", false},
		{"ANANANA", false},
		{"ABCD", false},
		{"ABA", false},
	}

	for _, v := range testTable {
		res := bwt.QueryExistence(v.seq)
		if res != v.expected {
			t.Fatalf("Test=%s ExpectedQueryExistence=%v Received=%v", v.seq, v.expected, res)
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
		bwt.QueryExistence(seq)
	}
}
