package bwt

import (
	"strings"
	"testing"
)

type BWTCountTestCase struct {
	seq      string
	expected int
}

const augmentedQuickBrownFoxTest = "thequickbrownfoxjumpsoverthelazydogwithanovertfrownafterfumblingitsparallelogramshapedbananagramallarounddowntown"

var threeAugmentedQuickBrownFoxTest = strings.Join([]string{augmentedQuickBrownFoxTest, augmentedQuickBrownFoxTest, augmentedQuickBrownFoxTest}, "")

func TestBWT_Count(t *testing.T) {
	bwt := New(threeAugmentedQuickBrownFoxTest)

	testTable := []BWTCountTestCase{
		{"uick", 3},
		{"the", 6},
		{"over", 6},
		{"own", 12},
		{"ana", 6},
		{"an", 9},
		{"na", 9},
		{"rown", 6},
		{"townthe", 2},
		{"zzz", 0},
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
