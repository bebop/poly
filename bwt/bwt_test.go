package bwt

import (
	"strings"
	"testing"

	"golang.org/x/exp/slices"
)

type BWTCountTestCase struct {
	seq      string
	expected int
}

func TestBWT_Count(t *testing.T) {
	baseTestStr := "thequickbrownfoxjumpsoverthelazydogwithanovertfrownafterfumblingitsparallelogramshapedbananagramallarounddowntown"
	testStr := strings.Join([]string{baseTestStr, baseTestStr, baseTestStr}, "")

	bwt := New(testStr)

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

type BWTLocateTestCase struct {
	seq      string
	expected []int
}

func TestBWT_Locate(t *testing.T) {
	baseTestStr := "thequickbrownfoxjumpsoverthelazydogwithanovertfrownafterfumblingitsparallelogramshapedbananagramallarounddowntown" // len == 112
	testStr := strings.Join([]string{baseTestStr, baseTestStr, baseTestStr}, "")

	bwt := New(testStr)

	testTable := []BWTLocateTestCase{
		{"uick", []int{4, 117, 230}},
		{"the", []int{0, 25, 113, 138, 226, 251}},
		{"over", []int{21, 41, 134, 154, 247, 267}},
		{"own", []int{10, 48, 106, 110, 123, 161, 219, 223, 236, 274, 332, 336}},
		{"ana", []int{87, 89, 200, 202, 313, 315}},
		{"an", []int{39, 87, 89, 152, 200, 202, 265, 313, 315}},
		{"na", []int{50, 88, 90, 163, 201, 203, 276, 314, 316}},
		{"rown", []int{9, 47, 122, 160, 235, 273}},
		{"townthe", []int{109, 222}},
		{"zzz", nil},
	}

	for _, v := range testTable {
		offsets := bwt.Locate(v.seq)
		slices.Sort(offsets)
		if len(offsets) != len(v.expected) {
			t.Fatalf("seq=%s expectedOffsets=%v actualOffsets=%v", v.seq, v.expected, offsets)
		}
		for i := range offsets {
			if offsets[i] != v.expected[i] {
				t.Fatalf("seq=%s expectedOffsets=%v actualOffsets=%v", v.seq, v.expected, offsets)
			}
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
