package bwt

import (
	"strings"
	"testing"
)

type WaveletTreeAccessTestCase struct {
	pos      int
	expected string
}

func TestWaveletTree_Access(t *testing.T) {
	testStr := "AAAACCCCTTTTGGGG" + "ACTG" + "TGCA" + "TTAA" + "CCGG" + "GGGGTTTTCCCCAAAA"
	wt := NewWaveletTreeFromString(testStr)

	testCases := []WaveletTreeAccessTestCase{
		{0, "A"},
		{3, "A"},
		{4, "C"},
		{7, "C"},
		{8, "T"},
		{9, "T"},
		{11, "T"},
		{12, "G"},
		{13, "G"},
		{15, "G"},

		{16, "A"},
		{17, "C"},
		{18, "T"},
		{19, "G"},

		{20, "T"},
		{21, "G"},
		{22, "C"},
		{23, "A"},

		{24, "T"},
		{25, "T"},
		{26, "A"},
		{27, "A"},

		{28, "C"},
		{29, "C"},
		{30, "G"},
		{31, "G"},

		{32, "G"},
		{35, "G"},
		{36, "T"},
		{39, "T"},
		{40, "C"},
		{41, "C"},
		{43, "C"},
		{44, "A"},
		{46, "A"},
		{47, "A"},
	}

	for _, tc := range testCases {
		actual := string(wt.Access(tc.pos))
		if actual != tc.expected {
			t.Fatalf("expected access(%d) to be %s but got %s", tc.pos, tc.expected, actual)
		}
	}
}

type WaveletTreeRankTestCase struct {
	char     string
	pos      int
	expected int
}

func TestWaveletTree_Rank_Genomic(t *testing.T) {
	testStr := "AAAACCCCTTTTGGGG" + "ACTG" + "TGCA" + "TTAA" + "CCGG" + "GGGGTTTTCCCCAAAA"
	wt := NewWaveletTreeFromString(testStr)

	testCases := []WaveletTreeRankTestCase{
		{"A", 0, 0},
		{"A", 2, 2},
		{"A", 3, 3},
		{"A", 8, 4},
		{"C", 4, 0},
		{"C", 6, 2},
		{"C", 12, 4},
		{"T", 2, 0},
		{"T", 8, 0},
		{"T", 12, 4},
		{"T", 15, 4},
		{"G", 15, 3},

		{"A", 16, 4},
		{"A", 17, 5},
		{"G", 16, 4},

		{"T", 20, 5},
		{"A", 23, 5},

		{"T", 24, 6},
		{"T", 27, 8},

		{"C", 28, 6},
		{"G", 31, 7},

		{"G", 32, 8},
		{"G", 33, 9},
		{"T", 36, 8},
		{"T", 38, 10},
		{"C", 40, 8},
		{"C", 43, 11},
		{"A", 44, 8},
		{"A", 47, 11},
	}

	for _, tc := range testCases {
		actual := wt.Rank(tc.char[0], tc.pos)
		if actual != tc.expected {
			t.Fatalf("expected rank(%s, %d) to be %d but got %d", tc.char, tc.pos, tc.expected, actual)
		}
	}
}

type WaveletTreeSelectTestCase struct {
	char     string
	rank     int
	expected int
}

func TestWaveletTree_Select(t *testing.T) {
	testStr := "AAAACCCCTTTTGGGG" + "ACTG" + "TGCA" + "TTAA" + "CCGG" + "GGGGTTTTCCCCAAAA"
	wt := NewWaveletTreeFromString(testStr)

	testCases := []WaveletTreeSelectTestCase{
		{"A", 0, 0},
		{"A", 1, 1},
		{"A", 2, 2},
		{"A", 3, 3},
		{"C", 0, 4},
		{"C", 3, 7},

		{"A", 4, 16},
		{"C", 4, 17},
		{"T", 4, 18},
		{"G", 4, 19},

		{"T", 5, 20},
		{"G", 5, 21},
		{"C", 5, 22},
		{"A", 5, 23},

		{"T", 6, 24},
		{"T", 7, 25},
		{"A", 6, 26},

		{"C", 6, 28},
		{"G", 6, 30},
		{"G", 7, 31},

		{"G", 8, 32},
		{"A", 11, 47},
	}

	for _, tc := range testCases {
		actual := wt.Select(tc.char[0], tc.rank)
		if actual != tc.expected {
			t.Fatalf("expected select(%s, %d) to be %d but got %d", tc.char, tc.rank, tc.expected, actual)
		}
	}
}

func TestWaveletTree_Access_Reconstruction(t *testing.T) {
	enhancedQuickBrownFox := "the quick brown fox jumps over the lazy dog with an overt frown after fumbling its parallelogram shaped bananagram all around downtown"
	enhancedQuickBrownFoxRepeated := strings.Join([]string{enhancedQuickBrownFox, enhancedQuickBrownFox, enhancedQuickBrownFox, enhancedQuickBrownFox, enhancedQuickBrownFox}, " ")

	testCases := []string{
		"the quick brown fox jumped over the lazy dog",
		enhancedQuickBrownFox,
		enhancedQuickBrownFoxRepeated,
	}

	for _, str := range testCases {
		wt := NewWaveletTreeFromString(str)
		actual := ""
		for i := 0; i < len(str); i++ {
			actual += string(wt.Access(i))
		}
		if actual != str {
			t.Fatalf("expected to rebuild:\n%s\nbut instead got:\n%s", str, actual)
		}
	}
}
