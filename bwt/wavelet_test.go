package bwt

import "testing"

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

func TestWaveletTree_Rank(t *testing.T) {
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
