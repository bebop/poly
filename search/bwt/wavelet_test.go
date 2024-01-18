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
	wt, err := newWaveletTreeFromString(testStr)
	if err != nil {
		t.Fatal(err)
	}

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
	wt, err := newWaveletTreeFromString(testStr)
	if err != nil {
		t.Fatal(err)
	}

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
	wt, err := newWaveletTreeFromString(testStr)
	if err != nil {
		t.Fatal(err)
	}

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

// TestWaveletTree_Access_Reconstruction these tests are to ensure that the wavelet tree is formed correctly. If we can reconstruct the string, we can be
// fairly confident that the WaveletTree is well formed.
func TestWaveletTree_Access_Reconstruction(t *testing.T) {
	// Build with a fair sized alphabet
	enhancedQuickBrownFox := "the quick brown fox jumps over the lazy dog with an overt frown after fumbling its parallelogram shaped bananagram all around downtown"
	enhancedQuickBrownFoxRepeated := strings.Join([]string{enhancedQuickBrownFox, enhancedQuickBrownFox, enhancedQuickBrownFox, enhancedQuickBrownFox, enhancedQuickBrownFox}, " ")
	// Make it very large to account for any succinct data structures being used under the hood. For example, this helped uncover and errors
	// diagnose issues with the Jacobson's Rank used under the hood.
	enhancedQuickBrownFoxSuperLarge := ""
	for i := 0; i < 100; i++ {
		enhancedQuickBrownFoxSuperLarge += enhancedQuickBrownFoxRepeated
	}

	testCases := []string{
		"the quick brown fox jumped over the lazy dog",
		"the quick brown fox jumped over the lazy dog!", // odd numbered alphabet
		enhancedQuickBrownFox,
		enhancedQuickBrownFoxRepeated,
		enhancedQuickBrownFoxSuperLarge,
	}

	for _, str := range testCases {
		wt, err := newWaveletTreeFromString(str)
		if err != nil {
			t.Fatal(err)
		}
		actual := wt.reconstruct()
		if actual != str {
			t.Fatalf("expected to rebuild:\n%s\nbut instead got:\n%s", str, actual)
		}
	}
}

func TestWaveletTreeEmptyStr(t *testing.T) {
	str := ""
	_, err := newWaveletTreeFromString(str)
	if err == nil {
		t.Fatal("expected error but got nil")
	}
}

func TestWaveletTreeSingleChar(t *testing.T) {
	char := "l"
	wt, err := newWaveletTreeFromString(char)
	if err != nil {
		t.Fatal(err)
	}
	r := wt.Rank(char[0], 1)
	s := wt.Select(char[0], 0)
	a := wt.Access(0)

	if r != 1 {
		t.Fatalf("expected Rank(%s, %d) to be %d but got %d", char, 1, 1, r)
	}
	if s != 0 {
		t.Fatalf("expected Select(%s, %d) to be %d but got %d", char, 0, 0, s)
	}
	if a != char[0] {
		t.Fatalf("expected Access(%d) to be %d but got %d", 1, 1, s)
	}
}

func TestWaveletTreeSingleAlpha(t *testing.T) {
	str := "lll"
	wt, err := newWaveletTreeFromString(str)
	if err != nil {
		t.Fatal(err)
	}
	r := wt.Rank(str[0], 1)
	s := wt.Select(str[0], 1)
	a := wt.Access(0)

	if r != 1 {
		t.Fatalf("expected Rank(%s, %d) to be %d but got %d", str, 1, 1, r)
	}
	if s != 1 {
		t.Fatalf("expected Select(%s, %d) to be %d but got %d", str, 1, 1, s)
	}
	if a != str[0] {
		t.Fatalf("expected Access(%d) to be %d but got %d", 1, 1, s)
	}
}
func TestBuildWaveletTree_ZeroAlpha(t *testing.T) {
	bytes := []byte("AAAACCCCTTTTGGGG")
	alpha := []charInfo{}

	root := buildWaveletTree(0, alpha, bytes)

	if root != nil {
		t.Fatalf("expected root to be nil but got %v", root)
	}
}
func TestWaveletTree_LookupCharInfo_Panic(t *testing.T) {
	wt := waveletTree{
		alpha: []charInfo{},
	}

	defer func() {
		if r := recover(); r == nil {
			t.Errorf("expected panic but got nil")
		}
	}()

	wt.lookupCharInfo('B')
}
