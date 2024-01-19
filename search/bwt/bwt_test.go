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

	bwt, err := New(testStr)
	if err != nil {
		t.Fatal(err)
	}

	testTable := []BWTCountTestCase{
		{"uick", 3},
		{"over", 6},
		{"own", 12},
		{"ana", 6},
		{"an", 9},
		{"na", 9},
		{"rown", 6},
		{"frown", 3},
		{"brown", 3},
		{"all", 6},
		{"alle", 3},
		{"alla", 3},
		{"l", 21},
		{"the", 6},
		{"town", 3},
		{"townthe", 2},
		{"nt", 5},
		// patterns that should not exist
		{"@", 0},
		{"zzz", 0},
		{"clown", 0},
		{"crown", 0},
		{"spark", 0},
		{"brawn", 0},
		{"overtly", 0},
	}

	for _, v := range testTable {
		count, err := bwt.Count(v.seq)
		if err != nil {
			t.Fatalf("seq=%s unexpectedError=%s", v.seq, err)
		}
		if count != v.expected {
			t.Fatalf("seq=%s expectedCount=%v actualCount=%v", v.seq, v.expected, count)
		}
	}
}

func TestBWT_Count_EmptyPattern(t *testing.T) {
	testStr := "banana"
	bwt, err := New(testStr)
	if err != nil {
		t.Fatal(err)
	}
	_, err = bwt.Count("")
	if err == nil {
		t.Fatal("Expected error for empty pattern but got nil")
	}
}

type BWTLocateTestCase struct {
	seq      string
	expected []int
}

func TestBWT_Locate(t *testing.T) {
	baseTestStr := "thequickbrownfoxjumpsoverthelazydogwithanovertfrownafterfumblingitsparallelogramshapedbananagramallarounddowntown" // len == 112
	testStr := strings.Join([]string{baseTestStr, baseTestStr, baseTestStr}, "")

	bwt, err := New(testStr)
	if err != nil {
		t.Fatal(err)
	}

	testTable := []BWTLocateTestCase{
		{"uick", []int{4, 117, 230}},
		{"over", []int{21, 41, 134, 154, 247, 267}},
		{"own", []int{10, 48, 106, 110, 123, 161, 219, 223, 236, 274, 332, 336}},
		{"ana", []int{87, 89, 200, 202, 313, 315}},
		{"an", []int{39, 87, 89, 152, 200, 202, 265, 313, 315}},
		{"na", []int{50, 88, 90, 163, 201, 203, 276, 314, 316}},
		{"rown", []int{9, 47, 122, 160, 235, 273}},
		{"frown", []int{46, 159, 272}},
		{"brown", []int{8, 121, 234}},
		{"all", []int{70, 96, 183, 209, 296, 322}},
		{"alle", []int{70, 183, 296}},
		{"alla", []int{96, 209, 322}},
		{"l", []int{28, 60, 71, 72, 74, 97, 98, 141, 173, 184, 185, 187, 210, 211, 254, 286, 297, 298, 300, 323, 324}},
		{"the", []int{0, 25, 113, 138, 226, 251}},
		{"town", []int{109, 222, 335}},
		{"townthe", []int{109, 222}},
		{"nt", []int{108, 112, 221, 225, 334}},
		{"overtly", nil},

		// patterns that should not exist
		{"zzz", nil},
		{"@", nil},
		{"clown", nil},
		{"crown", nil},
		{"spark", nil},
		{"brawn", nil},
	}

	for _, v := range testTable {
		offsets, err := bwt.Locate(v.seq)
		if err != nil {
			t.Fatalf("seq=%s unexpectedError=%s", v.seq, err)
		}
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

func TestBWT_Locate_EmptyPattern(t *testing.T) {
	testStr := "banana"
	bwt, err := New(testStr)
	if err != nil {
		t.Fatal(err)
	}
	_, err = bwt.Locate("")
	if err == nil {
		t.Fatal("Expected error for empty pattern but got nil")
	}
}

type BWTExtractTestCase struct {
	start    int
	end      int
	expected string
}

func TestBWT_Extract(t *testing.T) {
	baseTestStr := "thequickbrownfoxjumpsoverthelazydogwithanovertfrownafterfumblingitsparallelogramshapedbananagramallarounddowntown" // len == 112
	testStr := strings.Join([]string{baseTestStr, baseTestStr, baseTestStr}, "")

	bwt, err := New(testStr)
	if err != nil {
		t.Fatal(err)
	}

	testTable := []BWTExtractTestCase{
		{4, 8, "uick"},
		{117, 121, "uick"},
		{230, 234, "uick"},
		{0, 3, "the"},
		{25, 28, "the"},
		{113, 116, "the"},
		{138, 141, "the"},
		{226, 229, "the"},
		{251, 254, "the"},
		{21, 25, "over"},
		{41, 45, "over"},
		{134, 138, "over"},
		{154, 158, "over"},
		{247, 251, "over"},
		{267, 271, "over"},
		{10, 13, "own"},
		{48, 51, "own"},
		{106, 109, "own"},
		{123, 126, "own"},
		{161, 164, "own"},
		{219, 222, "own"},
		{223, 226, "own"},
		{236, 239, "own"},
		{274, 277, "own"},
		{332, 335, "own"},
		{336, 339, "own"},
		{87, 90, "ana"},
		{89, 92, "ana"},
		{200, 203, "ana"},
		{202, 205, "ana"},
		{313, 316, "ana"},
		{315, 318, "ana"},
		{39, 41, "an"},
		{87, 89, "an"},
		{152, 154, "an"},
		{200, 202, "an"},
		{202, 204, "an"},
		{265, 267, "an"},
		{313, 315, "an"},
		{50, 52, "na"},
		{88, 90, "na"},
		{163, 165, "na"},
		{201, 203, "na"},
		{203, 205, "na"},
		{276, 278, "na"},
		{314, 316, "na"},
		{316, 318, "na"},
		{9, 13, "rown"},
		{47, 51, "rown"},
		{122, 126, "rown"},
		{160, 164, "rown"},
		{235, 239, "rown"},
		{273, 277, "rown"},
		{109, 116, "townthe"},
		{222, 229, "townthe"},
	}

	for _, v := range testTable {
		str, err := bwt.Extract(v.start, v.end)
		if err != nil {
			t.Fatalf("extractRange=(%d, %d) unexpectedError=%s", v.start, v.end, err)
		}
		if str != v.expected {
			t.Fatalf("extractRange=(%d, %d) expected=%s actual=%s", v.start, v.end, v.expected, str)
		}
	}
}

func TestBWT_Extract_InvalidRanges(t *testing.T) {
	testStr := "banana"
	bwt, err := New(testStr)
	if err != nil {
		t.Fatal(err)
	}
	_, err = bwt.Extract(5, 4)
	if err == nil {
		t.Fatal("Expected error for invalid range but got nil")
	}
	_, err = bwt.Extract(4, 4)
	if err == nil {
		t.Fatal("Expected error for invalid range but got nil")
	}
}

func TestBWT_Extract_DoNotAllowExtractionOfLastNullChar(t *testing.T) {
	testStr := "banana"

	bwt, err := New(testStr)
	if err != nil {
		t.Fatal(err)
	}

	str, err := bwt.Extract(0, 6)
	if err != nil {
		t.Fatalf("extractRange=(%d, %d) unexpectedError=%s", 0, 6, err)
	}
	if str != testStr {
		t.Fatalf("extractRange=(%d, %d) expected=%s actual=%s", 0, 6, testStr, str)
	}

	_, err = bwt.Extract(0, 7)

	if err == nil {
		t.Fatalf("extractRange=(%d, %d) expected err but was nil", 0, 7)
	}

	if !strings.Contains(err.Error(), "exceeds the max range") {
		t.Fatalf("expected error to contain \"exceeds the max range\" but received \"%s\"", err)
	}
}

func TestBWT_GetTransform(t *testing.T) {
	baseTestStr := "thequickbrownfoxjumpsoverthelazydogwithanovertfrownafterfumblingitsparallelogramshapedbananagramallarounddowntown" // len == 112
	testStr := strings.Join([]string{baseTestStr, baseTestStr, baseTestStr}, "")

	bwt, err := New(testStr)
	if err != nil {
		t.Fatal(err)
	}

	expected := "nnnnnnnmmmrrrrrrrrrnnnbbbhhhhhhppplllllldddmmmkkkiiieeennnyyydddppphhhlllhhhtttvvvvvvnnntttaaarrrnnnaaaooooootttsssttttttuuulllwwwgggxxxcccllleeelllbbbaaaaaaeeeaaauuuuuuaaawwwwaaaaaauuuwwwiiiaaawwwwwllldddrrrnnnssstrrrrrrttdddfffsssaaammmeeeaaaggggggeeeaaafffbbbeeeeeemmmppptttfffrrriiirrrnn$nnniiiqqqfffjjjooooooooogggooooooooooooooozzzaaa"
	actual := bwt.GetTransform()
	if expected != actual {
		t.Fatalf("expected did not match actual\nexpected:\t%s\nactual:\t%s", expected, actual)
	}
}

func TestBWT_Len(t *testing.T) {
	testStr := "banana"

	bwt, err := New(testStr)
	if err != nil {
		t.Fatal(err)
	}

	if bwt.Len() != len(testStr) {
		t.Fatalf("expected Len to be %d but got %d", len(testStr), bwt.Len())
	}
}

type sparseOnesTestCase struct {
	pos      int
	expected int
}

func TestSparseOnesRS_RankOnes(t *testing.T) {
	runs := runInfo{
		6,
		12,
		33,
		99,
		204,
		205,
		300,
		302,
		305,
		306,
		999,
	}

	testCases := []sparseOnesTestCase{
		{0, 0},
		{4, 0},
		{5, 0},
		{6, 0},

		{7, 0},
		{8, 0},
		{9, 0},
		{11, 0},
		{12, 1},

		{13, 1},
		{15, 1},
		{22, 1},
		{32, 1},
		{33, 2},

		{56, 2},
		{64, 2},
		{65, 2},
		{79, 2},
		{98, 2},
		{99, 3},

		{100, 3},
		{112, 3},
		{168, 3},
		{197, 3},
		{199, 3},
		{203, 3},
		{204, 4},

		{205, 5},

		{206, 5},
		{271, 5},
		{299, 5},
		{300, 6},

		{301, 6},
		{302, 7},

		{302, 7},
		{303, 7},
		{304, 7},
		{305, 8},

		{306, 9},

		{999, 10},
	}

	for _, v := range testCases {
		actual := runs.Rank(v.pos)
		if actual != v.expected {
			t.Fatalf("expected RankOnes(%d) to be %d but got %d", v.pos, v.expected, actual)
		}
	}
}

func TestNewBWTWithSequenceContainingNullChar(t *testing.T) {
	nc := nullChar
	testStr := "banana" + nc

	_, err := New(testStr)
	if err == nil {
		t.Fatal("expected error but got nil")
	}
}

func TestNewBWTEmptySequence(t *testing.T) {
	testStr := ""

	_, err := New(testStr)
	if err == nil {
		t.Fatal("expected error but got nil")
	}
}

// TestBWTReconstruction this helps us ensure that the LF mapping is correct and that the suffix array lookup
// must be well formed. Otherwise, we would not be able to recreate the original sequence.
func TestBWTReconstruction(t *testing.T) {
	baseTestStr := "thequickbrownfoxjumpsoverthelazydogwithanovertfrownafterfumblingitsparallelogramshapedbananagramallarounddowntown"
	testStr := strings.Join([]string{baseTestStr, baseTestStr, baseTestStr}, "")

	bwt, err := New(testStr)
	if err != nil {
		t.Fatal(err)
	}

	extracted, err := bwt.Extract(0, bwt.Len())
	if err != nil {
		t.Fatal(err)
	}
	if extracted != testStr {
		t.Log("Reconstruction failed")
		t.Log("Expected:\t", testStr)
		t.Log("Actual:\t", extracted)
		t.Fail()
	}

	// This will either result in an even or all alphabet. The alphabet matters.
	testStrWithOneMoreAlpha := testStr + "!"
	bwt, err = New(testStrWithOneMoreAlpha)
	if err != nil {
		t.Fatal(err)
	}
	extracted, err = bwt.Extract(0, bwt.Len())
	if err != nil {
		t.Fatal(err)
	}
	if extracted != testStrWithOneMoreAlpha {
		t.Log("Reconstruction failed with extra alpha character")
		t.Log("Expected:\t", testStr)
		t.Log("Actual:\t", extracted)
		t.Fail()
	}
}

func TestBWTStartError(t *testing.T) {
	testStr := "banana"

	bwt, err := New(testStr)
	if err != nil {
		t.Fatal(err)
	}

	_, err = bwt.Extract(-1, 6)
	if err == nil {
		t.Fatal("expected error but got nil")
	}
}
func TestBWT_GetFCharPosFromOriginalSequenceCharPos_Panic(t *testing.T) {
	testStr := "banana"
	bwt, err := New(testStr)
	if err != nil {
		t.Fatal(err)
	}

	// Call the function with an invalid original position
	originalPos := -1
	defer func() {
		if r := recover(); r == nil {
			t.Errorf("Expected panic, but it did not occur")
		}
	}()
	bwt.getFCharPosFromOriginalSequenceCharPos(originalPos)
}
func TestBWT_LFSearch_InvalidChar(t *testing.T) {
	testStr := "banana"
	bwt, err := New(testStr)
	if err != nil {
		t.Fatal(err)
	}

	pattern := "x" // Invalid character

	result := bwt.lfSearch(pattern)

	if result.start != 0 || result.end != 0 {
		t.Fatalf("Expected search range to be (0, 0), but got (%d, %d)", result.start, result.end)
	}
}
func TestBWT_LookupSkipByOffset_PanicOffsetExceedsMaxBound(t *testing.T) {
	testStr := "banana"
	bwt, err := New(testStr)
	if err != nil {
		t.Fatal(err)
	}

	offset := bwt.getLenOfOriginalStringWithNullChar()
	defer func() {
		if r := recover(); r == nil {
			t.Errorf("Expected panic, but it did not occur")
		}
	}()
	bwt.lookupSkipByOffset(offset)
}

func TestBWT_LookupSkipByOffset_PanicOffsetExceedsMinBound(t *testing.T) {
	testStr := "banana"
	bwt, err := New(testStr)
	if err != nil {
		t.Fatal(err)
	}

	offset := -1
	defer func() {
		if r := recover(); r == nil {
			t.Errorf("Expected panic, but it did not occur")
		}
	}()
	bwt.lookupSkipByOffset(offset)
}

func TestBWTRecovery(t *testing.T) {
	// Test panic recovery for bwtRecovery function
	var err error
	operation := "test operation"

	defer func() {
		if err == nil {
			t.Fatal("expected bwtRecovery to recover from the panic and set an error message, but got nil")
		}
	}()
	defer bwtRecovery(operation, &err)
	doPanic()
}

func doPanic() {
	panic("test panic")
}
