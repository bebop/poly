package bwt

import (
	"fmt"
	"log"
	"strings"
	"testing"

	"golang.org/x/exp/slices"
)

func ExampleBWT_Count() {
	inputSequence := "AACCTGCCGTCGGGGCTGCCCGTCGCGGGACGTCGAAACGTGGGGCGAAACGTG"

	bwt, err := New(inputSequence)
	if err != nil {
		log.Fatal(err)
	}

	fmt.Println(bwt.Count("CG"))
	// Output: 10
}

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
		count, err := bwt.Count(v.seq)
		if err != nil {
			t.Fatalf("seq=%s unexpectedError=%s", v.seq, err)
		}
		if count != v.expected {
			t.Fatalf("seq=%s expectedCount=%v actualCount=%v", v.seq, v.expected, count)
		}
	}
}

func ExampleBWT_Locate() {
	inputSequence := "AACCTGCCGTCGGGGCTGCCCGTCGCGGGACGTCGAAACGTGGGGCGAAACGTG"

	bwt, err := New(inputSequence)
	if err != nil {
		log.Fatal(err)
	}

	offsets, err := bwt.Locate("CG")
	if err != nil {
		log.Fatal(err)
	}
	slices.Sort(offsets)
	fmt.Println(offsets)
	// Output: [7 10 20 23 25 30 33 38 45 50]
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

func ExampleBWT_Extract() {
	inputSequence := "AACCTGCCGTCGGGGCTGCCCGTCGCGGGACGTCGAAACGTGGGGCGAAACGTG"

	bwt, err := New(inputSequence)
	if err != nil {
		log.Fatal(err)
	}

	fmt.Println(bwt.Extract(48, 54))
	// Output: AACGTG
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
