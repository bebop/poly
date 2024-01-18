package align_test

import (
	"testing"

	"github.com/bebop/poly/alphabet"
	"github.com/bebop/poly/search/align"
	"github.com/bebop/poly/search/align/matrix"
)

func TestNeedlemanWunsch(t *testing.T) {
	mat := [][]int{
		/*       A C G T U */
		/* A */ {1, -1, -1, -1, -1},
		/* C */ {-1, 1, -1, -1, -1},
		/* G */ {-1, -1, 1, -1, -1},
		/* T */ {-1, -1, -1, 1, -1},
		/* U */ {-1, -1, -1, -1, 1},
	}

	alphabet := alphabet.NewAlphabet([]string{"A", "C", "G", "T", "U"})
	subMatrix, err := matrix.NewSubstitutionMatrix(alphabet, alphabet, mat)
	if err != nil {
		t.Errorf("error: %s", err)
	}
	scoring, err := align.NewScoring(subMatrix, -1)
	if err != nil {
		t.Errorf("error: %s", err)
	}

	a := "GATTACA"
	b := "GCATGCU"
	score, alignA, alignB, err := align.NeedlemanWunsch(a, b, scoring)

	if err != nil {
		t.Errorf("error: %s", err)
	}

	if score != 0 {
		t.Errorf("score: %d, A: %s, B: %s", score, alignA, alignB)
	}

	c := "GATTACA"
	d := "GATTACA"

	score, alignC, alignD, err := align.NeedlemanWunsch(c, d, scoring)

	if err != nil {
		t.Errorf("error: %s", err)
	}

	if score != 7 {
		t.Errorf("score: %d, A: %s, B: %s", score, alignC, alignD)
	}

	e := "GATTACA"
	f := "GAT"

	score, alignE, alignF, err := align.NeedlemanWunsch(e, f, scoring)

	if err != nil {
		t.Errorf("error: %s", err)
	}

	if score != -1 {
		t.Errorf("score: %d, A: %s, B: %s", score, alignE, alignF)
	}

	// check zero length string
	g := ""
	h := "GAT"

	score, alignG, alignH, err := align.NeedlemanWunsch(g, h, scoring)

	if err != nil {
		t.Errorf("error: %s", err)
	}

	if score != -3 {
		t.Errorf("score: %d, A: %s, B: %s", score, alignG, alignH)
	}

	// check zero length strings
	i := ""
	j := ""

	score, alignI, alignJ, err := align.NeedlemanWunsch(i, j, scoring)

	if err != nil {
		t.Errorf("error: %s", err)
	}

	if score != 0 {
		t.Errorf("score: %d, A: %s, B: %s", score, alignI, alignJ)
	}

	// check 1 length strings
	k := "G"
	l := "A"

	score, alignK, alignL, err := align.NeedlemanWunsch(k, l, scoring)

	if err != nil {
		t.Errorf("error: %s", err)
	}

	if score != -1 {
		t.Errorf("score: %d, A: %s, B: %s", score, alignK, alignL)
	}

	// check 1 length strings of same letter
	m := "G"
	n := "G"

	score, alignM, alignN, err := align.NeedlemanWunsch(m, n, scoring)

	if err != nil {
		t.Errorf("error: %s", err)
	}
	if score != 1 {
		t.Errorf("score: %d, A: %s, B: %s", score, alignM, alignN)
	}

	// check 1 length strings with n length string
	o := "G"
	p := "GATTACA"

	score, alignO, alignP, err := align.NeedlemanWunsch(o, p, scoring)

	if err != nil {
		t.Errorf("error: %s", err)
	}

	if score != -5 {
		t.Errorf("score: %d, A: %s, B: %s", score, alignO, alignP)
	}
}

func TestSmithWaterman(t *testing.T) {
	mat := [][]int{
		/*       - A C G T */
		/* - */ {0, 0, 0, 0, 0},
		/* A */ {0, 3, -3, -3, -3},
		/* C */ {0, -3, 3, -3, -3},
		/* G */ {0, -3, -3, 3, -3},
		/* T */ {0, -3, -3, -3, 3},
	}

	alphabet := alphabet.NewAlphabet([]string{"-", "A", "C", "G", "T"})
	subMatrix, _ := matrix.NewSubstitutionMatrix(alphabet, alphabet, mat)
	scoring, err := align.NewScoring(subMatrix, -2)

	if err != nil {
		t.Errorf("error: %s", err)
	}

	// Wikipedia example: https://en.wikipedia.org/wiki/Smith-Waterman_algorithm#Example
	a := "TGTTACGG"
	b := "GGTTGACTA"

	score, alignA, alignB, err := align.SmithWaterman(a, b, scoring)

	if err != nil {
		t.Errorf("error: %s", err)
	}

	if score != 13 {
		t.Errorf("score: %d, A: %s, B: %s", score, alignA, alignB)
	}
	if alignA != "GTT-AC" {
		t.Errorf("Alignment is %s, expected GTT-AC", alignB)
	}
	if alignB != "GTTGAC" {
		t.Errorf("Alignment is %s, expected GTTGAC", alignA)
	}

	c := "ACACACTA"
	d := "AGCACACA"

	score, alignC, alignD, err := align.SmithWaterman(c, d, scoring)

	if err != nil {
		t.Errorf("error: %s", err)
	}

	if score != 17 {
		t.Errorf("score: %d, A: %s, B: %s", score, alignC, alignD)
	}
	if alignC != "A-CACACTA" {
		t.Errorf("Alignment is %s, expected A-CACACTA", alignC)
	}
	if alignD != "AGCACAC-A" {
		t.Errorf("Alignment is %s, expected AGCACAC-A", alignD)
	}

	// Test edge cases

	// check zero length string
	e := ""
	f := "GAT"

	score, alignE, alignF, err := align.SmithWaterman(e, f, scoring)

	if err != nil {
		t.Errorf("error: %s", err)
	}
	if score != 0 {
		t.Errorf("score: %d, A: %s, B: %s", score, alignE, alignF)
	}
	if alignE != "" {
		t.Errorf("Alignment is %s, expected ''", alignE)
	}
	if alignF != "" {
		t.Errorf("Alignment is %s, expected ''", alignF)
	}

	// check zero length strings
	g := ""
	h := ""

	score, alignG, alignH, err := align.SmithWaterman(g, h, scoring)

	if err != nil {
		t.Errorf("error: %s", err)
	}
	if score != 0 {
		t.Errorf("score: %d, A: %s, B: %s", score, alignG, alignH)
	}
	if alignG != "" {
		t.Errorf("Alignment is %s, expected ''", alignG)
	}
	if alignH != "" {
		t.Errorf("Alignment is %s, expected ''", alignH)
	}

	// check 1 length strings
	i := "G"
	j := "A"

	score, alignI, alignJ, err := align.SmithWaterman(i, j, scoring)

	if err != nil {
		t.Errorf("error: %s", err)
	}
	if score != 0 {
		t.Errorf("score: %d, A: %s, B: %s", score, alignI, alignJ)
	}
	if alignI != "" {
		t.Errorf("Alignment is %s, expected ''", alignI)
	}
	if alignJ != "" {
		t.Errorf("Alignment is %s, expected ''", alignJ)
	}

	// check 1 length strings of same letter
	k := "G"
	l := "G"

	score, alignK, alignL, err := align.SmithWaterman(k, l, scoring)

	if err != nil {
		t.Errorf("error: %s", err)
	}
	if score != 3 {
		t.Errorf("score: %d, A: %s, B: %s", score, alignK, alignL)
	}
	if alignK != "G" {
		t.Errorf("Alignment is %s, expected G", alignK)
	}
	if alignL != "G" {
		t.Errorf("Alignment is %s, expected G", alignL)
	}

	// check 1 length strings with n length string
	m := "G"
	n := "GATTACA"

	score, alignM, alignN, err := align.SmithWaterman(m, n, scoring)

	if err != nil {
		t.Errorf("error: %s", err)
	}
	if score != 3 {
		t.Errorf("score: %d, A: %s, B: %s", score, alignM, alignN)
	}
	if alignM != "G" {
		t.Errorf("Alignment is %s, expected G", alignM)
	}
	if alignN != "G" {
		t.Errorf("Alignment is %s, expected G", alignN)
	}
}
