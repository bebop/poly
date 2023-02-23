package align_test

import (
	"testing"

	"github.com/TimothyStiles/poly/align"
)

func TestNeedlemanWunsch(t *testing.T) {

	scoring := align.NewScoring()

	a := "GATTACA"
	b := "GCATGCU"
	score, alignA, alignB := align.NeedlemanWunsch(a, b, scoring)
	if score != 0 {
		t.Errorf("score: %d, A: %s, B: %s", score, alignA, alignB)
	}

	c := "GATTACA"
	d := "GATTACA"

	score, alignC, alignD := align.NeedlemanWunsch(c, d, scoring)
	if score != 7 {
		t.Errorf("score: %d, A: %s, B: %s", score, alignC, alignD)
	}

	e := "GATTACA"
	f := "GAT"

	score, alignE, alignF := align.NeedlemanWunsch(e, f, scoring)
	if score != -1 {
		t.Errorf("score: %d, A: %s, B: %s", score, alignE, alignF)
	}

	// check zero length string
	g := ""
	h := "GAT"

	score, alignG, alignH := align.NeedlemanWunsch(g, h, scoring)
	if score != -3 {
		t.Errorf("score: %d, A: %s, B: %s", score, alignG, alignH)
	}

	// check zero length strings
	i := ""
	j := ""

	score, alignI, alignJ := align.NeedlemanWunsch(i, j, scoring)
	if score != 0 {
		t.Errorf("score: %d, A: %s, B: %s", score, alignI, alignJ)
	}

	// check 1 length strings
	k := "G"
	l := "A"

	score, alignK, alignL := align.NeedlemanWunsch(k, l, scoring)
	if score != -1 {
		t.Errorf("score: %d, A: %s, B: %s", score, alignK, alignL)
	}

	// check 1 length strings of same letter
	m := "G"
	n := "G"

	score, alignM, alignN := align.NeedlemanWunsch(m, n, scoring)
	if score != 1 {
		t.Errorf("score: %d, A: %s, B: %s", score, alignM, alignN)
	}

	// check 1 length strings with n length string
	o := "G"
	p := "GATTACA"

	score, alignO, alignP := align.NeedlemanWunsch(o, p, scoring)
	if score != -5 {
		t.Errorf("score: %d, A: %s, B: %s", score, alignO, alignP)
	}
}


func TestSmithWaterman(t *testing.T) {

	scoring := align.Scoring{
		Match:      3,
		Mismatch:   -3,
		GapPenalty: -2,
	}

	// Wikipedia example: https://en.wikipedia.org/wiki/Smith%E2%80%93Waterman_algorithm#Example
	a := "TGTTACGG"
	b := "GGTTGACTA"

	score, alignA, alignB := align.SmithWaterman(a, b, scoring)
	if score != 13 {
		t.Errorf("score: %d, A: %s, B: %s", score, alignA, alignB)
	}
	if alignA != "GTT-AC" {
		t.Errorf("Alignment is %s, expected GTT-AC", alignB)
	}
	if alignB != "GTTGAC" {
		t.Errorf("Alignment is %s, expected GTTGAC", alignA)
	}

	c := "GATTACA"
	d := "GCATGCU"

	score, alignC, alignD := align.SmithWaterman(c, d, scoring)

	if score != 7 {
		t.Errorf("score: %d, A: %s, B: %s", score, alignC, alignD)
	}
	if alignC != "G-AT" {
		t.Errorf("Alignment is %s, expected G-AT", alignC)
	}
	if alignD != "GCAT" {
		t.Errorf("Alignment is %s, expected GCAT", alignD)
	}

	// Test edge cases

	// check zero length string
	e := ""
	f := "GAT"

	score, alignE, alignF := align.SmithWaterman(e, f, scoring)
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

	score, alignG, alignH := align.SmithWaterman(g, h, scoring)
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

	score, alignI, alignJ := align.SmithWaterman(i, j, scoring)
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

	score, alignK, alignL := align.SmithWaterman(k, l, scoring)
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

	score, alignM, alignN := align.SmithWaterman(m, n, scoring)
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
