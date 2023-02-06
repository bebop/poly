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
	scoring := align.NewScoring()

	a := "GATTACA"
	b := "AGCAT"

	score, alignA, alignB := align.SmithWaterman(a, b, scoring)
	if score != 3 {
		t.Errorf("score: %d, A: %s, B: %s", score, alignA, alignB)
	}

	c := "GATTACA"
	d := ""

	score, alignC, alignD := align.SmithWaterman(c, d, scoring)
	if score != 0 {
		t.Errorf("score: %d. A: %s, B: %s", score, alignC, alignD)

	}

	e := "GATTACA"
	f := "GATTACA"

	score, alignE, alignF := align.SmithWaterman(e, f, scoring)
	if score != 7 {
		t.Errorf("score: %d. A: %s, B: %s", score, alignE, alignF)
	}

	g := "GATTACA"
	h := "G"

	score, alignG, alignH := align.SmithWaterman(g, h, scoring)
	if score != 1 {
		t.Errorf("score: %d. A: %s, B: %s", score, alignG, alignH)
	}
}
