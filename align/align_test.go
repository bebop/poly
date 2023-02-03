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
}
