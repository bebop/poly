package align_test

import (
	"testing"

	"github.com/TimothyStiles/poly/align"
)

func TestNeedlemanWunsch(t *testing.T) {
	a := "GATTACA"
	b := "GCATGCU"
	score, alignA, alignB := align.NeedlemanWunsch(a, b)
	if score != 0 {
		t.Errorf("score: %d, A: %s, B: %s", score, alignA, alignB)
	}

	c := "GATTACA"
	d := "GATTACA"

	score, alignC, alignD := align.NeedlemanWunsch(c, d)
	if score != 7 {
		t.Errorf("score: %d, A: %s, B: %s", score, alignC, alignD)
	}
}
