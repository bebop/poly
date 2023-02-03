package align_test

import (
	"fmt"

	"github.com/TimothyStiles/poly/align"
)

func ExampleNeedlemanWunsch() {
	a := "GATTACA"
	b := "GCATGCU"
	scoring := align.NewScoring()
	score, alignA, alignB := align.NeedlemanWunsch(a, b, scoring)
	fmt.Printf("score: %d, A: %s, B: %s", score, alignA, alignB)
}
