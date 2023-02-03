package align_test

import (
	"fmt"

	"github.com/TimothyStiles/poly/align"
)

func ExampleNeedlemanWunsch() {
	a := "GATTACA"
	b := "GCATGCU"
	score, alignA, alignB := align.NeedlemanWunsch(a, b)
	fmt.Printf("score: %d, A: %s, B: %s", score, alignA, alignB)
}
