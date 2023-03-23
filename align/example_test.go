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

	// Output: score: 0, A: G-ATTACA, B: GCA-TGCU
}

func ExampleSmithWaterman() {
	a := "GATTACA"
	b := "GCATGCU"

	scoring := align.NewScoring()
	score, alignA, alignB := align.SmithWaterman(a, b, scoring)

	fmt.Printf("score: %d, A: %s, B: %s", score, alignA, alignB)

	// Output: score: 2, A: AT, B: AT
}
