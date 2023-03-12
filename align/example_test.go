package align_test

import (
	"fmt"

	"github.com/TimothyStiles/poly/align"
	"github.com/TimothyStiles/poly/align/matrix"
	"github.com/TimothyStiles/poly/alphabet"
)

func ExampleNeedlemanWunsch() {
	a := "GATTACA"
	b := "GCATGCU"

	alphabet := alphabet.NewAlphabet([]string{"-", "A", "C", "G", "T"})
	subMatrix, err := matrix.NewSubstitutionMatrix(alphabet, alphabet, matrix.NUC_4)
	if err != nil {
		fmt.Println(err)
	}

	scoring := align.NewScoring(subMatrix, -1)
	score, alignA, alignB := align.NeedlemanWunsch(a, b, scoring)

	fmt.Printf("score: %d, A: %s, B: %s", score, alignA, alignB)

	// Output: score: 0, A: G-ATTACA, B: GCA-TGCU
}

func ExampleSmithWaterman() {
	a := "GATTACA"
	b := "GCATGCU"

	alphabet := alphabet.NewAlphabet([]string{"-", "A", "C", "G", "T"})
	subMatrix, err := matrix.NewSubstitutionMatrix(alphabet, alphabet, matrix.NUC_4)
	if err != nil {
		fmt.Println(err)
	}
	scoring := align.NewScoring(subMatrix, -1)
	score, alignA, alignB := align.SmithWaterman(a, b, scoring)

	fmt.Printf("score: %d, A: %s, B: %s", score, alignA, alignB)

	// Output: score: 2, A: AT, B: AT
}
