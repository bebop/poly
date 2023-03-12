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

	m := [][]int{
		/*       A C G T U */
		/* A */ {1, -1, -1, -1, -1},
		/* C */ {-1, 1, -1, -1, -1},
		/* G */ {-1, -1, 1, -1, -1},
		/* T */ {-1, -1, -1, 1, -1},
		/* U */ {-1, -1, -1, -1, 1},
	}

	alphabet := alphabet.NewAlphabet([]string{"A", "C", "G", "T", "U"})
	subMatrix, err := matrix.NewSubstitutionMatrix(alphabet, alphabet, m)
	if err != nil {
		fmt.Println(err)
		return
	}

	scoring := align.NewScoring(subMatrix, -1)
	score, alignA, alignB := align.NeedlemanWunsch(a, b, scoring)

	fmt.Printf("score: %d, A: %s, B: %s", score, alignA, alignB)

	// Output: score: 0, A: G-ATTACA, B: GCA-TGCU
}

func ExampleSmithWaterman() {
	a := "GATTACA"
	b := "GCATGCU"

	m := [][]int{
		/*       A C G T U */
		/* A */ {1, -1, -1, -1, -1},
		/* C */ {-1, 1, -1, -1, -1},
		/* G */ {-1, -1, 1, -1, -1},
		/* T */ {-1, -1, -1, 1, -1},
		/* U */ {-1, -1, -1, -1, 1},
	}

	alphabet := alphabet.NewAlphabet([]string{"A", "C", "G", "T", "U"})
	subMatrix, err := matrix.NewSubstitutionMatrix(alphabet, alphabet, m)
	if err != nil {
		fmt.Println(err)
		return
	}
	scoring := align.NewScoring(subMatrix, -1)
	score, alignA, alignB := align.SmithWaterman(a, b, scoring)

	fmt.Printf("score: %d, A: %s, B: %s", score, alignA, alignB)

	// Output: score: 2, A: AT, B: AT
}

func ExampleSmithWaterman_matrix_nuc_4() {
	a := "GATTACA"
	b := "GCATGCT"

	alphabet := alphabet.NewAlphabet([]string{"A", "C", "G", "T", "-"})
	subMatrix, err := matrix.NewSubstitutionMatrix(alphabet, alphabet, matrix.NUC_4)
	if err != nil {
		fmt.Println(err)
		return
	}
	scoring := align.NewScoring(subMatrix, -1)
	score, alignA, alignB := align.SmithWaterman(a, b, scoring)

	fmt.Printf("score: %d, A: %s, B: %s", score, alignA, alignB)

	// Output: score: 15, A: GATTAC, B: GCATGC
}
