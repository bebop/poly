package matrix_test

import (
	"testing"

	"github.com/bebop/poly/alphabet"
	"github.com/bebop/poly/search/align/matrix"
	"github.com/stretchr/testify/assert"
)

func TestSubstitutionMatrix(t *testing.T) {
	alpha1 := alphabet.NewAlphabet([]string{"-", "A", "C", "G", "T"})
	alpha2 := alphabet.NewAlphabet([]string{"-", "A", "C", "G", "T"})
	NUC_4 := [][]int{ //nolint:stylecheck
		/*       - A C G T */
		/* - */ {0, 0, 0, 0, 0},
		/* A */ {0, 5, -4, -4, -4},
		/* C */ {0, -4, 5, -4, -4},
		/* G */ {0, -4, -4, 5, -4},
		/* T */ {0, -4, -4, -4, 5},
	}
	subMat, err := matrix.NewSubstitutionMatrix(alpha1, alpha2, NUC_4)

	assert.Nil(t, err)

	testCases := []struct {
		symbol1 string
		symbol2 string
		score   int
	}{
		{"A", "A", 5},
		{"A", "C", -4},
		{"C", "T", -4},
		{"-", "-", 0},
	}

	for _, tc := range testCases {
		sym1, _ := alpha1.Encode(tc.symbol1)
		sym2, _ := alpha2.Encode(tc.symbol2)
		score, err := subMat.Score(tc.symbol1, tc.symbol2)

		assert.Nil(t, err)
		assert.Equal(t, NUC_4[sym1][sym2], score)

		if score != tc.score {
			t.Errorf("Expected score %d for symbols %s and %s, but got %d", tc.score, tc.symbol1, tc.symbol2, score)
		}
	}
}
