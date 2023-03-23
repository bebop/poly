/*
Package matrix provides a struct for substitution matrices and a struct for scoring matrices.
*/

package matrix

import (
	"fmt"

	"github.com/TimothyStiles/poly/alphabet"
)

// SubstitutionMatrix is a struct that holds a substitution matrix and the two alphabets that the matrix is defined over.
type SubstitutionMatrix struct {
	FirstAlphabet  *alphabet.Alphabet
	SecondAlphabet *alphabet.Alphabet
	scores         [][]int
}

// NewSubstitutionMatrix creates a new substitution matrix from two alphabets and a 2D array of scores.
func NewSubstitutionMatrix(firstAlphabet, secondAlphabet *alphabet.Alphabet, scores [][]int) (*SubstitutionMatrix, error) {
	if len(firstAlphabet.Symbols()) != len(scores) || len(secondAlphabet.Symbols()) != len(scores[0]) {
		return nil, fmt.Errorf("invalid dimensions of substitution matrix")
	}
	return &SubstitutionMatrix{firstAlphabet, secondAlphabet, scores}, nil
}

// Score returns the score of two symbols in the substitution matrix.
func (matrix *SubstitutionMatrix) Score(a, b string) (int, error) {
	firstSymbolIndex, err := matrix.FirstAlphabet.Encode(a)
	if err != nil {
		return 0, err
	}
	secondSymbolIndex, err := matrix.SecondAlphabet.Encode(b)
	if err != nil {
		return 0, err
	}
	return matrix.scores[firstSymbolIndex][secondSymbolIndex], nil
}
