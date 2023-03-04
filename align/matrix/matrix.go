package matrix

import (
	"fmt"

	"github.com/TimothyStiles/poly/alphabet"
)

type SubstitutionMatrix struct {
	alpha1 *alphabet.Alphabet
	alpha2 *alphabet.Alphabet
	scores [][]int
}

func NewSubstitutionMatrix(alpha1, alpha2 *alphabet.Alphabet, scores [][]int) (*SubstitutionMatrix, error) {
	if len(alpha1.Symbols()) != len(scores) || len(alpha2.Symbols()) != len(scores[0]) {
		return nil, fmt.Errorf("invalid dimensions of substitution matrix")
	}
	return &SubstitutionMatrix{alpha1, alpha2, scores}, nil
}

func (s *SubstitutionMatrix) Score(a, b string) (int, error) {
	i, err := s.alpha1.Encode(a)
	if err != nil {
		return 0, err
	}
	j, err := s.alpha2.Encode(b)
	if err != nil {
		return 0, err
	}
	return s.scores[i][j], nil
}
