package variants

import (
	"fmt"
	"sort"
	"testing"

	"github.com/stretchr/testify/assert"
)

func TestIUPAC(t *testing.T) {
	testSeq := "ATN"
	testVariants := []string{"ATG", "ATA", "ATT", "ATC"}
	testVariantsIUPAC, err := AllVariantsIUPAC(testSeq)
	assert.Nil(t, err)

	sort.Strings(testVariants)
	sort.Strings(testVariantsIUPAC)

	for index := range testVariants {
		assert.Equal(t, testVariants[index], testVariantsIUPAC[index], "IUPAC variant has changed")
	}
}

func TestIUPAC_errors(t *testing.T) {
	testSeq := "ATX"
	seqVariants, err := AllVariantsIUPAC(testSeq)
	assert.NotNil(t, err, "expected error for unsupported IUPAC character, got nil")
	assert.Equal(t, seqVariants, []string{})
}

func ExampleAllVariantsIUPAC() {
	// AllVariantsIUPAC takes a string as input
	// and returns all iupac variants as output
	mendelIUPAC := "ATGGARAAYGAYGARCTN"
	// ambiguous IUPAC codes for most of the sequences that code for the protein MENDEL
	mendelIUPACvariants, _ := AllVariantsIUPAC(mendelIUPAC)
	fmt.Println(mendelIUPACvariants)
	// Output: [ATGGAGAATGATGAGCTG ATGGAGAATGATGAGCTA ATGGAGAATGATGAGCTT ATGGAGAATGATGAGCTC ATGGAGAATGATGAACTG ATGGAGAATGATGAACTA ATGGAGAATGATGAACTT ATGGAGAATGATGAACTC ATGGAGAATGACGAGCTG ATGGAGAATGACGAGCTA ATGGAGAATGACGAGCTT ATGGAGAATGACGAGCTC ATGGAGAATGACGAACTG ATGGAGAATGACGAACTA ATGGAGAATGACGAACTT ATGGAGAATGACGAACTC ATGGAGAACGATGAGCTG ATGGAGAACGATGAGCTA ATGGAGAACGATGAGCTT ATGGAGAACGATGAGCTC ATGGAGAACGATGAACTG ATGGAGAACGATGAACTA ATGGAGAACGATGAACTT ATGGAGAACGATGAACTC ATGGAGAACGACGAGCTG ATGGAGAACGACGAGCTA ATGGAGAACGACGAGCTT ATGGAGAACGACGAGCTC ATGGAGAACGACGAACTG ATGGAGAACGACGAACTA ATGGAGAACGACGAACTT ATGGAGAACGACGAACTC ATGGAAAATGATGAGCTG ATGGAAAATGATGAGCTA ATGGAAAATGATGAGCTT ATGGAAAATGATGAGCTC ATGGAAAATGATGAACTG ATGGAAAATGATGAACTA ATGGAAAATGATGAACTT ATGGAAAATGATGAACTC ATGGAAAATGACGAGCTG ATGGAAAATGACGAGCTA ATGGAAAATGACGAGCTT ATGGAAAATGACGAGCTC ATGGAAAATGACGAACTG ATGGAAAATGACGAACTA ATGGAAAATGACGAACTT ATGGAAAATGACGAACTC ATGGAAAACGATGAGCTG ATGGAAAACGATGAGCTA ATGGAAAACGATGAGCTT ATGGAAAACGATGAGCTC ATGGAAAACGATGAACTG ATGGAAAACGATGAACTA ATGGAAAACGATGAACTT ATGGAAAACGATGAACTC ATGGAAAACGACGAGCTG ATGGAAAACGACGAGCTA ATGGAAAACGACGAGCTT ATGGAAAACGACGAGCTC ATGGAAAACGACGAACTG ATGGAAAACGACGAACTA ATGGAAAACGACGAACTT ATGGAAAACGACGAACTC]
}

func ExampleAllVariantsIUPAC_error() {
	// AllVariantsIUPAC takes a string as input
	// and returns all iupac variants as output
	mendelIUPAC := "ATGGARAAYGAYGARXYZ"
	// ambiguous IUPAC codes for most of the sequences that code for the protein MENDEL
	mendelIUPACvariants, _ := AllVariantsIUPAC(mendelIUPAC)
	fmt.Println(mendelIUPACvariants)
	// Output: []
}

func Test_cartRunes(t *testing.T) {
	for _, tc := range []struct {
		in  [][]rune
		out [][]rune
	}{
		{
			in: [][]rune{
				{1, 2, 3},
			},
			out: [][]rune{
				{1},
				{2},
				{3},
			},
		},
		{
			in:  [][]rune{nil},
			out: nil,
		},
		{
			in: [][]rune{
				{'A', 'T', 'G'},
				{'A', 'T', 'A', 'G'},
			},
			out: [][]rune{
				{'A', 'A'},
				{'A', 'T'},
				{'A', 'A'},
				{'A', 'G'},
				{'T', 'A'},
				{'T', 'T'},
				{'T', 'A'},
				{'T', 'G'},
				{'G', 'A'},
				{'G', 'T'},
				{'G', 'A'},
				{'G', 'G'},
			},
		},
	} {
		runes := cartRune(tc.in...)
		assert.Equal(t, runes, tc.out)
	}
}
