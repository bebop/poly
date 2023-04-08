package fold_test

import (
	"math"
	"strings"
	"testing"

	"github.com/TimothyStiles/poly/fold"
	"github.com/stretchr/testify/assert"
	"github.com/stretchr/testify/require"
)

func TestFold(t *testing.T) {
	t.Run("FoldCache", func(t *testing.T) {
		seq := "ATGGATTTAGATAGAT"
		fc, err := fold.NewFoldingContext(seq, 37.0)
		require.NoError(t, err)
		seq_dg, err := fold.MinimumFreeEnergy(seq, 37.0)
		require.NoError(t, err)

		assert.InDelta(t, seq_dg, fc.W[0][len(seq)-1].E, 1)
	})
	t.Run("FoldDNA", func(t *testing.T) {
		// unafold's estimates for free energy estimates of DNA oligos
		unafold_dgs := map[string]float64{
			"GGGAGGTCGTTACATCTGGGTAACACCGGTACTGATCCGGTGACCTCCC":                         -10.94, // three branched structure
			"GGGAGGTCGCTCCAGCTGGGAGGAGCGTTGGGGGTATATACCCCCAACACCGGTACTGATCCGGTGACCTCCC": -23.4,  // four branched structure
			"CGCAGGGAUACCCGCG":                         -3.8,
			"TAGCTCAGCTGGGAGAGCGCCTGCTTTGCACGCAGGAGGT": -6.85,
			"GGGGGCATAGCTCAGCTGGGAGAGCGCCTGCTTTGCACGCAGGAGGTCTGCGGTTCGATCCCGCGCGCTCCCACCA": -15.50,
			"TGAGACGGAAGGGGATGATTGTCCCCTTCCGTCTCA":                                         -18.10,
			"ACCCCCTCCTTCCTTGGATCAAGGGGCTCAA":                                              -3.65,
		}

		for seq, ufold := range unafold_dgs {
			d, err := fold.MinimumFreeEnergy(seq, 37.0)
			require.NoError(t, err)
			// accepting a 60% difference
			delta := math.Abs(0.6 * math.Min(d, ufold))
			assert.InDeltaf(t, d, ufold, delta, seq)
		}
	})
	t.Run("FoldRNA", func(t *testing.T) {
		// unafold's estimates for free energy estimates of RNA oligos
		// most tests available at https://github.com/jaswindersingh2/SPOT-RNA/blob/master/sample_inputs/batch_seq.fasta
		unafold_dgs := map[string]float64{
			"ACCCCCUCCUUCCUUGGAUCAAGGGGCUCAA":        -9.5,
			"AAGGGGUUGGUCGCCUCGACUAAGCGGCUUGGAAUUCC": -10.1,
			"UUGGAGUACACAACCUGUACACUCUUUC":           -4.3,
			"AGGGAAAAUCCC":                           -3.3,
			"GCUUACGAGCAAGUUAAGCAAC":                 -4.6,
			"UGGGAGGUCGUCUAACGGUAGGACGGCGGACUCUGGAUCCGCUGGUGGAGGUUCGAGUCCUCCCCUCCCAGCCA":   -32.8,
			"GGGCGAUGAGGCCCGCCCAAACUGCCCUGAAAAGGGCUGAUGGCCUCUACUG":                         -20.7,
			"GGGGGCAUAGCUCAGCUGGGAGAGCGCCUGCUUUGCACGCAGGAGGUCUGCGGUUCGAUCCCGCGCGCUCCCACCA": -31.4,
		}

		for seq, ufold := range unafold_dgs {
			d, err := fold.MinimumFreeEnergy(seq, 37.0)
			require.NoError(t, err)

			// accepting a 50% difference
			delta := math.Abs(0.5 * math.Min(d, ufold))
			assert.InDeltaf(t, d, ufold, delta, seq)
		}
	})
	t.Run("DotBracket", func(t *testing.T) {
		seq := "GGGAGGTCGTTACATCTGGGTAACACCGGTACTGATCCGGTGACCTCCC"
		structs, err := fold.Fold(seq, 37.0)
		require.NoError(t, err)

		assert.Equal(t, "((((((((.((((......))))..((((.......)))).))))))))", fold.DotBracket(structs))
	})
	t.Run("Multibranch", func(t *testing.T) {
		seq := "GGGAGGTCGTTACATCTGGGTAACACCGGTACTGATCCGGTGACCTCCC" // three branch

		structs, err := fold.Fold(seq, 37.0)
		require.NoError(t, err)

		found := false
		foundIJ := fold.Subsequence{7, 41}
		for _, s := range structs {
			if strings.Contains(s.Desc, "BIFURCATION") {
				for _, ij := range s.Inner {
					if ij == foundIJ {
						found = true
					}
				}
			}
		}
		assert.True(t, found, "not found a  BIFURCATION with (7, 41) in ij")
	})
	t.Run("Pair", func(t *testing.T) {
		seq := "ATGGAATAGTG"
		assert.Equal(t, fold.Pair(seq, 0, 1, 9, 10), "AT/TG")
	})
	t.Run("Stack", func(t *testing.T) {
		seq := "GCUCAGCUGGGAGAGC"
		temp := 37.0
		fc, err := fold.NewFoldingContext(seq, temp)
		require.NoError(t, err)

		e := fold.Stack(1, 2, 14, 13, fc)
		assert.InDelta(t, e, -2.1, 0.1)
	})
	t.Run("Bulge", func(t *testing.T) {
		// mock bulge of CAT on one side and AG on other
		// from pg 429 of SantaLucia, 2004
		seq := "ACCCCCATCCTTCCTTGAGTCAAGGGGCTCAA"
		fc, err := fold.NewFoldingContext(seq, 37)
		require.NoError(t, err)

		pair_dg, err := fold.Bulge(5, 7, 18, 17, fc)
		require.NoError(t, err)
		assert.InDelta(t, pair_dg, 3.22, 0.4)
	})
	t.Run("Hairpin", func(t *testing.T) {
		// hairpin = "CCTTGG"
		seq := "ACCCCCTCCTTCCTTGGATCAAGGGGCTCAA"
		i := 11
		j := 16

		fc, err := fold.NewFoldingContext(seq, 37)
		require.NoError(t, err)

		hairpin_dg, err := fold.Hairpin(i, j, fc)
		require.NoError(t, err)
		// this differs from Unafold
		assert.InDelta(t, hairpin_dg, 4.3, 1.0)

		// from page 428 of SantaLucia, 2004
		// hairpin = "CGCAAG"
		seq = "ACCCGCAAGCCCTCCTTCCTTGGATCAAGGGGCTCAA"
		i = 3
		j = 8

		fc, err = fold.NewFoldingContext(seq, 37)
		require.NoError(t, err)

		hairpin_dg, err = fold.Hairpin(i, j, fc)
		require.NoError(t, err)
		assert.InDelta(t, hairpin_dg, 0.67, 0.1)

		seq = "CUUUGCACG"
		i = 0
		j = 8

		fc, err = fold.NewFoldingContext(seq, 37)
		require.NoError(t, err)

		hairpin_dg, err = fold.Hairpin(i, j, fc)
		require.NoError(t, err)
		assert.InDelta(t, hairpin_dg, 4.5, 0.2)
	})
	t.Run("InternalLoop", func(t *testing.T) {
		seq := "ACCCCCTCCTTCCTTGGATCAAGGGGCTCAA"
		i := 6
		j := 21

		fc, err := fold.NewFoldingContext(seq, 37)
		require.NoError(t, err)

		dg, err := fold.InternalLoop(i, i+4, j, j-4, fc)
		require.NoError(t, err)
		assert.InDelta(t, dg, 3.5, 0.1)
	})
	t.Run("W", func(t *testing.T) {
		seq := "GCUCAGCUGGGAGAGC"
		i := 0
		j := 15

		fc, err := fold.NewFoldingContext(seq, 37)
		require.NoError(t, err)

		struc, err := fold.W(i, j, fc)
		require.NoError(t, err)
		assert.InDelta(t, struc.E, -3.8, 0.2)

		seq = "CCUGCUUUGCACGCAGG"
		i = 0
		j = 16

		fc, err = fold.NewFoldingContext(seq, 37)
		require.NoError(t, err)

		struc, err = fold.W(i, j, fc)
		require.NoError(t, err)
		assert.InDelta(t, struc.E, -6.4, 0.2)

		seq = "GCGGUUCGAUCCCGC"
		i = 0
		j = 14

		fc, err = fold.NewFoldingContext(seq, 37)
		require.NoError(t, err)

		struc, err = fold.W(i, j, fc)
		require.NoError(t, err)
		assert.InDelta(t, struc.E, -4.2, 0.2)
	})
}

func NewFoldingContext(seq string, f float64) {
	panic("unimplemented")
}
