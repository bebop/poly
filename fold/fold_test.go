package fold

import (
	"fmt"
	"math"
	"strings"
	"testing"

	"github.com/stretchr/testify/assert"
	"github.com/stretchr/testify/require"
)

func TestFold(t *testing.T) {
	t.Run("FoldCache", func(t *testing.T) {
		seq := "ATGGATTTAGATAGAT"
		foldContext, err := newFoldingContext(seq, 37.0)
		require.NoError(t, err)
		res, err := Zuker(seq, 37.0)
		require.NoError(t, err)
		seqDg := res.MinimumFreeEnergy()
		require.NoError(t, err)

		assert.InDelta(t, seqDg, foldContext.unpairedMinimumFreeEnergyW[0][len(seq)-1].energy, 1)
	})
	t.Run("FoldDNA", func(t *testing.T) {
		// unafold's estimates for free energy estimates of DNA oligos
		unafoldDgs := map[string]float64{
			"GGGAGGTCGTTACATCTGGGTAACACCGGTACTGATCCGGTGACCTCCC":                         -10.94, // three branched structure
			"GGGAGGTCGCTCCAGCTGGGAGGAGCGTTGGGGGTATATACCCCCAACACCGGTACTGATCCGGTGACCTCCC": -23.4,  // four branched structure
			"CGCAGGGAUACCCGCG":                         -3.8,
			"TAGCTCAGCTGGGAGAGCGCCTGCTTTGCACGCAGGAGGT": -6.85,
			"GGGGGCATAGCTCAGCTGGGAGAGCGCCTGCTTTGCACGCAGGAGGTCTGCGGTTCGATCCCGCGCGCTCCCACCA": -15.50,
			"TGAGACGGAAGGGGATGATTGTCCCCTTCCGTCTCA":                                         -18.10,
			"ACCCCCTCCTTCCTTGGATCAAGGGGCTCAA":                                              -3.65,
		}

		for seq, ufold := range unafoldDgs {
			res, err := Zuker(seq, 37.0)
			require.NoError(t, err)
			d := res.MinimumFreeEnergy()
			// accepting a 60% difference
			delta := math.Abs(0.6 * math.Min(d, ufold))
			assert.InDeltaf(t, d, ufold, delta, seq)
		}
	})
	t.Run("FoldRNA", func(t *testing.T) {
		// unafold's estimates for free energy estimates of RNA oligos
		// most tests available at https://github.com/jaswindersingh2/SPOT-RNA/blob/master/sample_inputs/batch_seq.fasta
		unafoldDgs := map[string]float64{
			"ACCCCCUCCUUCCUUGGAUCAAGGGGCUCAA":        -9.5,
			"AAGGGGUUGGUCGCCUCGACUAAGCGGCUUGGAAUUCC": -10.1,
			"UUGGAGUACACAACCUGUACACUCUUUC":           -4.3,
			"AGGGAAAAUCCC":                           -3.3,
			"GCUUACGAGCAAGUUAAGCAAC":                 -4.6,
			"UGGGAGGUCGUCUAACGGUAGGACGGCGGACUCUGGAUCCGCUGGUGGAGGUUCGAGUCCUCCCCUCCCAGCCA":   -32.8,
			"GGGCGAUGAGGCCCGCCCAAACUGCCCUGAAAAGGGCUGAUGGCCUCUACUG":                         -20.7,
			"GGGGGCAUAGCUCAGCUGGGAGAGCGCCUGCUUUGCACGCAGGAGGUCUGCGGUUCGAUCCCGCGCGCUCCCACCA": -31.4,
		}

		for seq, ufold := range unafoldDgs {
			res, err := Zuker(seq, 37.0)
			require.NoError(t, err)
			d := res.MinimumFreeEnergy()

			// accepting a 50% difference
			delta := math.Abs(0.5 * math.Min(d, ufold))
			assert.InDeltaf(t, d, ufold, delta, seq)
		}
	})
	t.Run("DotBracket", func(t *testing.T) {
		seq := "GGGAGGTCGTTACATCTGGGTAACACCGGTACTGATCCGGTGACCTCCC"
		res, err := Zuker(seq, 37.0)
		require.NoError(t, err)

		assert.Equal(t, "((((((((.((((......))))..((((.......)))).))))))))", res.DotBracket())
	})
	t.Run("multibranch", func(t *testing.T) {
		seq := "GGGAGGTCGTTACATCTGGGTAACACCGGTACTGATCCGGTGACCTCCC" // three branch

		res, err := Zuker(seq, 37.0)
		require.NoError(t, err)

		found := false
		foundIJ := subsequence{7, 41}
		for _, s := range res.structs {
			if strings.Contains(s.description, "BIFURCATION") {
				for _, ij := range s.inner {
					if ij == foundIJ {
						found = true
					}
				}
			}
		}
		assert.True(t, found, "not found a  BIFURCATION with (7, 41) in ij")
	})
	t.Run("pair", func(t *testing.T) {
		seq := "ATGGAATAGTG"
		assert.Equal(t, pair(seq, 0, 1, 9, 10), "AT/TG")
	})
	t.Run("stack", func(t *testing.T) {
		seq := "GCUCAGCUGGGAGAGC"
		temp := 37.0
		foldContext, err := newFoldingContext(seq, temp)
		require.NoError(t, err)

		e := stack(1, 2, 14, 13, foldContext)
		assert.InDelta(t, e, -2.1, 0.1)
	})
	t.Run("Bulge", func(t *testing.T) {
		// mock bulge of CAT on one side and AG on other
		// from pg 429 of SantaLucia, 2004
		seq := "ACCCCCATCCTTCCTTGAGTCAAGGGGCTCAA"
		foldContext, err := newFoldingContext(seq, 37)
		require.NoError(t, err)

		pairDg, err := Bulge(5, 7, 18, 17, foldContext)
		require.NoError(t, err)
		assert.InDelta(t, pairDg, 3.22, 0.4)
	})
	t.Run("hairpin", func(t *testing.T) {
		// hairpin = "CCTTGG"
		seq := "ACCCCCTCCTTCCTTGGATCAAGGGGCTCAA"
		i := 11
		j := 16

		foldContext, err := newFoldingContext(seq, 37)
		require.NoError(t, err)

		hairpinDg, err := hairpin(i, j, foldContext)
		require.NoError(t, err)
		// this differs from Unafold
		assert.InDelta(t, hairpinDg, 4.3, 1.0)

		// from page 428 of SantaLucia, 2004
		// hairpin = "CGCAAG"
		seq = "ACCCGCAAGCCCTCCTTCCTTGGATCAAGGGGCTCAA"
		i = 3
		j = 8

		foldContext, err = newFoldingContext(seq, 37)
		require.NoError(t, err)

		hairpinDg, err = hairpin(i, j, foldContext)
		require.NoError(t, err)
		assert.InDelta(t, hairpinDg, 0.67, 0.1)

		seq = "CUUUGCACG"
		i = 0
		j = 8

		foldContext, err = newFoldingContext(seq, 37)
		require.NoError(t, err)

		hairpinDg, err = hairpin(i, j, foldContext)
		require.NoError(t, err)
		assert.InDelta(t, hairpinDg, 4.5, 0.2)
	})
	t.Run("internalLoop", func(t *testing.T) {
		seq := "ACCCCCTCCTTCCTTGGATCAAGGGGCTCAA"
		i := 6
		j := 21

		foldContext, err := newFoldingContext(seq, 37)
		require.NoError(t, err)

		dg, err := internalLoop(i, i+4, j, j-4, foldContext)
		require.NoError(t, err)
		assert.InDelta(t, dg, 3.5, 0.1)
	})
	t.Run("W", func(t *testing.T) {
		seq := "GCUCAGCUGGGAGAGC"
		i := 0
		j := 15

		foldContext, err := newFoldingContext(seq, 37)
		require.NoError(t, err)

		struc, err := unpairedMinimumFreeEnergyW(i, j, foldContext)
		require.NoError(t, err)
		assert.InDelta(t, struc.energy, -3.8, 0.2)

		seq = "CCUGCUUUGCACGCAGG"
		i = 0
		j = 16

		foldContext, err = newFoldingContext(seq, 37)
		require.NoError(t, err)

		struc, err = unpairedMinimumFreeEnergyW(i, j, foldContext)
		require.NoError(t, err)
		assert.InDelta(t, struc.energy, -6.4, 0.2)

		seq = "GCGGUUCGAUCCCGC"
		i = 0
		j = 14

		foldContext, err = newFoldingContext(seq, 37)
		require.NoError(t, err)

		struc, err = unpairedMinimumFreeEnergyW(i, j, foldContext)
		require.NoError(t, err)
		assert.InDelta(t, struc.energy, -4.2, 0.2)
	})
}
func TestZuker_ErrorCreatingFoldingContext(t *testing.T) {
	seq := "ATGGATTTAGATAGATADFQ#(RSDOFIA)"
	temp := 4000.0

	expectedErr := fmt.Errorf("error creating folding context: the sequence ATGGATTTAGATAGATADFQ#(RSDOFIA) is not RNA or DNA")

	_, err := Zuker(seq, temp)
	require.Error(t, err)
	assert.Equal(t, expectedErr.Error(), err.Error())
}
