package fold

import (
	"fmt"
	"math"
	"strings"

	"github.com/TimothyStiles/poly/transform"
)

// Fold the DNA sequence and return the lowest free energy score.
//
// Based on the approach described in:
// Zuker and Stiegler, 1981
// https://www.ncbi.nlm.nih.gov/pmc/articles/PMC326673/pdf/nar00394-0137.pdf
//
// If the sequence is 50 or more bp long, "isolated" matching bp
// are ignored in V(i,j). This is based on an approach described in:
// Mathews, Sabina, Zuker and Turner, 1999
// https://www.ncbi.nlm.nih.gov/pubmed/10329189
// Args:
//
//	seq: The sequence to Fold
//	temp: The temperature the Fold takes place in, in Celsius
//
// Returns a slice of NucleicAcidStructure with the energy and description,
// i.e. stacks, bulges, hairpins, etc.
func Fold(seq string, temp float64) ([]NucleicAcidStructure, error) {
	fc, err := NewFoldingContext(seq, temp)
	if err != nil {
		return nil, fmt.Errorf("error creating folding context: %w", err)
	}

	// get the minimum free energy structure out of the cache
	return Traceback(0, len(seq)-1, fc), nil
}

// MinimumFreeEnergy folds the sequence and return just the delta G of the structure
// Args:
//
//	seq: The sequence to fold
//	temp: The temperature to fold at
//
// Returns: the minimum free energy of the folded sequence
func MinimumFreeEnergy(seq string, temp float64) (float64, error) {
	NucleicAcidStructures, err := Fold(seq, temp)
	if err != nil {
		return 0, fmt.Errorf("error folding: %w", err)
	}
	summedEnergy := 0.0
	for _, structure := range NucleicAcidStructures {
		summedEnergy += structure.E
	}
	return RoundFloat(summedEnergy, 2), nil
}

// Get the dot bracket notation for a secondary structure.
// Args:
//
//	NucleicAcidStructures: A list of NucleicAcidStructure, usually from the fold function
//
// # Returns the dot-bracket notation of the secondary structure
//
// Dot-bracket notation, consisting in a balanced parentheses string composed
// by a three-character alphabet {.,(,)}, that can be unambiguously converted
// in the RNA secondary structure. See example_test.go for a small example.
func DotBracket(NucleicAcidStructures []NucleicAcidStructure) string {
	maxj := 0
	for _, structure := range NucleicAcidStructures {
		for _, ij := range structure.Inner {
			if ij.End > maxj {
				maxj = ij.End
			}
		}
	}
	maxj += 1
	result := make([]byte, maxj)
	for i := range result {
		result[i] = '.'
	}
	for _, structure := range NucleicAcidStructures {
		if len(structure.Inner) == 1 {
			ij := structure.Inner[0]
			result[ij.Start] = '('
			result[ij.End] = ')'
		}
	}
	return string(result)
}

// Find and return the lowest free energy structure in the subsequence starting
// at i and terminating at j.
//
// Figure 2B in Zuker and Stiegler, 1981
// Args:
//
//		seq: The sequence being folded
//		i: The start index
//		j: The end index (inclusive)
//	 fc: The FoldingContext for this sequence
//
// Returns the free energy for the subsequence from i to j
func W(i, j int, fc FoldContext) (NucleicAcidStructure, error) {
	if !fc.W[i][j].Equal(STRUCT_DEFAULT) {
		return fc.W[i][j], nil
	}

	if j-i < 4 {
		fc.W[i][j] = STRUCT_NULL
		return fc.W[i][j], nil
	}

	w1, err := W(i+1, j, fc)
	if err != nil {
		return STRUCT_DEFAULT, fmt.Errorf("w: subsequence (%d, %d): %w", i, j, err)
	}
	w2, err := W(i, j-1, fc)
	if err != nil {
		return STRUCT_DEFAULT, fmt.Errorf("w: subsequence (%d, %d): %w", i, j, err)
	}
	w3, err := V(i, j, fc)
	if err != nil {
		return STRUCT_DEFAULT, fmt.Errorf("w: subsequence (%d, %d): %w", i, j, err)
	}

	w4 := STRUCT_NULL
	for k := i + 1; k < j-1; k++ {
		w4_test, err := MultiBranch(i, k, j, fc, false)
		if err != nil {
			return STRUCT_DEFAULT, fmt.Errorf("w: subsequence (%d, %d): %w", i, j, err)
		}

		if w4_test.Valid() && w4_test.E < w4.E {
			w4 = w4_test
		}
	}

	wret := minStruct(w1, w2, w3, w4)
	fc.W[i][j] = wret
	return wret, nil
}

// Find, store and return the minimum free energy of the structure between i
// and j.
//
// If i and j don't bp, store and return INF.
// See: Figure 2B of Zuker, 1981
// Args:
//
//		i: The start index
//		j: The end index (inclusive)
//	 fc: The FoldingContext for this sequence
//
// Returns the minimum energy folding structure possible between i and j on seq
func V(i, j int, fc FoldContext) (NucleicAcidStructure, error) {
	if !fc.V[i][j].Equal(STRUCT_DEFAULT) {
		return fc.V[i][j], nil
	}

	// the ends must basepair for V(i,j)
	if fc.Energies.COMPLEMENT[fc.Seq[i]] != fc.Seq[j] {
		fc.V[i][j] = STRUCT_NULL
		return fc.V[i][j], nil
	}
	// if the basepair is isolated, and the seq large, penalize at 1,600 kcal/mol
	// heuristic for speeding this up
	// from https://www.ncbi.nlm.nih.gov/pubmed/10329189
	isolated_outer := true
	if i > 0 && j < len(fc.Seq)-1 {
		isolated_outer = fc.Energies.COMPLEMENT[fc.Seq[i-1]] != fc.Seq[j+1]
	}
	isolated_inner := fc.Energies.COMPLEMENT[fc.Seq[i+1]] != fc.Seq[j-1]

	if isolated_outer && isolated_inner {
		fc.V[i][j] = NucleicAcidStructure{E: 1600}
		return fc.V[i][j], nil
	}

	p := Pair(fc.Seq, i, i+1, j, j-1)
	hp, err := Hairpin(i, j, fc)
	if err != nil {
		return STRUCT_DEFAULT, fmt.Errorf("v: subsequence (%d, %d): %w", i, j, err)
	}
	e1 := NucleicAcidStructure{E: hp, Desc: "HAIRPIN:" + p}
	if j-i == 4 { // small hairpin; 4bp
		fc.V[i][j] = e1
		fc.W[i][j] = e1
		return fc.V[i][j], nil
	}

	n := len(fc.Seq)
	e2 := NucleicAcidStructure{E: math.Inf(1)}
	for i1 := i + 1; i1 < j-4; i1++ {
		for j1 := i1 + 4; j1 < j; j1++ {
			// i1 and j1 must match
			if fc.Energies.COMPLEMENT[fc.Seq[i1]] != fc.Seq[j1] {
				continue
			}

			p := Pair(fc.Seq, i, i1, j, j1)
			pair_left := Pair(fc.Seq, i, i+1, j, j-1)
			pair_right := Pair(fc.Seq, i1-1, i1, j1+1, j1)
			_, plin := fc.Energies.NN[pair_left]
			_, prin := fc.Energies.NN[pair_right]
			pair_inner := plin || prin

			stck := i1 == i+1 && j1 == j-1
			bulge_left := i1 > i+1
			bulge_right := j1 < j-1

			var (
				e2_test      float64
				e2_test_type string
				err          error
			)
			switch {
			case stck:
				// it's a neighboring/stacking pair in a helix
				e2_test = Stack(i, i1, j, j1, fc)
				e2_test_type = fmt.Sprintf("STACK:%s", p)

				if i > 0 && j == n-1 || i == 0 && j < n-1 {
					// there's a dangling end
					e2_test_type = fmt.Sprintf("STACK_DE:%s", p)
				}
			case bulge_left && bulge_right && !pair_inner:
				// it's an interior loop
				il, err := InternalLoop(i, i1, j, j1, fc)
				if err != nil {
					return STRUCT_DEFAULT, fmt.Errorf("v: subsequence (%d, %d): %w", i, j, err)
				}
				e2_test = il
				e2_test_type = fmt.Sprintf("INTERIOR_LOOP:%d/%d", i1-i, j-j1)

				if i1-i == 2 && j-j1 == 2 {
					loop_left := fc.Seq[i : i1+1]
					loop_right := fc.Seq[j1 : j+1]
					// technically an interior loop of 1. really 1bp mismatch
					e2_test_type = fmt.Sprintf("STACK:%s/%s", loop_left, transform.Reverse(loop_right))
				}
			case bulge_left && !bulge_right:
				// it's a bulge on the left side
				e2_test, err = Bulge(i, i1, j, j1, fc)
				if err != nil {
					return STRUCT_DEFAULT, fmt.Errorf("v: subsequence (%d, %d): %w", i, j, err)
				}
				e2_test_type = fmt.Sprintf("BULGE:%d", i1-i)
			case !bulge_left && bulge_right:
				// it's a bulge on the right side
				e2_test, err = Bulge(i, i1, j, j1, fc)
				if err != nil {
					return STRUCT_DEFAULT, fmt.Errorf("v: subsequence (%d, %d): %w", i, j, err)
				}
				e2_test_type = fmt.Sprintf("BULGE:%d", j-j1)
			default:
				// it's basically a hairpin, only outside bp match
				continue
			}

			// add V(i', j')
			tv, err := V(i1, j1, fc)
			if err != nil {
				return STRUCT_DEFAULT, fmt.Errorf("v: subsequence (%d, %d): %w", i, j, err)
			}
			e2_test += tv.E
			if e2_test != math.Inf(-1) && e2_test < e2.E {
				e2 = NucleicAcidStructure{E: e2_test, Desc: e2_test_type, Inner: []Subsequence{{i1, j1}}}
			}
		}
	}

	e3 := STRUCT_NULL
	if !isolated_outer || i == 0 || j == len(fc.Seq)-1 {
		for k := i + 1; k < j-1; k++ {
			e3_test, err := MultiBranch(i, k, j, fc, true)
			if err != nil {
				return STRUCT_DEFAULT, fmt.Errorf("v: subsequence (%d, %d): %w", i, j, err)
			}

			if e3_test.Valid() && e3_test.E < e3.E {
				e3 = e3_test
			}
		}
	}
	e := minStruct(e1, e2, e3)
	fc.V[i][j] = e
	return e, nil
}

// Bulge calculates the free energy associated with a bulge.
//
// Args:
//
//		i: The start index of the bulge
//		i1: The index to the right of i
//		j: The end index of the bulge
//		j1: The index to the left of j
//	 fc: The FoldingContext for this sequence
//
// Returns the increment in free energy from the bulge
func Bulge(i, i1, j, j1 int, fc FoldContext) (float64, error) {
	loop_len := max(i1-i-1, j-j1-1)
	if loop_len <= 0 {
		return 0, fmt.Errorf("bulge: the length of the bulge at (%d, %d) is %d", i, j, loop_len)
	}

	var dG float64

	// add penalty based on size
	if en, ok := fc.Energies.BULGE_LOOPS[loop_len]; ok {
		d_h, d_s := en.H, en.S
		dG = DeltaG(d_h, d_s, fc.T)
	} else {
		// it's too large for pre-calculated list, extrapolate
		en := fc.Energies.BULGE_LOOPS[30]
		d_h, d_s := en.H, en.S
		dG = DeltaG(d_h, d_s, fc.T)
		dG = JacobsonStockmayer(loop_len, 30, dG, fc.T)
	}

	if loop_len == 1 {
		// if len 1, include the delta G of intervening NN (SantaLucia 2004)
		pair := Pair(fc.Seq, i, i1, j, j1)
		if _, ok := fc.Energies.NN[pair]; !ok {
			return 0, fmt.Errorf("bulge: pair %q not in the NN energies", pair)
		}
		dG += Stack(i, i1, j, j1, fc)
	}

	// penalize AT terminal bonds
	for _, k := range []int{i, i1, j, j1} {
		if fc.Seq[k] == 'A' {
			dG += 0.5
		}
	}

	return dG, nil
}

func addBranch(structure NucleicAcidStructure, branches *[]Subsequence, fc FoldContext) error {
	if !structure.Valid() || len(structure.Inner) == 0 {
		return nil
	}
	if len(structure.Inner) == 1 {
		*branches = append(*branches, structure.Inner[0])
		return nil
	}
	for _, subsq := range structure.Inner {
		str, err := W(subsq.Start, subsq.End, fc)
		if err != nil {
			return err
		}
		err = addBranch(str, branches, fc)
		if err != nil {
			return err
		}
	}
	return nil
}

// MultiBranch calculats a multi-branch energy penalty using a linear formula.
//
// From Jaeger, Turner, and Zuker, 1989.
// Found to be better than logarithmic in Ward, et al. 2017
// Args:
//
//		i: The left starting index
//		k: The mid-point in the search
//		j: The right ending index
//	 fc: The FoldingContext for this sequence
//		helix: Whether this multibranch is enclosed by a helix
//		helix: Whether V(i, j) bond with one another in a helix
//
// Returns a multi-branch structure
func MultiBranch(i, k, j int, fc FoldContext, helix bool) (NucleicAcidStructure, error) {
	var (
		left, right NucleicAcidStructure
		err         error
	)
	if helix {
		left, err = W(i+1, k, fc)
		if err != nil {
			return STRUCT_DEFAULT, fmt.Errorf("multibranch: subsequence (%d, %d, %d): %w", i, j, k, err)
		}
		right, err = W(k+1, j-1, fc)
		if err != nil {
			return STRUCT_DEFAULT, fmt.Errorf("multibranch: subsequence (%d, %d, %d): %w", i, j, k, err)
		}
	} else {
		left, err = W(i, k, fc)
		if err != nil {
			return STRUCT_DEFAULT, fmt.Errorf("multibranch: subsequence (%d, %d, %d): %w", i, j, k, err)
		}
		right, err = W(k+1, j, fc)
		if err != nil {
			return STRUCT_DEFAULT, fmt.Errorf("multibranch: subsequence (%d, %d, %d): %w", i, j, k, err)
		}
	}

	if !left.Valid() || !right.Valid() {
		return STRUCT_NULL, nil
	}

	// gather all branches of this multi-branch structure
	var branches []Subsequence

	// in python this was a recursive closure, in Go this is not possible so
	// we pull it out and pass all the parameters
	err = addBranch(left, &branches, fc)
	if err != nil {
		return STRUCT_DEFAULT, fmt.Errorf("multibranch: subsequence (%d, %d, %d): %w", i, j, k, err)
	}
	err = addBranch(right, &branches, fc)
	if err != nil {
		return STRUCT_DEFAULT, fmt.Errorf("multibranch: subsequence (%d, %d, %d): %w", i, j, k, err)
	}

	// this isn't multi-branched
	if len(branches) < 2 {
		return STRUCT_NULL, nil
	}

	// if there's a helix, i,j counts as well
	if helix {
		branches = append(branches, Subsequence{i, j})
	}

	// count up unpaired bp and asymmetry
	branches_count := len(branches)
	unpaired := 0
	e_sum := 0.0
	ij := Subsequence{i, j}
	for index, ij2 := range branches {
		i2, j2 := ij2.Start, ij2.End
		ij1 := branches[abs((index-1)%len(branches))]
		j1 := ij1.End
		ij3 := branches[abs((index+1)%len(branches))]
		i3, j3 := ij3.Start, ij3.End

		// add energy from unpaired bp to the right
		// of the helix as though it was a dangling end
		// if there's only one bp, it goes to whichever
		// helix (this or the next) has the more favorable energy
		unpaired_left := 0
		unpaired_right := 0
		de := 0.0
		if index == len(branches)-1 && !helix {
			// pass
		} else if ij3 == ij {
			unpaired_left = i2 - j1 - 1
			unpaired_right = j3 - j2 - 1

			if unpaired_left != 0 && unpaired_right != 0 {
				de = Stack(i2-1, i2, j2+1, j2, fc)
			} else if unpaired_right != 0 {
				de = Stack(-1, i2, j2+1, j2, fc)
				if unpaired_right == 1 {
					de = min(Stack(i3, -1, j3, j3-1, fc), de)
				}
			}
		} else if ij2 == ij {
			unpaired_left = j2 - j1 - 1
			unpaired_right = i3 - i2 - 1

			if unpaired_left != 0 && unpaired_right != 0 {
				de = Stack(i2-1, i2, j2+1, j2, fc)
			} else if unpaired_right != 0 {
				de = Stack(i2, i2+1, j2, -1, fc)
				if unpaired_right == 1 {
					de = min(Stack(i3-1, i3, -1, j3, fc), de)
				}
			}
		} else {
			unpaired_left = i2 - j1 - 1
			unpaired_right = i3 - j2 - 1

			if unpaired_left != 0 && unpaired_right != 0 {
				de = Stack(i2-1, i2, j2+1, j2, fc)
			} else if unpaired_right != 0 {
				de = Stack(-1, i2, j2+1, j2, fc)
				if unpaired_right == 1 {
					de = min(Stack(i2-1, i2, j2+1, j2, fc), de)
				}
			}
		}
		e_sum += de
		unpaired += unpaired_right
		if unpaired_right < 0 {
			return STRUCT_DEFAULT, fmt.Errorf("multibranch: subsequence (%d, %d, %d): unpaired_right < 0", i, j, k)
		}

		if ij2 != ij { // add energy
			w, err := W(i2, j2, fc)
			if err != nil {
				return STRUCT_DEFAULT, fmt.Errorf("multibranch: subsequence (%d, %d, %d): %w", i, j, k, err)
			}
			e_sum += w.E
		}

	}

	if unpaired < 0 {
		return STRUCT_DEFAULT, fmt.Errorf("multibranch: subsequence (%d, %d, %d): unpaired < 0", i, j, k)
	}

	// penalty for unmatched bp and multi-branch
	a, b, c, d := fc.Energies.MULTIBRANCH.A, fc.Energies.MULTIBRANCH.B, fc.Energies.MULTIBRANCH.C, fc.Energies.MULTIBRANCH.D
	e_multibranch := a + b*float64(len(branches)) + c*float64(unpaired)

	if unpaired == 0 {
		e_multibranch = a + d
	}

	// energy of min-energy neighbors
	e := e_multibranch + e_sum

	// pointer to next structures
	if helix {
		// branches.pop()
		branches = branches[:len(branches)-1]
	}

	return NucleicAcidStructure{E: e, Desc: fmt.Sprintf("BIFURCATION:%dn/%dh", unpaired, branches_count), Inner: branches}, nil
}

// InternalLoop calculates the free energy of an internal loop.
//
// The first and last bp of both left and right sequences
// are not themselves parts of the loop, but are the terminal
// bp on either side of it. They are needed for when there's
// a single internal looping bp (where just the mismatching
// free energies are used)
// Note that both left and right sequences are in 5' to 3' direction
// This is adapted from the "Internal Loops" section of SantaLucia/Hicks, 2004
// Args:
//
//		i:  The index of the start of structure on left side
//		i1: The index to the right of i
//		j:  The index of the end of structure on right side
//		j1: The index to the left of j
//	 fc: The FoldingContext for this sequence
//
// Returns the free energy associated with the internal loop
func InternalLoop(i, i1, j, j1 int, fc FoldContext) (float64, error) {
	loop_left := i1 - i - 1
	loop_right := j - j1 - 1
	loop_len := loop_left + loop_right

	if loop_left < 1 {
		return 0, fmt.Errorf("internal loop: subsequence (%d, %d, %d, %d): missing left part of the loop", i, i1, j, j1)

	}
	if loop_right < 1 {
		return 0, fmt.Errorf("internal loop: subsequence (%d, %d, %d, %d): missing right part of the loop", i, i1, j, j1)
	}

	// single bp mismatch, sum up the two single mismatch pairs
	if loop_left == 1 && loop_right == 1 {
		mm_left := Stack(i, i1, j, j1, fc)
		mm_right := Stack(i1-1, i1, j1+1, j1, fc)
		return mm_left + mm_right, nil
	}
	var d_h, d_s, dG float64
	// apply a penalty based on loop size
	if en, ok := fc.Energies.INTERNAL_LOOPS[loop_len]; ok {
		d_h, d_s = en.H, en.S
		dG = DeltaG(d_h, d_s, fc.T)
	} else {
		// it's too large an internal loop, extrapolate
		en := fc.Energies.INTERNAL_LOOPS[30]
		d_h, d_s = en.H, en.S
		dG = DeltaG(d_h, d_s, fc.T)
		dG = JacobsonStockmayer(loop_len, 30, dG, fc.T)
	}

	// apply an asymmetry penalty
	loop_asymmetry := math.Abs(float64(loop_left - loop_right))
	dG += 0.3 * loop_asymmetry

	// apply penalty based on the mismatching pairs on either side of the loop
	pair_left_mm := Pair(fc.Seq, i, i+1, j, j-1)
	en := fc.Energies.TERMINAL_MM[pair_left_mm]
	d_h, d_s = en.H, en.S
	dG += DeltaG(d_h, d_s, fc.T)

	pair_right_mm := Pair(fc.Seq, i1-1, i1, j1+1, j1)
	en = fc.Energies.TERMINAL_MM[pair_right_mm]
	d_h, d_s = en.H, en.S
	dG += DeltaG(d_h, d_s, fc.T)

	return dG, nil
}

// Stack return the free energy of a stack
//
// Using the indexes i and j, check whether it's at the end of
// the sequence or internal. Then check whether it's a match
// or mismatch, and return.
// Two edge-cases are terminal mismatches and dangling ends.
// The energy of a dangling end is added to the energy of a pair
// where i XOR j is at the sequence's end.
// Args:
//
//		i: The start index on left side of the pair/stack
//		i1: The index to the right of i
//		j: The end index on right side of the pair/stack
//		j1: The index to the left of j
//	 fc: The FoldingContext for this sequence
//
// Returns the free energy of the NN pairing
func Stack(i, i1, j, j1 int, fc FoldContext) float64 {
	// if any(x >= len(seq) for x in [i, i1, j, j1]):
	//    return 0.0
	for _, x := range []int{i, i1, j, j1} {
		if x >= len(fc.Seq) {
			return 0
		}
	}

	pair := Pair(fc.Seq, i, i1, j, j1)
	// if any(x == -1 for x in [i, i1, j, j1]):
	for _, x := range []int{i, i1, j, j1} {
		if x == -1 {
			// it's a dangling end
			en := fc.Energies.DE[pair]
			d_h, d_s := en.H, en.S
			return DeltaG(d_h, d_s, fc.T)
		}
	}

	if i > 0 && j < len(fc.Seq)-1 {
		// it's internal
		en, ok := fc.Energies.NN[pair]
		if !ok {
			en = fc.Energies.INTERNAL_MM[pair]
		}
		d_h, d_s := en.H, en.S
		return DeltaG(d_h, d_s, fc.T)
	}
	if i == 0 && j == len(fc.Seq)-1 {
		// it's terminal
		en, ok := fc.Energies.NN[pair]
		if !ok {
			en = fc.Energies.INTERNAL_MM[pair]
		}
		d_h, d_s := en.H, en.S
		return DeltaG(d_h, d_s, fc.T)
	}

	if i > 0 && j == len(fc.Seq)-1 {
		// it's dangling on left
		en, ok := fc.Energies.NN[pair]
		if !ok {
			en = fc.Energies.INTERNAL_MM[pair]
		}
		d_h, d_s := en.H, en.S
		dG := DeltaG(d_h, d_s, fc.T)

		pair_de := fmt.Sprintf("%c%c/.%c", fc.Seq[i-1], fc.Seq[i], fc.Seq[j])
		if en, ok := fc.Energies.DE[pair_de]; ok {
			d_h, d_s := en.H, en.S
			dG += DeltaG(d_h, d_s, fc.T)
		}
		return dG
	}

	if i == 0 && j < len(fc.Seq)-1 {
		// it's dangling on right
		en, ok := fc.Energies.NN[pair]
		if !ok {
			en = fc.Energies.INTERNAL_MM[pair]
		}
		d_h, d_s := en.H, en.S
		dG := DeltaG(d_h, d_s, fc.T)

		pair_de := fmt.Sprintf(".%c/%c%c", +fc.Seq[i], fc.Seq[j+1], fc.Seq[j])
		if en, ok := fc.Energies.DE[pair_de]; ok {
			d_h, d_s := en.H, en.S
			dG += DeltaG(d_h, d_s, fc.T)
			return dG
		}
	}
	return 0
}

// Hairpin calculates the free energy of a hairpin.
// Args:
//
//		i:  The index of start of hairpin
//		j:  The index of end of hairpin
//	 fc: The FoldingContext for this sequence
//
// Returns the free energy increment from the hairpin structure
func Hairpin(i, j int, fc FoldContext) (float64, error) {
	if j-i < 4 {
		return math.Inf(1), nil
	}

	hairpinSeq := fc.Seq[i : j+1]
	hairpin_len := len(hairpinSeq) - 2
	pair := Pair(fc.Seq, i, i+1, j, j-1)

	if fc.Energies.COMPLEMENT[hairpinSeq[0]] != hairpinSeq[len(hairpinSeq)-1] {
		// not known terminal pair, nothing to close "hairpin"
		return 0, fmt.Errorf("hairpin: subsequence (%d, %d): unknown hairpin terminal pairing %c - %c", i, j, hairpinSeq[0], hairpinSeq[len(hairpinSeq)-1])
	}

	dG := 0.0
	if fc.Energies.TRI_TETRA_LOOPS != nil {
		if en, ok := fc.Energies.TRI_TETRA_LOOPS[hairpinSeq]; ok {
			// it's a pre-known hairpin with known value
			d_h, d_s := en.H, en.S
			dG = DeltaG(d_h, d_s, fc.T)
		}
	}

	// add penalty based on size
	if en, ok := fc.Energies.HAIRPIN_LOOPS[hairpin_len]; ok {
		d_h, d_s := en.H, en.S
		dG += DeltaG(d_h, d_s, fc.T)
	} else {
		// it's too large, extrapolate
		en := fc.Energies.HAIRPIN_LOOPS[30]
		d_h, d_s := en.H, en.S
		d_g_inc := DeltaG(d_h, d_s, fc.T)
		dG += JacobsonStockmayer(hairpin_len, 30, d_g_inc, fc.T)
	}

	// add penalty for a terminal mismatch
	en, ok := fc.Energies.TERMINAL_MM[pair]
	if hairpin_len > 3 && ok {
		d_h, d_s := en.H, en.S
		dG += DeltaG(d_h, d_s, fc.T)
	}

	// add penalty if length 3 and AT closing, formula 8 from SantaLucia, 2004
	if hairpin_len == 3 && (hairpinSeq[0] == 'A' || hairpinSeq[len(hairpinSeq)-1] == 'A') {
		dG += 0.5 // convert to entropy
	}

	return dG, nil
}

// Find the free energy given delta h, s and temp
// Args:
//
//	d_h: The enthalpy increment in kcal / mol
//	d_s: The entropy increment in cal / mol
//	temp: The temperature in Kelvin
//
// Returns the free energy increment in kcal / (mol x K)
func DeltaG(d_h, d_s, temp float64) float64 {
	return d_h - temp*(d_s/1000.0)
}

// Estimate the free energy of length query_len based on one of length known_len.
//
// The Jacobson-Stockmayer entry extrapolation formula is used
// for bulges, hairpins, etc that fall outside the 30nt upper limit
// for pre-calculated free-energies. See SantaLucia and Hicks (2004).
// Args:
//
//	query_len: Length of element without known free energy value
//	known_len: Length of element with known free energy value (d_g_x)
//	d_g_x: The free energy of the element known_len
//	temp: Temperature in Kelvin
//
// Returns the free energy for a structure of length query_len
func JacobsonStockmayer(query_len, known_len int, d_g_x, temp float64) float64 {
	gas_constant := 1.9872e-3
	return d_g_x + 2.44*gas_constant*temp*math.Log(float64(query_len)/float64(known_len))
}

// Pair Returns a stack representation, a key for the NN maps
// Args:
//
//	s: Sequence being folded
//	i: leftmost index
//	i1: index to right of i
//	j: rightmost index
//	j1: index to left of j
//
// Returns string representation of the Pair
func Pair(s string, i, i1, j, j1 int) string {
	ss := []rune(s)
	ret := []rune{'.', '.', '/', '.', '.'}
	if i >= 0 {
		ret[0] = ss[i]
	}
	if i1 >= 0 {
		ret[1] = ss[i1]
	}
	if j >= 0 {
		ret[3] = ss[j]
	}
	if j1 >= 0 {
		ret[4] = ss[j1]
	}
	return string(ret)
}

// Traceback thru the V(i,j) and W(i,j) caches to find the structure
// For each step, get to the lowest energy W(i,j) within that block
// Store the structure in W(i,j)
// Inc i and j
// If the next structure is viable according to V(i,j), store as well
// Repeat
// Args:
//
//		i: The leftmost index to start searching in
//		j: The rightmost index to start searching in
//	 fc: The FoldingContext for this sequence
//
// Returns a list of NucleicAcidStructure in the final secondary structure
func Traceback(i, j int, fc FoldContext) []NucleicAcidStructure {
	// move i,j down-left to start coordinates
	s := fc.W[i][j]
	if !strings.Contains(s.Desc, "HAIRPIN") {
		for fc.W[i+1][j].Equal(s) {
			i += 1
		}
		for fc.W[i][j-1].Equal(s) {
			j -= 1
		}
	}

	NucleicAcidStructures := []NucleicAcidStructure{}
	for {
		s = fc.V[i][j]

		NucleicAcidStructures = append(NucleicAcidStructures, NucleicAcidStructure{E: s.E, Desc: s.Desc, Inner: []Subsequence{{Start: i, End: j}}})

		// it's a hairpin, end of structure
		if len(s.Inner) == 0 {
			// set the energy of everything relative to the hairpin
			return trackbackEnergy(NucleicAcidStructures)
		}

		// it's a stack, bulge, etc
		// there's another single structure beyond this
		if len(s.Inner) == 1 {
			i, j = s.Inner[0].Start, s.Inner[0].End
			// ij = s.IJs[0]
			continue
		}

		// it's a multibranch
		e_sum := 0.0
		NucleicAcidStructures = trackbackEnergy(NucleicAcidStructures)
		branches := []NucleicAcidStructure{}
		for _, ij1 := range s.Inner {
			i1, j1 := ij1.Start, ij1.End
			tb := Traceback(i1, j1, fc)
			if len(tb) > 0 && len(tb[0].Inner) > 0 {
				ij2 := tb[0].Inner[0]
				i2, j2 := ij2.Start, ij2.End
				e_sum += fc.W[i2][j2].E
				branches = append(branches, tb...)
			}
		}

		NucleicAcidStructures[len(NucleicAcidStructures)-1].E -= e_sum
		return append(NucleicAcidStructures, branches...)
	}
}

// Return the struct with the lowest free energy that isn't -inf
// Args:
//
//	structures: NucleicAcidStructure being compared
//
// Returns the min free energy structure
func minStruct(structures ...NucleicAcidStructure) NucleicAcidStructure {
	minimumStructure := STRUCT_NULL
	for _, str := range structures {
		if str.E != math.Inf(-1) && str.E < minimumStructure.E {
			minimumStructure = str
		}
	}
	return minimumStructure
}

// trackbackEnergy add energy to each structure, based on how it's
// W(i,j) differs from the one after
// Args:
//
//	structures: The NucleicAcidStructure for whom energy is being calculated
//
// Returns a slice of NucleicAcidStructure in the folded DNA with energy
func trackbackEnergy(structures []NucleicAcidStructure) []NucleicAcidStructure {
	for i := 0; i < len(structures)-1; i++ {
		structures[i].E -= structures[i+1].E
	}
	return structures
}
