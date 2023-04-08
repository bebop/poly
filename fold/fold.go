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
// are ignored in V(start,end). This is based on an approach described in:
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
	foldContext, err := NewFoldingContext(seq, temp)
	if err != nil {
		return nil, fmt.Errorf("error creating folding context: %w", err)
	}

	// get the minimum free energy structure out of the cache
	return Traceback(0, len(seq)-1, foldContext), nil
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
		summedEnergy += structure.Energy
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
// at start and terminating at end.
//
// Figure 2B in Zuker and Stiegler, 1981
// Args:
//
//		seq: The sequence being folded
//		start: The start index
//		end: The end index (inclusive)
//	 foldContext: The FoldContext for this sequence
//
// Returns the free energy for the subsequence from start to end
func W(start, end int, foldContext FoldContext) (NucleicAcidStructure, error) {
	if !foldContext.W[start][end].Equal(StructDefault) {
		return foldContext.W[start][end], nil
	}

	if end-start < 4 {
		foldContext.W[start][end] = StructInvalid
		return foldContext.W[start][end], nil
	}

	w1, err := W(start+1, end, foldContext)
	if err != nil {
		return StructDefault, fmt.Errorf("w: subsequence (%d, %d): %w", start, end, err)
	}
	w2, err := W(start, end-1, foldContext)
	if err != nil {
		return StructDefault, fmt.Errorf("w: subsequence (%d, %d): %w", start, end, err)
	}
	w3, err := V(start, end, foldContext)
	if err != nil {
		return StructDefault, fmt.Errorf("w: subsequence (%d, %d): %w", start, end, err)
	}

	w4 := StructInvalid
	for k := start + 1; k < end-1; k++ {
		w4_test, err := Multibranch(start, k, end, foldContext, false)
		if err != nil {
			return StructDefault, fmt.Errorf("w: subsequence (%d, %d): %w", start, end, err)
		}

		if w4_test.Valid() && w4_test.Energy < w4.Energy {
			w4 = w4_test
		}
	}

	wret := minStruct(w1, w2, w3, w4)
	foldContext.W[start][end] = wret
	return wret, nil
}

// Find, store and return the minimum free energy of the structure between start
// and end.
//
// If start and end don't bp, store and return INF.
// See: Figure 2B of Zuker, 1981
// Args:
//
//		start: The start index
//		end: The end index (inclusive)
//	 foldContext: The FoldingContext for this sequence
//
// Returns the minimum energy folding structure possible between start and end on seq
func V(start, end int, foldContext FoldContext) (NucleicAcidStructure, error) {
	if !foldContext.V[start][end].Equal(StructDefault) {
		return foldContext.V[start][end], nil
	}

	// the ends must basepair for V(start,end)
	if foldContext.Energies.Complement[foldContext.Seq[start]] != foldContext.Seq[end] {
		foldContext.V[start][end] = StructInvalid
		return foldContext.V[start][end], nil
	}
	// if the basepair is isolated, and the seq large, penalize at 1,600 kcal/mol
	// heuristic for speeding this up
	// from https://www.ncbi.nlm.nih.gov/pubmed/10329189
	isolated_outer := true
	if start > 0 && end < len(foldContext.Seq)-1 {
		isolated_outer = foldContext.Energies.Complement[foldContext.Seq[start-1]] != foldContext.Seq[end+1]
	}
	isolated_inner := foldContext.Energies.Complement[foldContext.Seq[start+1]] != foldContext.Seq[end-1]

	if isolated_outer && isolated_inner {
		foldContext.V[start][end] = NucleicAcidStructure{Energy: 1600}
		return foldContext.V[start][end], nil
	}

	p := Pair(foldContext.Seq, start, start+1, end, end-1)
	hp, err := Hairpin(start, end, foldContext)
	if err != nil {
		return StructDefault, fmt.Errorf("v: subsequence (%d, %d): %w", start, end, err)
	}
	e1 := NucleicAcidStructure{Energy: hp, Desc: "HAIRPIN:" + p}
	if end-start == 4 { // small hairpin; 4bp
		foldContext.V[start][end] = e1
		foldContext.W[start][end] = e1
		return foldContext.V[start][end], nil
	}

	n := len(foldContext.Seq)
	e2 := NucleicAcidStructure{Energy: math.Inf(1)}
	for i1 := start + 1; i1 < end-4; i1++ {
		for j1 := i1 + 4; j1 < end; j1++ {
			// i1 and j1 must match
			if foldContext.Energies.Complement[foldContext.Seq[i1]] != foldContext.Seq[j1] {
				continue
			}

			p := Pair(foldContext.Seq, start, i1, end, j1)
			pair_left := Pair(foldContext.Seq, start, start+1, end, end-1)
			pair_right := Pair(foldContext.Seq, i1-1, i1, j1+1, j1)
			_, plin := foldContext.Energies.NN[pair_left]
			_, prin := foldContext.Energies.NN[pair_right]
			pair_inner := plin || prin

			stck := i1 == start+1 && j1 == end-1
			bulge_left := i1 > start+1
			bulge_right := j1 < end-1

			var (
				e2_test      float64
				e2_test_type string
				err          error
			)
			switch {
			case stck:
				// it's a neighboring/stacking pair in a helix
				e2_test = Stack(start, i1, end, j1, foldContext)
				e2_test_type = fmt.Sprintf("STACK:%s", p)

				if start > 0 && end == n-1 || start == 0 && end < n-1 {
					// there's a dangling end
					e2_test_type = fmt.Sprintf("STACK_DE:%s", p)
				}
			case bulge_left && bulge_right && !pair_inner:
				// it's an interior loop
				il, err := InternalLoop(start, i1, end, j1, foldContext)
				if err != nil {
					return StructDefault, fmt.Errorf("v: subsequence (%d, %d): %w", start, end, err)
				}
				e2_test = il
				e2_test_type = fmt.Sprintf("INTERIOR_LOOP:%d/%d", i1-start, end-j1)

				if i1-start == 2 && end-j1 == 2 {
					loop_left := foldContext.Seq[start : i1+1]
					loop_right := foldContext.Seq[j1 : end+1]
					// technically an interior loop of 1. really 1bp mismatch
					e2_test_type = fmt.Sprintf("STACK:%s/%s", loop_left, transform.Reverse(loop_right))
				}
			case bulge_left && !bulge_right:
				// it's a bulge on the left side
				e2_test, err = Bulge(start, i1, end, j1, foldContext)
				if err != nil {
					return StructDefault, fmt.Errorf("v: subsequence (%d, %d): %w", start, end, err)
				}
				e2_test_type = fmt.Sprintf("BULGE:%d", i1-start)
			case !bulge_left && bulge_right:
				// it's a bulge on the right side
				e2_test, err = Bulge(start, i1, end, j1, foldContext)
				if err != nil {
					return StructDefault, fmt.Errorf("v: subsequence (%d, %d): %w", start, end, err)
				}
				e2_test_type = fmt.Sprintf("BULGE:%d", end-j1)
			default:
				// it's basically a hairpin, only outside bp match
				continue
			}

			// add V(start', end')
			tv, err := V(i1, j1, foldContext)
			if err != nil {
				return StructDefault, fmt.Errorf("v: subsequence (%d, %d): %w", start, end, err)
			}
			e2_test += tv.Energy
			if e2_test != math.Inf(-1) && e2_test < e2.Energy {
				e2 = NucleicAcidStructure{Energy: e2_test, Desc: e2_test_type, Inner: []Subsequence{{i1, j1}}}
			}
		}
	}

	e3 := StructInvalid
	if !isolated_outer || start == 0 || end == len(foldContext.Seq)-1 {
		for k := start + 1; k < end-1; k++ {
			e3_test, err := Multibranch(start, k, end, foldContext, true)
			if err != nil {
				return StructDefault, fmt.Errorf("v: subsequence (%d, %d): %w", start, end, err)
			}

			if e3_test.Valid() && e3_test.Energy < e3.Energy {
				e3 = e3_test
			}
		}
	}
	e := minStruct(e1, e2, e3)
	foldContext.V[start][end] = e
	return e, nil
}

// Bulge calculates the free energy associated with a bulge.
//
// Args:
//
//		start: The start index of the bulge
//		i1: The index to the right of start
//		end: The end index of the bulge
//		j1: The index to the left of end
//	 foldContext: The FoldingContext for this sequence
//
// Returns the increment in free energy from the bulge
func Bulge(start, i1, end, j1 int, foldContext FoldContext) (float64, error) {
	loop_len := max(i1-start-1, end-j1-1)
	if loop_len <= 0 {
		return 0, fmt.Errorf("bulge: the length of the bulge at (%d, %d) is %d", start, end, loop_len)
	}

	var dG float64

	// add penalty based on size
	if en, ok := foldContext.Energies.BulgeLoops[loop_len]; ok {
		d_h, d_s := en.EnthalpyH, en.EntropyS
		dG = DeltaG(d_h, d_s, foldContext.T)
	} else {
		// it's too large for pre-calculated list, extrapolate
		en := foldContext.Energies.BulgeLoops[30]
		d_h, d_s := en.EnthalpyH, en.EntropyS
		dG = DeltaG(d_h, d_s, foldContext.T)
		dG = JacobsonStockmayer(loop_len, 30, dG, foldContext.T)
	}

	if loop_len == 1 {
		// if len 1, include the delta G of intervening NN (SantaLucia 2004)
		pair := Pair(foldContext.Seq, start, i1, end, j1)
		if _, ok := foldContext.Energies.NN[pair]; !ok {
			return 0, fmt.Errorf("bulge: pair %q not in the NN energies", pair)
		}
		dG += Stack(start, i1, end, j1, foldContext)
	}

	// penalize AT terminal bonds
	for _, k := range []int{start, i1, end, j1} {
		if foldContext.Seq[k] == 'A' {
			dG += 0.5
		}
	}

	return dG, nil
}

func addBranch(structure NucleicAcidStructure, branches *[]Subsequence, foldContext FoldContext) error {
	if !structure.Valid() || len(structure.Inner) == 0 {
		return nil
	}
	if len(structure.Inner) == 1 {
		*branches = append(*branches, structure.Inner[0])
		return nil
	}
	for _, subsq := range structure.Inner {
		str, err := W(subsq.Start, subsq.End, foldContext)
		if err != nil {
			return err
		}
		err = addBranch(str, branches, foldContext)
		if err != nil {
			return err
		}
	}
	return nil
}

// Multibranch calculats a multi-branch energy penalty using a linear formula.
//
// From Jaeger, Turner, and Zuker, 1989.
// Found to be better than logarithmic in Ward, et al. 2017
// Args:
//
//		start: The left starting index
//		k: The mid-point in the search
//		end: The right ending index
//	 foldContext: The FoldingContext for this sequence
//		helix: Whether this Multibranch is enclosed by a helix
//		helix: Whether V(start, end) bond with one another in a helix
//
// Returns a multi-branch structure
func Multibranch(start, k, end int, foldContext FoldContext, helix bool) (NucleicAcidStructure, error) {
	var (
		left, right NucleicAcidStructure
		err         error
	)
	if helix {
		left, err = W(start+1, k, foldContext)
		if err != nil {
			return StructDefault, fmt.Errorf("Multibranch: subsequence (%d, %d, %d): %w", start, end, k, err)
		}
		right, err = W(k+1, end-1, foldContext)
		if err != nil {
			return StructDefault, fmt.Errorf("Multibranch: subsequence (%d, %d, %d): %w", start, end, k, err)
		}
	} else {
		left, err = W(start, k, foldContext)
		if err != nil {
			return StructDefault, fmt.Errorf("Multibranch: subsequence (%d, %d, %d): %w", start, end, k, err)
		}
		right, err = W(k+1, end, foldContext)
		if err != nil {
			return StructDefault, fmt.Errorf("Multibranch: subsequence (%d, %d, %d): %w", start, end, k, err)
		}
	}

	if !left.Valid() || !right.Valid() {
		return StructInvalid, nil
	}

	// gather all branches of this multi-branch structure
	var branches []Subsequence

	// in python this was a recursive closure, in Go this is not possible so
	// we pull it out and pass all the parameters
	err = addBranch(left, &branches, foldContext)
	if err != nil {
		return StructDefault, fmt.Errorf("Multibranch: subsequence (%d, %d, %d): %w", start, end, k, err)
	}
	err = addBranch(right, &branches, foldContext)
	if err != nil {
		return StructDefault, fmt.Errorf("Multibranch: subsequence (%d, %d, %d): %w", start, end, k, err)
	}

	// this isn't multi-branched
	if len(branches) < 2 {
		return StructInvalid, nil
	}

	// if there's a helix, start,end counts as well
	if helix {
		branches = append(branches, Subsequence{start, end})
	}

	// count up unpaired bp and asymmetry
	branches_count := len(branches)
	unpaired := 0
	e_sum := 0.0
	ij := Subsequence{start, end}
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
				de = Stack(i2-1, i2, j2+1, j2, foldContext)
			} else if unpaired_right != 0 {
				de = Stack(-1, i2, j2+1, j2, foldContext)
				if unpaired_right == 1 {
					de = min(Stack(i3, -1, j3, j3-1, foldContext), de)
				}
			}
		} else if ij2 == ij {
			unpaired_left = j2 - j1 - 1
			unpaired_right = i3 - i2 - 1

			if unpaired_left != 0 && unpaired_right != 0 {
				de = Stack(i2-1, i2, j2+1, j2, foldContext)
			} else if unpaired_right != 0 {
				de = Stack(i2, i2+1, j2, -1, foldContext)
				if unpaired_right == 1 {
					de = min(Stack(i3-1, i3, -1, j3, foldContext), de)
				}
			}
		} else {
			unpaired_left = i2 - j1 - 1
			unpaired_right = i3 - j2 - 1

			if unpaired_left != 0 && unpaired_right != 0 {
				de = Stack(i2-1, i2, j2+1, j2, foldContext)
			} else if unpaired_right != 0 {
				de = Stack(-1, i2, j2+1, j2, foldContext)
				if unpaired_right == 1 {
					de = min(Stack(i2-1, i2, j2+1, j2, foldContext), de)
				}
			}
		}
		e_sum += de
		unpaired += unpaired_right
		if unpaired_right < 0 {
			return StructDefault, fmt.Errorf("Multibranch: subsequence (%d, %d, %d): unpaired_right < 0", start, end, k)
		}

		if ij2 != ij { // add energy
			w, err := W(i2, j2, foldContext)
			if err != nil {
				return StructDefault, fmt.Errorf("Multibranch: subsequence (%d, %d, %d): %w", start, end, k, err)
			}
			e_sum += w.Energy
		}

	}

	if unpaired < 0 {
		return StructDefault, fmt.Errorf("Multibranch: subsequence (%d, %d, %d): unpaired < 0", start, end, k)
	}

	// penalty for unmatched bp and multi-branch
	a, b, c, d := foldContext.Energies.Multibranch.A, foldContext.Energies.Multibranch.B, foldContext.Energies.Multibranch.C, foldContext.Energies.Multibranch.D
	e_Multibranch := a + b*float64(len(branches)) + c*float64(unpaired)

	if unpaired == 0 {
		e_Multibranch = a + d
	}

	// energy of min-energy neighbors
	e := e_Multibranch + e_sum

	// pointer to next structures
	if helix {
		// branches.pop()
		branches = branches[:len(branches)-1]
	}

	return NucleicAcidStructure{Energy: e, Desc: fmt.Sprintf("BIFURCATION:%dn/%dh", unpaired, branches_count), Inner: branches}, nil
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
//		start:  The index of the start of structure on left side
//		i1: The index to the right of start
//		end:  The index of the end of structure on right side
//		j1: The index to the left of end
//	 foldContext: The FoldingContext for this sequence
//
// Returns the free energy associated with the internal loop
func InternalLoop(start, i1, end, j1 int, foldContext FoldContext) (float64, error) {
	loop_left := i1 - start - 1
	loop_right := end - j1 - 1
	loop_len := loop_left + loop_right

	if loop_left < 1 {
		return 0, fmt.Errorf("internal loop: subsequence (%d, %d, %d, %d): missing left part of the loop", start, i1, end, j1)

	}
	if loop_right < 1 {
		return 0, fmt.Errorf("internal loop: subsequence (%d, %d, %d, %d): missing right part of the loop", start, i1, end, j1)
	}

	// single bp mismatch, sum up the two single mismatch pairs
	if loop_left == 1 && loop_right == 1 {
		mm_left := Stack(start, i1, end, j1, foldContext)
		mm_right := Stack(i1-1, i1, j1+1, j1, foldContext)
		return mm_left + mm_right, nil
	}
	var d_h, d_s, dG float64
	// apply a penalty based on loop size
	if en, ok := foldContext.Energies.InternalLoops[loop_len]; ok {
		d_h, d_s = en.EnthalpyH, en.EntropyS
		dG = DeltaG(d_h, d_s, foldContext.T)
	} else {
		// it's too large an internal loop, extrapolate
		en := foldContext.Energies.InternalLoops[30]
		d_h, d_s = en.EnthalpyH, en.EntropyS
		dG = DeltaG(d_h, d_s, foldContext.T)
		dG = JacobsonStockmayer(loop_len, 30, dG, foldContext.T)
	}

	// apply an asymmetry penalty
	loop_asymmetry := math.Abs(float64(loop_left - loop_right))
	dG += 0.3 * loop_asymmetry

	// apply penalty based on the mismatching pairs on either side of the loop
	pair_left_mm := Pair(foldContext.Seq, start, start+1, end, end-1)
	en := foldContext.Energies.TERMINAL_MM[pair_left_mm]
	d_h, d_s = en.EnthalpyH, en.EntropyS
	dG += DeltaG(d_h, d_s, foldContext.T)

	pair_right_mm := Pair(foldContext.Seq, i1-1, i1, j1+1, j1)
	en = foldContext.Energies.TERMINAL_MM[pair_right_mm]
	d_h, d_s = en.EnthalpyH, en.EntropyS
	dG += DeltaG(d_h, d_s, foldContext.T)

	return dG, nil
}

// Stack return the free energy of a stack
//
// Using the indexes start and end, check whether it's at the end of
// the sequence or internal. Then check whether it's a match
// or mismatch, and return.
// Two edge-cases are terminal mismatches and dangling ends.
// The energy of a dangling end is added to the energy of a pair
// where start XOR end is at the sequence's end.
// Args:
//
//		start: The start index on left side of the pair/stack
//		i1: The index to the right of start
//		end: The end index on right side of the pair/stack
//		j1: The index to the left of end
//	 foldContext: The FoldingContext for this sequence
//
// Returns the free energy of the NN pairing
func Stack(start, i1, end, j1 int, foldContext FoldContext) float64 {
	// if any(x >= len(seq) for x in [start, i1, end, j1]):
	//    return 0.0
	for _, x := range []int{start, i1, end, j1} {
		if x >= len(foldContext.Seq) {
			return 0
		}
	}

	pair := Pair(foldContext.Seq, start, i1, end, j1)
	// if any(x == -1 for x in [start, i1, end, j1]):
	for _, x := range []int{start, i1, end, j1} {
		if x == -1 {
			// it's a dangling end
			en := foldContext.Energies.DE[pair]
			d_h, d_s := en.EnthalpyH, en.EntropyS
			return DeltaG(d_h, d_s, foldContext.T)
		}
	}

	if start > 0 && end < len(foldContext.Seq)-1 {
		// it's internal
		en, ok := foldContext.Energies.NN[pair]
		if !ok {
			en = foldContext.Energies.INTERNAL_MM[pair]
		}
		d_h, d_s := en.EnthalpyH, en.EntropyS
		return DeltaG(d_h, d_s, foldContext.T)
	}
	if start == 0 && end == len(foldContext.Seq)-1 {
		// it's terminal
		en, ok := foldContext.Energies.NN[pair]
		if !ok {
			en = foldContext.Energies.INTERNAL_MM[pair]
		}
		d_h, d_s := en.EnthalpyH, en.EntropyS
		return DeltaG(d_h, d_s, foldContext.T)
	}

	if start > 0 && end == len(foldContext.Seq)-1 {
		// it's dangling on left
		en, ok := foldContext.Energies.NN[pair]
		if !ok {
			en = foldContext.Energies.INTERNAL_MM[pair]
		}
		d_h, d_s := en.EnthalpyH, en.EntropyS
		dG := DeltaG(d_h, d_s, foldContext.T)

		pair_de := fmt.Sprintf("%c%c/.%c", foldContext.Seq[start-1], foldContext.Seq[start], foldContext.Seq[end])
		if en, ok := foldContext.Energies.DE[pair_de]; ok {
			d_h, d_s := en.EnthalpyH, en.EntropyS
			dG += DeltaG(d_h, d_s, foldContext.T)
		}
		return dG
	}

	if start == 0 && end < len(foldContext.Seq)-1 {
		// it's dangling on right
		en, ok := foldContext.Energies.NN[pair]
		if !ok {
			en = foldContext.Energies.INTERNAL_MM[pair]
		}
		d_h, d_s := en.EnthalpyH, en.EntropyS
		dG := DeltaG(d_h, d_s, foldContext.T)

		pair_de := fmt.Sprintf(".%c/%c%c", +foldContext.Seq[start], foldContext.Seq[end+1], foldContext.Seq[end])
		if en, ok := foldContext.Energies.DE[pair_de]; ok {
			d_h, d_s := en.EnthalpyH, en.EntropyS
			dG += DeltaG(d_h, d_s, foldContext.T)
			return dG
		}
	}
	return 0
}

// Hairpin calculates the free energy of a hairpin.
// Args:
//
//		start:  The index of start of hairpin
//		end:  The index of end of hairpin
//	 foldContext: The FoldingContext for this sequence
//
// Returns the free energy increment from the hairpin structure
func Hairpin(start, end int, foldContext FoldContext) (float64, error) {
	if end-start < 4 {
		return math.Inf(1), nil
	}

	hairpinSeq := foldContext.Seq[start : end+1]
	hairpin_len := len(hairpinSeq) - 2
	pair := Pair(foldContext.Seq, start, start+1, end, end-1)

	if foldContext.Energies.Complement[hairpinSeq[0]] != hairpinSeq[len(hairpinSeq)-1] {
		// not known terminal pair, nothing to close "hairpin"
		return 0, fmt.Errorf("hairpin: subsequence (%d, %d): unknown hairpin terminal pairing %c - %c", start, end, hairpinSeq[0], hairpinSeq[len(hairpinSeq)-1])
	}

	dG := 0.0
	if foldContext.Energies.TriTetraLoops != nil {
		if en, ok := foldContext.Energies.TriTetraLoops[hairpinSeq]; ok {
			// it's a pre-known hairpin with known value
			d_h, d_s := en.EnthalpyH, en.EntropyS
			dG = DeltaG(d_h, d_s, foldContext.T)
		}
	}

	// add penalty based on size
	if en, ok := foldContext.Energies.HairpinLoops[hairpin_len]; ok {
		d_h, d_s := en.EnthalpyH, en.EntropyS
		dG += DeltaG(d_h, d_s, foldContext.T)
	} else {
		// it's too large, extrapolate
		en := foldContext.Energies.HairpinLoops[30]
		d_h, d_s := en.EnthalpyH, en.EntropyS
		d_g_inc := DeltaG(d_h, d_s, foldContext.T)
		dG += JacobsonStockmayer(hairpin_len, 30, d_g_inc, foldContext.T)
	}

	// add penalty for a terminal mismatch
	en, ok := foldContext.Energies.TERMINAL_MM[pair]
	if hairpin_len > 3 && ok {
		d_h, d_s := en.EnthalpyH, en.EntropyS
		dG += DeltaG(d_h, d_s, foldContext.T)
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
//	start: leftmost index
//	i1: index to right of start
//	end: rightmost index
//	j1: index to left of end
//
// Returns string representation of the Pair
func Pair(s string, start, i1, end, j1 int) string {
	ss := []rune(s)
	ret := []rune{'.', '.', '/', '.', '.'}
	if start >= 0 {
		ret[0] = ss[start]
	}
	if i1 >= 0 {
		ret[1] = ss[i1]
	}
	if end >= 0 {
		ret[3] = ss[end]
	}
	if j1 >= 0 {
		ret[4] = ss[j1]
	}
	return string(ret)
}

// Traceback thru the V(start,end) and W(start,end) caches to find the structure
// For each step, get to the lowest energy W(start,end) within that block
// Store the structure in W(start,end)
// Inc start and end
// If the next structure is viable according to V(start,end), store as well
// Repeat
// Args:
//
//		start: The leftmost index to start searching in
//		end: The rightmost index to start searching in
//	 foldContext: The FoldingContext for this sequence
//
// Returns a list of NucleicAcidStructure in the final secondary structure
func Traceback(start, end int, foldContext FoldContext) []NucleicAcidStructure {
	// move start,end down-left to start coordinates
	s := foldContext.W[start][end]
	if !strings.Contains(s.Desc, "HAIRPIN") {
		for foldContext.W[start+1][end].Equal(s) {
			start += 1
		}
		for foldContext.W[start][end-1].Equal(s) {
			end -= 1
		}
	}

	NucleicAcidStructures := []NucleicAcidStructure{}
	for {
		s = foldContext.V[start][end]

		NucleicAcidStructures = append(NucleicAcidStructures, NucleicAcidStructure{Energy: s.Energy, Desc: s.Desc, Inner: []Subsequence{{Start: start, End: end}}})

		// it's a hairpin, end of structure
		if len(s.Inner) == 0 {
			// set the energy of everything relative to the hairpin
			return trackbackEnergy(NucleicAcidStructures)
		}

		// it's a stack, bulge, etc
		// there's another single structure beyond this
		if len(s.Inner) == 1 {
			start, end = s.Inner[0].Start, s.Inner[0].End
			// ij = s.IJs[0]
			continue
		}

		// it's a Multibranch
		e_sum := 0.0
		NucleicAcidStructures = trackbackEnergy(NucleicAcidStructures)
		branches := []NucleicAcidStructure{}
		for _, ij1 := range s.Inner {
			i1, j1 := ij1.Start, ij1.End
			tb := Traceback(i1, j1, foldContext)
			if len(tb) > 0 && len(tb[0].Inner) > 0 {
				ij2 := tb[0].Inner[0]
				i2, j2 := ij2.Start, ij2.End
				e_sum += foldContext.W[i2][j2].Energy
				branches = append(branches, tb...)
			}
		}

		NucleicAcidStructures[len(NucleicAcidStructures)-1].Energy -= e_sum
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
	minimumStructure := StructInvalid
	for _, str := range structures {
		if str.Energy != math.Inf(-1) && str.Energy < minimumStructure.Energy {
			minimumStructure = str
		}
	}
	return minimumStructure
}

// trackbackEnergy add energy to each structure, based on how it's
// W(start,end) differs from the one after
// Args:
//
//	structures: The NucleicAcidStructure for whom energy is being calculated
//
// Returns a slice of NucleicAcidStructure in the folded DNA with energy
func trackbackEnergy(structures []NucleicAcidStructure) []NucleicAcidStructure {
	for start := 0; start < len(structures)-1; start++ {
		structures[start].Energy -= structures[start+1].Energy
	}
	return structures
}
