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
// are ignored in pairedMinimumFreeEnergyV(start,end). This is based on an approach described in:
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
	foldContext, err := newFoldingContext(seq, temp)
	if err != nil {
		return nil, fmt.Errorf("error creating folding context: %w", err)
	}

	// get the minimum free energy structure out of the cache
	return traceback(0, len(seq)-1, foldContext), nil
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
	return summedEnergy, nil
}

// DotBracket returns the dot-bracket notation of a secondary nucleic acid structure
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
		for _, innerSubsequence := range structure.Inner {
			if innerSubsequence.End > maxj {
				maxj = innerSubsequence.End
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
			innerSubsequence := structure.Inner[0]
			result[innerSubsequence.Start] = '('
			result[innerSubsequence.End] = ')'
		}
	}
	return string(result)
}

// unpairedMinimumFreeEnergyW returns the minimum free energy of a subsequence
// at start and terminating at end.
//
// Figure 2B in Zuker and Stiegler, 1981
// Args:
//
//		seq: The sequence being folded
//		start: The start index
//		end: The end index (inclusive)
//	 foldContext: The context for this sequence
//
// Returns the free energy for the subsequence from start to end
func unpairedMinimumFreeEnergyW(start, end int, foldContext context) (NucleicAcidStructure, error) {
	if !foldContext.unpairedMinimumFreeEnergyW[start][end].Equal(defaultStructure) {
		return foldContext.unpairedMinimumFreeEnergyW[start][end], nil
	}

	if end-start < 4 {
		foldContext.unpairedMinimumFreeEnergyW[start][end] = invalidStructure
		return foldContext.unpairedMinimumFreeEnergyW[start][end], nil
	}

	w1, err := unpairedMinimumFreeEnergyW(start+1, end, foldContext)
	if err != nil {
		return defaultStructure, fmt.Errorf("w: subsequence (%d, %d): %w", start, end, err)
	}
	w2, err := unpairedMinimumFreeEnergyW(start, end-1, foldContext)
	if err != nil {
		return defaultStructure, fmt.Errorf("w: subsequence (%d, %d): %w", start, end, err)
	}
	w3, err := pairedMinimumFreeEnergyV(start, end, foldContext)
	if err != nil {
		return defaultStructure, fmt.Errorf("w: subsequence (%d, %d): %w", start, end, err)
	}

	w4 := invalidStructure
	for k := start + 1; k < end-1; k++ {
		w4_test, err := multibranch(start, k, end, foldContext, false)
		if err != nil {
			return defaultStructure, fmt.Errorf("w: subsequence (%d, %d): %w", start, end, err)
		}

		if w4_test.Valid() && w4_test.Energy < w4.Energy {
			w4 = w4_test
		}
	}

	wret := minimumStructure(w1, w2, w3, w4)
	foldContext.unpairedMinimumFreeEnergyW[start][end] = wret
	return wret, nil
}

// pairedMinimumFreeEnergyV returns the minimum free energy of a subsequence of paired bases
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
func pairedMinimumFreeEnergyV(start, end int, foldContext context) (NucleicAcidStructure, error) {
	if !foldContext.pairedMinimumFreeEnergyV[start][end].Equal(defaultStructure) {
		return foldContext.pairedMinimumFreeEnergyV[start][end], nil
	}

	// the ends must basepair for pairedMinimumFreeEnergyV(start,end)
	if foldContext.energies.complement[foldContext.Seq[start]] != foldContext.Seq[end] {
		foldContext.pairedMinimumFreeEnergyV[start][end] = invalidStructure
		return foldContext.pairedMinimumFreeEnergyV[start][end], nil
	}
	// if the basepair is isolated, and the seq large, penalize at 1,600 kcal/mol
	// heuristic for speeding this up
	// from https://www.ncbi.nlm.nih.gov/pubmed/10329189
	isolatedOuter := true
	if start > 0 && end < len(foldContext.Seq)-1 {
		isolatedOuter = foldContext.energies.complement[foldContext.Seq[start-1]] != foldContext.Seq[end+1]
	}
	isolatedInner := foldContext.energies.complement[foldContext.Seq[start+1]] != foldContext.Seq[end-1]

	if isolatedOuter && isolatedInner {
		foldContext.pairedMinimumFreeEnergyV[start][end] = NucleicAcidStructure{Energy: 1600}
		return foldContext.pairedMinimumFreeEnergyV[start][end], nil
	}

	paired := pair(foldContext.Seq, start, start+1, end, end-1)
	hairpin, err := hairpin(start, end, foldContext)
	if err != nil {
		return defaultStructure, fmt.Errorf("v: subsequence (%d, %d): %w", start, end, err)
	}
	e1 := NucleicAcidStructure{Energy: hairpin, Description: "HAIRPIN:" + paired}
	if end-start == 4 { // small hairpin; 4bp
		foldContext.pairedMinimumFreeEnergyV[start][end] = e1
		foldContext.unpairedMinimumFreeEnergyW[start][end] = e1
		return foldContext.pairedMinimumFreeEnergyV[start][end], nil
	}

	n := len(foldContext.Seq)
	e2 := NucleicAcidStructure{Energy: math.Inf(1)}
	for i1 := start + 1; i1 < end-4; i1++ {
		for j1 := i1 + 4; j1 < end; j1++ {
			// i1 and j1 must match
			if foldContext.energies.complement[foldContext.Seq[i1]] != foldContext.Seq[j1] {
				continue
			}

			paired := pair(foldContext.Seq, start, i1, end, j1)
			pairLeft := pair(foldContext.Seq, start, start+1, end, end-1)
			pairRight := pair(foldContext.Seq, i1-1, i1, j1+1, j1)
			_, pairLeftInner := foldContext.energies.nearestNeighbors[pairLeft]
			_, pairRightInner := foldContext.energies.nearestNeighbors[pairRight]
			pairInner := pairLeftInner || pairRightInner

			isStack := i1 == start+1 && j1 == end-1
			bulgeLeft := i1 > start+1
			bulgeRight := j1 < end-1

			var (
				e2_test      float64
				e2_test_type string
				err          error
			)
			switch {
			case isStack:
				// it's a neighboring/stacking pair in a helix
				e2_test = stack(start, i1, end, j1, foldContext)
				e2_test_type = fmt.Sprintf("STACK:%s", paired)

				if start > 0 && end == n-1 || start == 0 && end < n-1 {
					// there's a dangling end
					e2_test_type = fmt.Sprintf("STACKDanglingEnds:%s", paired)
				}
			case bulgeLeft && bulgeRight && !pairInner:
				// it's an interior loop
				il, err := internalLoop(start, i1, end, j1, foldContext)
				if err != nil {
					return defaultStructure, fmt.Errorf("v: subsequence (%d, %d): %w", start, end, err)
				}
				e2_test = il
				e2_test_type = fmt.Sprintf("INTERIOR_LOOP:%d/%d", i1-start, end-j1)

				if i1-start == 2 && end-j1 == 2 {
					loopLeftIndex := foldContext.Seq[start : i1+1]
					loopRightIndex := foldContext.Seq[j1 : end+1]
					// technically an interior loop of 1. really 1bp mismatch
					e2_test_type = fmt.Sprintf("STACK:%s/%s", loopLeftIndex, transform.Reverse(loopRightIndex))
				}
			case bulgeLeft && !bulgeRight:
				// it's a bulge on the left side
				e2_test, err = Bulge(start, i1, end, j1, foldContext)
				if err != nil {
					return defaultStructure, fmt.Errorf("v: subsequence (%d, %d): %w", start, end, err)
				}
				e2_test_type = fmt.Sprintf("BULGE:%d", i1-start)
			case !bulgeLeft && bulgeRight:
				// it's a bulge on the right side
				e2_test, err = Bulge(start, i1, end, j1, foldContext)
				if err != nil {
					return defaultStructure, fmt.Errorf("v: subsequence (%d, %d): %w", start, end, err)
				}
				e2_test_type = fmt.Sprintf("BULGE:%d", end-j1)
			default:
				// it's basically a hairpin, only outside bp match
				continue
			}

			// add pairedMinimumFreeEnergyV(start', end')
			tv, err := pairedMinimumFreeEnergyV(i1, j1, foldContext)
			if err != nil {
				return defaultStructure, fmt.Errorf("v: subsequence (%d, %d): %w", start, end, err)
			}
			e2_test += tv.Energy
			if e2_test != math.Inf(-1) && e2_test < e2.Energy {
				e2 = NucleicAcidStructure{Energy: e2_test, Description: e2_test_type, Inner: []Subsequence{{i1, j1}}}
			}
		}
	}

	e3 := invalidStructure
	if !isolatedOuter || start == 0 || end == len(foldContext.Seq)-1 {
		for k := start + 1; k < end-1; k++ {
			e3_test, err := multibranch(start, k, end, foldContext, true)
			if err != nil {
				return defaultStructure, fmt.Errorf("v: subsequence (%d, %d): %w", start, end, err)
			}

			if e3_test.Valid() && e3_test.Energy < e3.Energy {
				e3 = e3_test
			}
		}
	}
	e := minimumStructure(e1, e2, e3)
	foldContext.pairedMinimumFreeEnergyV[start][end] = e
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
func Bulge(start, i1, end, j1 int, foldContext context) (float64, error) {
	loopLength := max(i1-start-1, end-j1-1)
	if loopLength <= 0 {
		return 0, fmt.Errorf("bulge: the length of the bulge at (%d, %d) is %d", start, end, loopLength)
	}

	var dG float64

	// add penalty based on size
	if foldEnergy, ok := foldContext.energies.bulgeLoops[loopLength]; ok {
		enthalpyHDifference, entropySDifference := foldEnergy.EnthalpyH, foldEnergy.EntropyS
		dG = deltaG(enthalpyHDifference, entropySDifference, foldContext.T)
	} else {
		// it's too large for pre-calculated list, extrapolate
		foldEnergy := foldContext.energies.bulgeLoops[30]
		enthalpyHDifference, entropySDifference := foldEnergy.EnthalpyH, foldEnergy.EntropyS
		dG = deltaG(enthalpyHDifference, entropySDifference, foldContext.T)
		dG = jacobsonStockmayer(loopLength, 30, dG, foldContext.T)
	}

	if loopLength == 1 {
		// if len 1, include the delta G of intervening nearestNeighbors (SantaLucia 2004)
		paired := pair(foldContext.Seq, start, i1, end, j1)
		if _, ok := foldContext.energies.nearestNeighbors[paired]; !ok {
			return 0, fmt.Errorf("bulge: paired %q not in the nearestNeighbors energies", paired)
		}
		dG += stack(start, i1, end, j1, foldContext)
	}

	// penalize AT terminal bonds
	for _, k := range []int{start, i1, end, j1} {
		if foldContext.Seq[k] == 'A' {
			dG += 0.5
		}
	}

	return dG, nil
}

func addBranch(structure NucleicAcidStructure, branches *[]Subsequence, foldContext context) error {
	if !structure.Valid() || len(structure.Inner) == 0 {
		return nil
	}
	if len(structure.Inner) == 1 {
		*branches = append(*branches, structure.Inner[0])
		return nil
	}
	for _, inner := range structure.Inner {
		structure, err := unpairedMinimumFreeEnergyW(inner.Start, inner.End, foldContext)
		if err != nil {
			return err
		}
		err = addBranch(structure, branches, foldContext)
		if err != nil {
			return err
		}
	}
	return nil
}

// multibranch calculates a multi-branch foldEnergy penalty using a linear formula.
//
// From Jaeger, Turner, and Zuker, 1989.
// Found to be better than logarithmic in Ward, et al. 2017
// Args:
//
//		start: The left starting index
//		k: The mid-point in the search
//		end: The right ending index
//	 foldContext: The FoldingContext for this sequence
//		helix: Whether this multibranch is enclosed by a helix
//		helix: Whether pairedMinimumFreeEnergyV(start, end) bond with one another in a helix
//
// Returns a multi-branch structure
func multibranch(start, k, end int, foldContext context, helix bool) (NucleicAcidStructure, error) {
	var (
		left, right NucleicAcidStructure
		err         error
	)
	if helix {
		left, err = unpairedMinimumFreeEnergyW(start+1, k, foldContext)
		if err != nil {
			return defaultStructure, fmt.Errorf("multibranch: subsequence (%d, %d, %d): %w", start, end, k, err)
		}
		right, err = unpairedMinimumFreeEnergyW(k+1, end-1, foldContext)
		if err != nil {
			return defaultStructure, fmt.Errorf("multibranch: subsequence (%d, %d, %d): %w", start, end, k, err)
		}
	} else {
		left, err = unpairedMinimumFreeEnergyW(start, k, foldContext)
		if err != nil {
			return defaultStructure, fmt.Errorf("multibranch: subsequence (%d, %d, %d): %w", start, end, k, err)
		}
		right, err = unpairedMinimumFreeEnergyW(k+1, end, foldContext)
		if err != nil {
			return defaultStructure, fmt.Errorf("multibranch: subsequence (%d, %d, %d): %w", start, end, k, err)
		}
	}

	if !left.Valid() || !right.Valid() {
		return invalidStructure, nil
	}

	// gather all branches of this multi-branch structure
	var branches []Subsequence

	// in python this was a recursive closure, in Go this is not possible so
	// we pull it out and pass all the parameters
	err = addBranch(left, &branches, foldContext)
	if err != nil {
		return defaultStructure, fmt.Errorf("multibranch: subsequence (%d, %d, %d): %w", start, end, k, err)
	}
	err = addBranch(right, &branches, foldContext)
	if err != nil {
		return defaultStructure, fmt.Errorf("multibranch: subsequence (%d, %d, %d): %w", start, end, k, err)
	}

	// this isn't multi-branched
	if len(branches) < 2 {
		return invalidStructure, nil
	}

	// if there's a helix, start,end counts as well
	if helix {
		branches = append(branches, Subsequence{start, end})
	}

	// count up unpaired bp and asymmetry
	branchCount := len(branches)
	unpaired := 0
	summedEnergy := 0.0
	subsequence := Subsequence{start, end}
	for index, ij2 := range branches {
		i2, j2 := ij2.Start, ij2.End
		ij1 := branches[abs((index-1)%len(branches))]
		j1 := ij1.End
		ij3 := branches[abs((index+1)%len(branches))]
		i3, j3 := ij3.Start, ij3.End

		// add foldEnergy from unpaired bp to the right
		// of the helix as though it was a dangling end
		// if there's only one bp, it goes to whichever
		// helix (this or the next) has the more favorable energy
		unpairedLeft := 0
		unpairedRight := 0
		danglingEnergy := 0.0
		if index == len(branches)-1 && !helix {
			// pass
		} else if ij3 == subsequence {
			unpairedLeft = i2 - j1 - 1
			unpairedRight = j3 - j2 - 1

			if unpairedLeft != 0 && unpairedRight != 0 {
				danglingEnergy = stack(i2-1, i2, j2+1, j2, foldContext)
			} else if unpairedRight != 0 {
				danglingEnergy = stack(-1, i2, j2+1, j2, foldContext)
				if unpairedRight == 1 {
					danglingEnergy = min(stack(i3, -1, j3, j3-1, foldContext), danglingEnergy)
				}
			}
		} else if ij2 == subsequence {
			unpairedLeft = j2 - j1 - 1
			unpairedRight = i3 - i2 - 1

			if unpairedLeft != 0 && unpairedRight != 0 {
				danglingEnergy = stack(i2-1, i2, j2+1, j2, foldContext)
			} else if unpairedRight != 0 {
				danglingEnergy = stack(i2, i2+1, j2, -1, foldContext)
				if unpairedRight == 1 {
					danglingEnergy = min(stack(i3-1, i3, -1, j3, foldContext), danglingEnergy)
				}
			}
		} else {
			unpairedLeft = i2 - j1 - 1
			unpairedRight = i3 - j2 - 1

			if unpairedLeft != 0 && unpairedRight != 0 {
				danglingEnergy = stack(i2-1, i2, j2+1, j2, foldContext)
			} else if unpairedRight != 0 {
				danglingEnergy = stack(-1, i2, j2+1, j2, foldContext)
				if unpairedRight == 1 {
					danglingEnergy = min(stack(i2-1, i2, j2+1, j2, foldContext), danglingEnergy)
				}
			}
		}
		summedEnergy += danglingEnergy
		unpaired += unpairedRight
		if unpairedRight < 0 {
			return defaultStructure, fmt.Errorf("multibranch: subsequence (%d, %d, %d): unpairedRight < 0", start, end, k)
		}

		if ij2 != subsequence { // add energy
			w, err := unpairedMinimumFreeEnergyW(i2, j2, foldContext)
			if err != nil {
				return defaultStructure, fmt.Errorf("multibranch: subsequence (%d, %d, %d): %w", start, end, k, err)
			}
			summedEnergy += w.Energy
		}

	}

	if unpaired < 0 {
		return defaultStructure, fmt.Errorf("multibranch: subsequence (%d, %d, %d): unpaired < 0", start, end, k)
	}

	// penalty for unmatched bp and multi-branch
	a, b, c, d := foldContext.energies.multibranch.helicesCount, foldContext.energies.multibranch.unpairedCount, foldContext.energies.multibranch.coaxialStackCount, foldContext.energies.multibranch.terminalMismatchCount
	multibranchEnergy := a + b*float64(len(branches)) + c*float64(unpaired)

	if unpaired == 0 {
		multibranchEnergy = a + d
	}

	// energy of min-energy neighbors
	e := multibranchEnergy + summedEnergy

	// pointer to next structures
	if helix {
		// branches.pop()
		branches = branches[:len(branches)-1]
	}

	return NucleicAcidStructure{Energy: e, Description: fmt.Sprintf("BIFURCATION:%dn/%dh", unpaired, branchCount), Inner: branches}, nil
}

// internalLoop calculates the free energy of an internal loop.
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
func internalLoop(start, i1, end, j1 int, foldContext context) (float64, error) {
	loopLeftIndex := i1 - start - 1
	loopRightIndex := end - j1 - 1
	loopLength := loopLeftIndex + loopRightIndex

	if loopLeftIndex < 1 {
		return 0, fmt.Errorf("internal loop: subsequence (%d, %d, %d, %d): missing left part of the loop", start, i1, end, j1)

	}
	if loopRightIndex < 1 {
		return 0, fmt.Errorf("internal loop: subsequence (%d, %d, %d, %d): missing right part of the loop", start, i1, end, j1)
	}

	// single bp mismatch, sum up the two single mismatch pairs
	if loopLeftIndex == 1 && loopRightIndex == 1 {
		mismatchLeftEnergy := stack(start, i1, end, j1, foldContext)
		mismatchRightEnergy := stack(i1-1, i1, j1+1, j1, foldContext)
		return mismatchLeftEnergy + mismatchRightEnergy, nil
	}
	var enthalpyHDifference, entropySDifference, dG float64
	// apply a penalty based on loop size
	if foldEnergy, ok := foldContext.energies.internalLoops[loopLength]; ok {
		enthalpyHDifference, entropySDifference = foldEnergy.EnthalpyH, foldEnergy.EntropyS
		dG = deltaG(enthalpyHDifference, entropySDifference, foldContext.T)
	} else {
		// it's too large an internal loop, extrapolate
		foldEnergy := foldContext.energies.internalLoops[30]
		enthalpyHDifference, entropySDifference = foldEnergy.EnthalpyH, foldEnergy.EntropyS
		dG = deltaG(enthalpyHDifference, entropySDifference, foldContext.T)
		dG = jacobsonStockmayer(loopLength, 30, dG, foldContext.T)
	}

	// apply an asymmetry penalty
	loopAsymmetry := math.Abs(float64(loopLeftIndex - loopRightIndex))
	dG += 0.3 * loopAsymmetry

	// apply penalty based on the mismatching pairs on either side of the loop
	pairedMismatchLeftEnergy := pair(foldContext.Seq, start, start+1, end, end-1)
	foldEnergy := foldContext.energies.terminalMismatches[pairedMismatchLeftEnergy]
	enthalpyHDifference, entropySDifference = foldEnergy.EnthalpyH, foldEnergy.EntropyS
	dG += deltaG(enthalpyHDifference, entropySDifference, foldContext.T)

	pairedMismatchRightEnergy := pair(foldContext.Seq, i1-1, i1, j1+1, j1)
	foldEnergy = foldContext.energies.terminalMismatches[pairedMismatchRightEnergy]
	enthalpyHDifference, entropySDifference = foldEnergy.EnthalpyH, foldEnergy.EntropyS
	dG += deltaG(enthalpyHDifference, entropySDifference, foldContext.T)

	return dG, nil
}

// stack return the free energy of a stack
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
// Returns the free energy of the nearestNeighbors pairing
func stack(start, i1, end, j1 int, foldContext context) float64 {
	// if any(x >= len(seq) for x in [start, i1, end, j1]):
	//    return 0.0
	for _, x := range []int{start, i1, end, j1} {
		if x >= len(foldContext.Seq) {
			return 0
		}
	}

	paired := pair(foldContext.Seq, start, i1, end, j1)
	// if any(x == -1 for x in [start, i1, end, j1]):
	for _, x := range []int{start, i1, end, j1} {
		if x == -1 {
			// it's a dangling end
			foldEnergy := foldContext.energies.danglingEnds[paired]
			enthalpyHDifference, entropySDifference := foldEnergy.EnthalpyH, foldEnergy.EntropyS
			return deltaG(enthalpyHDifference, entropySDifference, foldContext.T)
		}
	}

	if start > 0 && end < len(foldContext.Seq)-1 {
		// it's internal
		foldEnergy, ok := foldContext.energies.nearestNeighbors[paired]
		if !ok {
			foldEnergy = foldContext.energies.internalMismatches[paired]
		}
		enthalpyHDifference, entropySDifference := foldEnergy.EnthalpyH, foldEnergy.EntropyS
		return deltaG(enthalpyHDifference, entropySDifference, foldContext.T)
	}
	if start == 0 && end == len(foldContext.Seq)-1 {
		// it's terminal
		foldEnergy, ok := foldContext.energies.nearestNeighbors[paired]
		if !ok {
			foldEnergy = foldContext.energies.internalMismatches[paired]
		}
		enthalpyHDifference, entropySDifference := foldEnergy.EnthalpyH, foldEnergy.EntropyS
		return deltaG(enthalpyHDifference, entropySDifference, foldContext.T)
	}

	if start > 0 && end == len(foldContext.Seq)-1 {
		// it's dangling on left
		foldEnergy, ok := foldContext.energies.nearestNeighbors[paired]
		if !ok {
			foldEnergy = foldContext.energies.internalMismatches[paired]
		}
		enthalpyHDifference, entropySDifference := foldEnergy.EnthalpyH, foldEnergy.EntropyS
		dG := deltaG(enthalpyHDifference, entropySDifference, foldContext.T)

		pairDanglingEnds := fmt.Sprintf("%c%c/.%c", foldContext.Seq[start-1], foldContext.Seq[start], foldContext.Seq[end])
		if foldEnergy, ok := foldContext.energies.danglingEnds[pairDanglingEnds]; ok {
			enthalpyHDifference, entropySDifference := foldEnergy.EnthalpyH, foldEnergy.EntropyS
			dG += deltaG(enthalpyHDifference, entropySDifference, foldContext.T)
		}
		return dG
	}

	if start == 0 && end < len(foldContext.Seq)-1 {
		// it's dangling on right
		foldEnergy, ok := foldContext.energies.nearestNeighbors[paired]
		if !ok {
			foldEnergy = foldContext.energies.internalMismatches[paired]
		}
		enthalpyHDifference, entropySDifference := foldEnergy.EnthalpyH, foldEnergy.EntropyS
		dG := deltaG(enthalpyHDifference, entropySDifference, foldContext.T)

		pairDanglingEnds := fmt.Sprintf(".%c/%c%c", +foldContext.Seq[start], foldContext.Seq[end+1], foldContext.Seq[end])
		if foldEnergy, ok := foldContext.energies.danglingEnds[pairDanglingEnds]; ok {
			enthalpyHDifference, entropySDifference := foldEnergy.EnthalpyH, foldEnergy.EntropyS
			dG += deltaG(enthalpyHDifference, entropySDifference, foldContext.T)
			return dG
		}
	}
	return 0
}

// hairpin calculates the free energy of a hairpin.
// Args:
//
//		start:  The index of start of hairpin
//		end:  The index of end of hairpin
//	 foldContext: The FoldingContext for this sequence
//
// Returns the free energy increment from the hairpin structure
func hairpin(start, end int, foldContext context) (float64, error) {
	if end-start < 4 {
		return math.Inf(1), nil
	}

	hairpinSeq := foldContext.Seq[start : end+1]
	hairpinLength := len(hairpinSeq) - 2
	paired := pair(foldContext.Seq, start, start+1, end, end-1)

	if foldContext.energies.complement[hairpinSeq[0]] != hairpinSeq[len(hairpinSeq)-1] {
		// not known terminal pair, nothing to close "hairpin"
		return 0, fmt.Errorf("hairpin: subsequence (%d, %d): unknown hairpin terminal pairing %c - %c", start, end, hairpinSeq[0], hairpinSeq[len(hairpinSeq)-1])
	}

	dG := 0.0
	if foldContext.energies.triTetraLoops != nil {
		if energy, ok := foldContext.energies.triTetraLoops[hairpinSeq]; ok {
			// it's a pre-known hairpin with known value
			enthalpyHDifference, entropySDifference := energy.EnthalpyH, energy.EntropyS
			dG = deltaG(enthalpyHDifference, entropySDifference, foldContext.T)
		}
	}

	// add penalty based on size
	if energy, ok := foldContext.energies.hairpinLoops[hairpinLength]; ok {
		enthalpyHDifference, entropySDifference := energy.EnthalpyH, energy.EntropyS
		dG += deltaG(enthalpyHDifference, entropySDifference, foldContext.T)
	} else {
		// it's too large, extrapolate
		energy := foldContext.energies.hairpinLoops[30]
		enthalpyHDifference, entropySDifference := energy.EnthalpyH, energy.EntropyS
		d_g_inc := deltaG(enthalpyHDifference, entropySDifference, foldContext.T)
		dG += jacobsonStockmayer(hairpinLength, 30, d_g_inc, foldContext.T)
	}

	// add penalty for a terminal mismatch
	energy, ok := foldContext.energies.terminalMismatches[paired]
	if hairpinLength > 3 && ok {
		enthalpyHDifference, entropySDifference := energy.EnthalpyH, energy.EntropyS
		dG += deltaG(enthalpyHDifference, entropySDifference, foldContext.T)
	}

	// add penalty if length 3 and AT closing, formula 8 from SantaLucia, 2004
	if hairpinLength == 3 && (hairpinSeq[0] == 'A' || hairpinSeq[len(hairpinSeq)-1] == 'A') {
		dG += 0.5 // convert to entropy
	}

	return dG, nil
}

// Find the free energy given delta h, s and temp
// Args:
//
//	enthalpyHDifference: The enthalpy increment in kcal / mol
//	entropySDifference: The entropy increment in cal / mol
//	temp: The temperature in Kelvin
//
// Returns the free energy increment in kcal / (mol x K)
func deltaG(enthalpyHDifference, entropySDifference, temp float64) float64 {
	return enthalpyHDifference - temp*(entropySDifference/1000.0)
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
func jacobsonStockmayer(query_len, known_len int, d_g_x, temp float64) float64 {
	gasConstant := 1.9872e-3
	return d_g_x + 2.44*gasConstant*temp*math.Log(float64(query_len)/float64(known_len))
}

// pair Returns a stack representation, a key for the nearestNeighbors maps
// Args:
//
//	s: Sequence being folded
//	start: leftmost index
//	i1: index to right of start
//	end: rightmost index
//	j1: index to left of end
//
// Returns string representation of the pair
func pair(s string, start, i1, end, j1 int) string {
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

// Traceback thru the pairedMinimumFreeEnergyV(start,end) and unpairedMinimumFreeEnergyW(start,end) caches to find the structure
// For each step, get to the lowest energy unpairedMinimumFreeEnergyW(start,end) within that block
// Store the structure in unpairedMinimumFreeEnergyW(start,end)
// Inc start and end
// If the next structure is viable according to pairedMinimumFreeEnergyV(start,end), store as well
// Repeat
// Args:
//
//		start: The leftmost index to start searching in
//		end: The rightmost index to start searching in
//	 foldContext: The FoldingContext for this sequence
//
// Returns a list of NucleicAcidStructure in the final secondary structure
func traceback(start, end int, foldContext context) []NucleicAcidStructure {
	// move start,end down-left to start coordinates
	structure := foldContext.unpairedMinimumFreeEnergyW[start][end]
	if !strings.Contains(structure.Description, "HAIRPIN") {
		for foldContext.unpairedMinimumFreeEnergyW[start+1][end].Equal(structure) {
			start += 1
		}
		for foldContext.unpairedMinimumFreeEnergyW[start][end-1].Equal(structure) {
			end -= 1
		}
	}

	NucleicAcidStructures := []NucleicAcidStructure{}
	for {
		structure = foldContext.pairedMinimumFreeEnergyV[start][end]

		NucleicAcidStructures = append(NucleicAcidStructures, NucleicAcidStructure{Energy: structure.Energy, Description: structure.Description, Inner: []Subsequence{{Start: start, End: end}}})

		// it's a hairpin, end of structure
		if len(structure.Inner) == 0 {
			// set the energy of everything relative to the hairpin
			return trackbackEnergy(NucleicAcidStructures)
		}

		// it's a stack, bulge, etc
		// there's another single structure beyond this
		if len(structure.Inner) == 1 {
			start, end = structure.Inner[0].Start, structure.Inner[0].End
			// subsequence = structure.IJs[0]
			continue
		}

		// it's a multibranch
		summedEnergy := 0.0
		NucleicAcidStructures = trackbackEnergy(NucleicAcidStructures)
		branches := []NucleicAcidStructure{}
		for _, ij1 := range structure.Inner {
			i1, j1 := ij1.Start, ij1.End
			tb := traceback(i1, j1, foldContext)
			if len(tb) > 0 && len(tb[0].Inner) > 0 {
				ij2 := tb[0].Inner[0]
				i2, j2 := ij2.Start, ij2.End
				summedEnergy += foldContext.unpairedMinimumFreeEnergyW[i2][j2].Energy
				branches = append(branches, tb...)
			}
		}

		NucleicAcidStructures[len(NucleicAcidStructures)-1].Energy -= summedEnergy
		return append(NucleicAcidStructures, branches...)
	}
}

// Return the struct with the lowest free energy that isn't -inf
// Args:
//
//	structures: NucleicAcidStructure being compared
//
// Returns the min free energy structure
func minimumStructure(structures ...NucleicAcidStructure) NucleicAcidStructure {
	minimumStructure := invalidStructure
	for _, structure := range structures {
		if structure.Energy != math.Inf(-1) && structure.Energy < minimumStructure.Energy {
			minimumStructure = structure
		}
	}
	return minimumStructure
}

// trackbackEnergy add energy to each structure, based on how it's
// unpairedMinimumFreeEnergyW(start,end) differs from the one after
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
