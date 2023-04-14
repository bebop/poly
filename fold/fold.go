/*
Package fold is a package for folding DNA and RNA sequences into secondary structures.

This package provides everything you need to fold a DNA or RNA sequence into a secondary structure
and get the minimum free energy of the structure. Most of the code was ported from the
python SeqFold package by Lattice Automation and Joshua Timmons but we hope to have a
linear fold algorithm in the near future.

Biological context:

DNA, RNA, and proteins all fold. Protein is a particularly tricky thing to predict
partially because there are so many more amino acids than there are nucleotides.

ACG(T/U) vs. ACDEFGHIKLMNPQRSTVWY (20 amino acids)

These folding predictions help us understand how to design primers, guide RNAs, and
other nucleic acid sequences that fold into a particular structure.

Fortunately for us, DNA and RNA are much easier to predict because there are only 4 nucleotides
and the rules for folding are much more well defined.

Each function has citations to the original papers that describe the algorithms used.
Most of the algorithms used in this package are based on the work of Zuker and Stiegler, 1981
but we're hoping to add more algorithms in the near future such as linear fold.

TTFN,
Tim
*/
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
func Fold(seq string, temp float64) (Result, error) {
	foldContext, err := newFoldingContext(seq, temp)
	if err != nil {
		return Result{}, fmt.Errorf("error creating folding context: %w", err)
	}

	// get the minimum free energy structure out of the cache

	return Result{
		structs: traceback(0, len(seq)-1, foldContext),
	}, nil
}

// unpairedMinimumFreeEnergyW returns the minimum free energy of a subsequence
// at start and terminating at end.
//
// From Zuker and Stiegler, 1981: let W(i,j) be the minimum free energy of all
// possible admissible structures formed from the subsequence Sij.
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
func unpairedMinimumFreeEnergyW(start, end int, foldContext context) (nucleicAcidStructure, error) {
	if !foldContext.unpairedMinimumFreeEnergyW[start][end].Equal(defaultStructure) {
		return foldContext.unpairedMinimumFreeEnergyW[start][end], nil
	}

	// pentanucleotide sequences form no stable structure
	if end-start < minLenForStruct {
		foldContext.unpairedMinimumFreeEnergyW[start][end] = invalidStructure
		return foldContext.unpairedMinimumFreeEnergyW[start][end], nil
	}

	endDanglingLeft, err := unpairedMinimumFreeEnergyW(start+1, end, foldContext)
	if err != nil {
		return defaultStructure, fmt.Errorf("w: subsequence (%d, %d): %w", start, end, err)
	}
	endDanglingRight, err := unpairedMinimumFreeEnergyW(start, end-1, foldContext)
	if err != nil {
		return defaultStructure, fmt.Errorf("w: subsequence (%d, %d): %w", start, end, err)
	}
	endsPaired, err := pairedMinimumFreeEnergyV(start, end, foldContext)
	if err != nil {
		return defaultStructure, fmt.Errorf("w: subsequence (%d, %d): %w", start, end, err)
	}

	endBifurcation := invalidStructure
	for k := start + 1; k < end-1; k++ {
		testBranch, err := multibranch(start, k, end, foldContext, false)
		if err != nil {
			return defaultStructure, fmt.Errorf("w: subsequence (%d, %d): %w", start, end, err)
		}

		if testBranch.Valid() && testBranch.energy < endBifurcation.energy {
			endBifurcation = testBranch
		}
	}

	minStuctEnergy := minimumStructure(endDanglingLeft, endDanglingRight, endsPaired, endBifurcation)
	foldContext.unpairedMinimumFreeEnergyW[start][end] = minStuctEnergy
	return minStuctEnergy, nil
}

// pairedMinimumFreeEnergyV returns the minimum free energy of a subsequence of
// paired bases and end.
// From Figure 2B of Zuker, 1981: let V(i,j) be the minimum free energy of all
// possible admissible structures formed from Sij in which Si and Sj base pair
// with each other. If Si and Sj cannot base pair, then V(i,j) = infinity
//
// If start and end don't bp, store and return INF.
// See: Figure 2B of Zuker, 1981
// Args:
//
//		start: The start index
//		end: The end index (inclusive)
//	 foldContext: The context for this sequence
//
// Returns the minimum energy folding structure possible between start and end on seq
func pairedMinimumFreeEnergyV(start, end int, foldContext context) (nucleicAcidStructure, error) {
	if !foldContext.pairedMinimumFreeEnergyV[start][end].Equal(defaultStructure) {
		return foldContext.pairedMinimumFreeEnergyV[start][end], nil
	}

	// the ends must basepair for pairedMinimumFreeEnergyV(start,end)
	if foldContext.energies.complement(rune(foldContext.seq[start])) != rune(foldContext.seq[end]) {
		foldContext.pairedMinimumFreeEnergyV[start][end] = invalidStructure
		return foldContext.pairedMinimumFreeEnergyV[start][end], nil
	}
	// if the basepair is isolated, and the seq large, penalize at 1,600 kcal/mol
	// heuristic for speeding this up
	// from https://www.ncbi.nlm.nih.gov/pubmed/10329189
	isolatedOuter := true
	if start > 0 && end < len(foldContext.seq)-1 {
		isolatedOuter = foldContext.energies.complement(rune(foldContext.seq[start-1])) != rune(foldContext.seq[end+1])
	}
	isolatedInner := foldContext.energies.complement(rune(foldContext.seq[start+1])) != rune(foldContext.seq[end-1])

	if isolatedOuter && isolatedInner {
		foldContext.pairedMinimumFreeEnergyV[start][end] = nucleicAcidStructure{energy: isolatedBasePairPenalty}
		return foldContext.pairedMinimumFreeEnergyV[start][end], nil
	}

	paired := pair(foldContext.seq, start, start+1, end, end-1)
	hairpin, err := hairpin(start, end, foldContext)
	if err != nil {
		return defaultStructure, fmt.Errorf("v: subsequence (%d, %d): %w", start, end, err)
	}
	e1 := nucleicAcidStructure{energy: hairpin, description: "HAIRPIN:" + paired}
	if end-start == 4 { // small hairpin; 4bp
		foldContext.pairedMinimumFreeEnergyV[start][end] = e1
		foldContext.unpairedMinimumFreeEnergyW[start][end] = e1
		return foldContext.pairedMinimumFreeEnergyV[start][end], nil
	}

	n := len(foldContext.seq)
	e2 := nucleicAcidStructure{energy: math.Inf(1)}
	for i1 := start + 1; i1 < end-4; i1++ {
		for j1 := i1 + 4; j1 < end; j1++ {
			// i1 and j1 must match
			if foldContext.energies.complement(rune(foldContext.seq[i1])) != rune(foldContext.seq[j1]) {
				continue
			}

			paired := pair(foldContext.seq, start, i1, end, j1)
			pairLeft := pair(foldContext.seq, start, start+1, end, end-1)
			pairRight := pair(foldContext.seq, i1-1, i1, j1+1, j1)
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
					loopLeftIndex := foldContext.seq[start : i1+1]
					loopRightIndex := foldContext.seq[j1 : end+1]
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
			e2_test += tv.energy
			if e2_test != math.Inf(-1) && e2_test < e2.energy {
				e2 = nucleicAcidStructure{energy: e2_test, description: e2_test_type, inner: []subsequence{{i1, j1}}}
			}
		}
	}

	e3 := invalidStructure
	if !isolatedOuter || start == 0 || end == len(foldContext.seq)-1 {
		for k := start + 1; k < end-1; k++ {
			e3_test, err := multibranch(start, k, end, foldContext, true)
			if err != nil {
				return defaultStructure, fmt.Errorf("v: subsequence (%d, %d): %w", start, end, err)
			}

			if e3_test.Valid() && e3_test.energy < e3.energy {
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
		enthalpyHDifference, entropySDifference := foldEnergy.enthalpyH, foldEnergy.entropyS
		dG = deltaG(enthalpyHDifference, entropySDifference, foldContext.temp)
	} else {
		// it's too large for pre-calculated list, extrapolate
		foldEnergy := foldContext.energies.bulgeLoops[maxLenPreCalulated]
		enthalpyHDifference, entropySDifference := foldEnergy.enthalpyH, foldEnergy.entropyS
		dG = deltaG(enthalpyHDifference, entropySDifference, foldContext.temp)
		dG = jacobsonStockmayer(loopLength, maxLenPreCalulated, dG, foldContext.temp)
	}

	if loopLength == 1 {
		// if len 1, include the delta G of intervening nearestNeighbors (SantaLucia 2004)
		paired := pair(foldContext.seq, start, i1, end, j1)
		if _, ok := foldContext.energies.nearestNeighbors[paired]; !ok {
			return 0, fmt.Errorf("bulge: paired %q not in the nearestNeighbors energies", paired)
		}
		dG += stack(start, i1, end, j1, foldContext)
	}

	// penalize AT terminal bonds
	for _, k := range []int{start, i1, end, j1} {
		if foldContext.seq[k] == 'A' {
			dG += 0.5
		}
	}

	return dG, nil
}

func addBranch(structure nucleicAcidStructure, branches *[]subsequence, foldContext context) error {
	if !structure.Valid() || len(structure.inner) == 0 {
		return nil
	}
	if len(structure.inner) == 1 {
		*branches = append(*branches, structure.inner[0])
		return nil
	}
	for _, inner := range structure.inner {
		structure, err := unpairedMinimumFreeEnergyW(inner.start, inner.end, foldContext)
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
func multibranch(start, k, end int, foldContext context, helix bool) (nucleicAcidStructure, error) {
	var (
		left, right nucleicAcidStructure
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
	var branches []subsequence

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
		branches = append(branches, subsequence{start, end})
	}

	// count up unpaired bp and asymmetry
	branchCount := len(branches)
	unpaired := 0
	summedEnergy := 0.0
	subsequence := subsequence{start, end}
	for index, ij2 := range branches {
		i2, j2 := ij2.start, ij2.end
		ij1 := branches[abs((index-1)%len(branches))]
		j1 := ij1.end
		ij3 := branches[abs((index+1)%len(branches))]
		i3, j3 := ij3.start, ij3.end

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
					danglingEnergy = math.Min(stack(i3, -1, j3, j3-1, foldContext), danglingEnergy)
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
					danglingEnergy = math.Min(stack(i3-1, i3, -1, j3, foldContext), danglingEnergy)
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
					danglingEnergy = math.Min(stack(i2-1, i2, j2+1, j2, foldContext), danglingEnergy)
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
			summedEnergy += w.energy
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

	return nucleicAcidStructure{energy: e, description: fmt.Sprintf("BIFURCATION:%dn/%dh", unpaired, branchCount), inner: branches}, nil
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
		enthalpyHDifference, entropySDifference = foldEnergy.enthalpyH, foldEnergy.entropyS
		dG = deltaG(enthalpyHDifference, entropySDifference, foldContext.temp)
	} else {
		// it's too large an internal loop, extrapolate
		foldEnergy := foldContext.energies.internalLoops[maxLenPreCalulated]
		enthalpyHDifference, entropySDifference = foldEnergy.enthalpyH, foldEnergy.entropyS
		dG = deltaG(enthalpyHDifference, entropySDifference, foldContext.temp)
		dG = jacobsonStockmayer(loopLength, maxLenPreCalulated, dG, foldContext.temp)
	}

	// apply an asymmetry penalty
	loopAsymmetry := math.Abs(float64(loopLeftIndex - loopRightIndex))
	dG += 0.3 * loopAsymmetry

	// apply penalty based on the mismatching pairs on either side of the loop
	pairedMismatchLeftEnergy := pair(foldContext.seq, start, start+1, end, end-1)
	foldEnergy := foldContext.energies.terminalMismatches[pairedMismatchLeftEnergy]
	enthalpyHDifference, entropySDifference = foldEnergy.enthalpyH, foldEnergy.entropyS
	dG += deltaG(enthalpyHDifference, entropySDifference, foldContext.temp)

	pairedMismatchRightEnergy := pair(foldContext.seq, i1-1, i1, j1+1, j1)
	foldEnergy = foldContext.energies.terminalMismatches[pairedMismatchRightEnergy]
	enthalpyHDifference, entropySDifference = foldEnergy.enthalpyH, foldEnergy.entropyS
	dG += deltaG(enthalpyHDifference, entropySDifference, foldContext.temp)

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
		if x >= len(foldContext.seq) {
			return 0
		}
	}

	paired := pair(foldContext.seq, start, i1, end, j1)
	// if any(x == -1 for x in [start, i1, end, j1]):
	for _, x := range []int{start, i1, end, j1} {
		if x == -1 {
			// it's a dangling end
			foldEnergy := foldContext.energies.danglingEnds[paired]
			enthalpyHDifference, entropySDifference := foldEnergy.enthalpyH, foldEnergy.entropyS
			return deltaG(enthalpyHDifference, entropySDifference, foldContext.temp)
		}
	}

	if start > 0 && end < len(foldContext.seq)-1 {
		// it's internal
		foldEnergy, ok := foldContext.energies.nearestNeighbors[paired]
		if !ok {
			foldEnergy = foldContext.energies.internalMismatches[paired]
		}
		enthalpyHDifference, entropySDifference := foldEnergy.enthalpyH, foldEnergy.entropyS
		return deltaG(enthalpyHDifference, entropySDifference, foldContext.temp)
	}
	if start == 0 && end == len(foldContext.seq)-1 {
		// it's terminal
		foldEnergy, ok := foldContext.energies.nearestNeighbors[paired]
		if !ok {
			foldEnergy = foldContext.energies.internalMismatches[paired]
		}
		enthalpyHDifference, entropySDifference := foldEnergy.enthalpyH, foldEnergy.entropyS
		return deltaG(enthalpyHDifference, entropySDifference, foldContext.temp)
	}

	if start > 0 && end == len(foldContext.seq)-1 {
		// it's dangling on left
		foldEnergy, ok := foldContext.energies.nearestNeighbors[paired]
		if !ok {
			foldEnergy = foldContext.energies.internalMismatches[paired]
		}
		enthalpyHDifference, entropySDifference := foldEnergy.enthalpyH, foldEnergy.entropyS
		dG := deltaG(enthalpyHDifference, entropySDifference, foldContext.temp)

		pairDanglingEnds := fmt.Sprintf("%c%c/.%c", foldContext.seq[start-1], foldContext.seq[start], foldContext.seq[end])
		if foldEnergy, ok := foldContext.energies.danglingEnds[pairDanglingEnds]; ok {
			enthalpyHDifference, entropySDifference := foldEnergy.enthalpyH, foldEnergy.entropyS
			dG += deltaG(enthalpyHDifference, entropySDifference, foldContext.temp)
		}
		return dG
	}

	if start == 0 && end < len(foldContext.seq)-1 {
		// it's dangling on right
		foldEnergy, ok := foldContext.energies.nearestNeighbors[paired]
		if !ok {
			foldEnergy = foldContext.energies.internalMismatches[paired]
		}
		enthalpyHDifference, entropySDifference := foldEnergy.enthalpyH, foldEnergy.entropyS
		dG := deltaG(enthalpyHDifference, entropySDifference, foldContext.temp)

		pairDanglingEnds := fmt.Sprintf(".%c/%c%c", +foldContext.seq[start], foldContext.seq[end+1], foldContext.seq[end])
		if foldEnergy, ok := foldContext.energies.danglingEnds[pairDanglingEnds]; ok {
			enthalpyHDifference, entropySDifference := foldEnergy.enthalpyH, foldEnergy.entropyS
			dG += deltaG(enthalpyHDifference, entropySDifference, foldContext.temp)
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
	if end-start < minLenForStruct {
		return math.Inf(1), nil
	}

	hairpinSeq := foldContext.seq[start : end+1]
	hairpinLength := len(hairpinSeq) - 2
	paired := pair(foldContext.seq, start, start+1, end, end-1)

	if foldContext.energies.complement(rune(hairpinSeq[0])) != rune(hairpinSeq[len(hairpinSeq)-1]) {
		// not known terminal pair, nothing to close "hairpin"
		return 0, fmt.Errorf("hairpin: subsequence (%d, %d): unknown hairpin terminal pairing %c - %c", start, end, hairpinSeq[0], hairpinSeq[len(hairpinSeq)-1])
	}

	dG := 0.0
	if foldContext.energies.triTetraLoops != nil {
		if energy, ok := foldContext.energies.triTetraLoops[hairpinSeq]; ok {
			// it's a pre-known hairpin with known value
			enthalpyHDifference, entropySDifference := energy.enthalpyH, energy.entropyS
			dG = deltaG(enthalpyHDifference, entropySDifference, foldContext.temp)
		}
	}

	// add penalty based on size
	if energy, ok := foldContext.energies.hairpinLoops[hairpinLength]; ok {
		enthalpyHDifference, entropySDifference := energy.enthalpyH, energy.entropyS
		dG += deltaG(enthalpyHDifference, entropySDifference, foldContext.temp)
	} else {
		// it's too large, extrapolate
		energy := foldContext.energies.hairpinLoops[maxLenPreCalulated]
		enthalpyHDifference, entropySDifference := energy.enthalpyH, energy.entropyS
		d_g_inc := deltaG(enthalpyHDifference, entropySDifference, foldContext.temp)
		dG += jacobsonStockmayer(hairpinLength, maxLenPreCalulated, d_g_inc, foldContext.temp)
	}

	// add penalty for a terminal mismatch
	energy, ok := foldContext.energies.terminalMismatches[paired]
	if hairpinLength > 3 && ok {
		enthalpyHDifference, entropySDifference := energy.enthalpyH, energy.entropyS
		dG += deltaG(enthalpyHDifference, entropySDifference, foldContext.temp)
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
//	queryLen: Length of element without known free energy value
//	knownLen: Length of element with known free energy value (d_g_x)
//	dGx: The free energy of the element known_len
//	temp: Temperature in Kelvin
//
// Returns the free energy for a structure of length query_len
func jacobsonStockmayer(queryLen, knownLen int, dGx, temp float64) float64 {
	gasConstant := 1.9872e-3
	return dGx + 2.44*gasConstant*temp*math.Log(float64(queryLen)/float64(knownLen))
}

// pair Returns a stack representation, a key for the nearestNeighbors maps
// Args:
//
//	s: Sequence being folded
//	start: leftmost index
//	startRigh: index to right of start
//	end: rightmost index
//	endLeft: index to left of end
//
// Returns string representation of the pair
func pair(s string, start, startRigh, end, endLeft int) string {
	ss := []rune(s)
	ret := []rune{'.', '.', '/', '.', '.'}
	if start >= 0 {
		ret[0] = ss[start]
	}
	if startRigh >= 0 {
		ret[1] = ss[startRigh]
	}
	if end >= 0 {
		ret[3] = ss[end]
	}
	if endLeft >= 0 {
		ret[4] = ss[endLeft]
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
func traceback(start, end int, foldContext context) []nucleicAcidStructure {
	// move start,end down-left to start coordinates
	structure := foldContext.unpairedMinimumFreeEnergyW[start][end]
	if !strings.Contains(structure.description, "HAIRPIN") {
		for foldContext.unpairedMinimumFreeEnergyW[start+1][end].Equal(structure) {
			start += 1
		}
		for foldContext.unpairedMinimumFreeEnergyW[start][end-1].Equal(structure) {
			end -= 1
		}
	}

	NucleicAcidStructures := []nucleicAcidStructure{}
	for {
		structure = foldContext.pairedMinimumFreeEnergyV[start][end]

		NucleicAcidStructures = append(NucleicAcidStructures, nucleicAcidStructure{energy: structure.energy, description: structure.description, inner: []subsequence{{start: start, end: end}}})

		// it's a hairpin, end of structure
		if len(structure.inner) == 0 {
			// set the energy of everything relative to the hairpin
			return trackbackEnergy(NucleicAcidStructures)
		}

		// it's a stack, bulge, etc
		// there's another single structure beyond this
		if len(structure.inner) == 1 {
			start, end = structure.inner[0].start, structure.inner[0].end
			// subsequence = structure.IJs[0]
			continue
		}

		// it's a multibranch
		summedEnergy := 0.0
		NucleicAcidStructures = trackbackEnergy(NucleicAcidStructures)
		branches := []nucleicAcidStructure{}
		for _, ij1 := range structure.inner {
			i1, j1 := ij1.start, ij1.end
			tb := traceback(i1, j1, foldContext)
			if len(tb) > 0 && len(tb[0].inner) > 0 {
				ij2 := tb[0].inner[0]
				i2, j2 := ij2.start, ij2.end
				summedEnergy += foldContext.unpairedMinimumFreeEnergyW[i2][j2].energy
				branches = append(branches, tb...)
			}
		}

		NucleicAcidStructures[len(NucleicAcidStructures)-1].energy -= summedEnergy
		return append(NucleicAcidStructures, branches...)
	}
}

// Return the struct with the lowest free energy that isn't -inf
// Args:
//
//	structures: NucleicAcidStructure being compared
//
// Returns the min free energy structure
func minimumStructure(structures ...nucleicAcidStructure) nucleicAcidStructure {
	minimumStructure := invalidStructure
	for _, structure := range structures {
		if structure.energy != math.Inf(-1) && structure.energy < minimumStructure.energy {
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
func trackbackEnergy(structures []nucleicAcidStructure) []nucleicAcidStructure {
	for start := 0; start < len(structures)-1; start++ {
		structures[start].energy -= structures[start+1].energy
	}
	return structures
}
