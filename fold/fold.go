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

// Zuker folds the DNA sequence and return the lowest free energy score.
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
func Zuker(seq string, temp float64) (Result, error) {
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
	if end-start == minLenForStruct { // small hairpin; 4bp
		foldContext.pairedMinimumFreeEnergyV[start][end] = e1
		foldContext.unpairedMinimumFreeEnergyW[start][end] = e1
		return foldContext.pairedMinimumFreeEnergyV[start][end], nil
	}

	n := len(foldContext.seq)
	e2 := nucleicAcidStructure{energy: math.Inf(1)}
	for rightOfStart := start + 1; rightOfStart < end-minLenForStruct; rightOfStart++ {
		for leftOfEnd := rightOfStart + minLenForStruct; leftOfEnd < end; leftOfEnd++ {
			// rightOfStart and leftOfEnd must match
			if foldContext.energies.complement(rune(foldContext.seq[rightOfStart])) != rune(foldContext.seq[leftOfEnd]) {
				continue
			}

			paired := pair(foldContext.seq, start, rightOfStart, end, leftOfEnd)
			pairLeft := pair(foldContext.seq, start, start+1, end, end-1)
			pairRight := pair(foldContext.seq, rightOfStart-1, rightOfStart, leftOfEnd+1, leftOfEnd)
			_, pairLeftInner := foldContext.energies.nearestNeighbors[pairLeft]
			_, pairRightInner := foldContext.energies.nearestNeighbors[pairRight]
			pairInner := pairLeftInner || pairRightInner

			isStack := rightOfStart == start+1 && leftOfEnd == end-1
			bulgeLeft := rightOfStart > start+1
			bulgeRight := leftOfEnd < end-1

			var (
				e2Test     float64
				e2TestType string
				err        error
			)
			switch {
			case isStack:
				// it's a neighboring/stacking pair in a helix
				e2Test = stack(start, rightOfStart, end, leftOfEnd, foldContext)
				e2TestType = fmt.Sprintf("STACK:%s", paired)

				if start > 0 && end == n-1 || start == 0 && end < n-1 {
					// there's a dangling end
					e2TestType = fmt.Sprintf("STACKDanglingEnds:%s", paired)
				}
			case bulgeLeft && bulgeRight && !pairInner:
				// it's an interior loop
				il, err := internalLoop(start, rightOfStart, end, leftOfEnd, foldContext)
				if err != nil {
					return defaultStructure, fmt.Errorf("v: subsequence (%d, %d): %w", start, end, err)
				}
				e2Test = il
				e2TestType = fmt.Sprintf("INTERIOR_LOOP:%d/%d", rightOfStart-start, end-leftOfEnd)

				if rightOfStart-start == 2 && end-leftOfEnd == 2 {
					loopLeftIndex := foldContext.seq[start : rightOfStart+1]
					loopRightIndex := foldContext.seq[leftOfEnd : end+1]
					// technically an interior loop of 1. really 1bp mismatch
					e2TestType = fmt.Sprintf("STACK:%s/%s", loopLeftIndex, transform.Reverse(loopRightIndex))
				}
			case bulgeLeft && !bulgeRight:
				// it's a bulge on the left side
				e2Test, err = Bulge(start, rightOfStart, end, leftOfEnd, foldContext)
				if err != nil {
					return defaultStructure, fmt.Errorf("v: subsequence (%d, %d): %w", start, end, err)
				}
				e2TestType = fmt.Sprintf("BULGE:%d", rightOfStart-start)
			case !bulgeLeft && bulgeRight:
				// it's a bulge on the right side
				e2Test, err = Bulge(start, rightOfStart, end, leftOfEnd, foldContext)
				if err != nil {
					return defaultStructure, fmt.Errorf("v: subsequence (%d, %d): %w", start, end, err)
				}
				e2TestType = fmt.Sprintf("BULGE:%d", end-leftOfEnd)
			default:
				// it's basically a hairpin, only outside bp match
				continue
			}

			// add pairedMinimumFreeEnergyV(start', end')
			tv, err := pairedMinimumFreeEnergyV(rightOfStart, leftOfEnd, foldContext)
			if err != nil {
				return defaultStructure, fmt.Errorf("v: subsequence (%d, %d): %w", start, end, err)
			}
			e2Test += tv.energy
			if e2Test != math.Inf(-1) && e2Test < e2.energy {
				e2 = nucleicAcidStructure{energy: e2Test, description: e2TestType, inner: []subsequence{{rightOfStart, leftOfEnd}}}
			}
		}
	}

	e3 := invalidStructure
	if !isolatedOuter || start == 0 || end == len(foldContext.seq)-1 {
		for k := start + 1; k < end-1; k++ {
			e3Test, err := multibranch(start, k, end, foldContext, true)
			if err != nil {
				return defaultStructure, fmt.Errorf("v: subsequence (%d, %d): %w", start, end, err)
			}

			if e3Test.Valid() && e3Test.energy < e3.energy {
				e3 = e3Test
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
//		rightOfStart: The index to the right of start
//		end: The end index of the bulge
//		leftOfEnd: The index to the left of end
//	 foldContext: The FoldingContext for this sequence
//
// Returns the increment in free energy from the bulge
func Bulge(start, rightOfStart, end, leftOfEnd int, foldContext context) (float64, error) {
	loopLength := max(rightOfStart-start-1, end-leftOfEnd-1)
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
		paired := pair(foldContext.seq, start, rightOfStart, end, leftOfEnd)
		if _, ok := foldContext.energies.nearestNeighbors[paired]; !ok {
			return 0, fmt.Errorf("bulge: paired %q not in the nearestNeighbors energies", paired)
		}
		dG += stack(start, rightOfStart, end, leftOfEnd, foldContext)
	}

	// penalize AT terminal bonds
	for _, k := range []int{start, rightOfStart, end, leftOfEnd} {
		if foldContext.seq[k] == 'A' {
			dG += closingATPenalty
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
//		mid: The mid-point in the search
//		end: The right ending index
//	 foldContext: The FoldingContext for this sequence
//		helix: Whether this multibranch is enclosed by a helix
//		helix: Whether pairedMinimumFreeEnergyV(start, end) bond with one another in a helix
//
// Returns a multi-branch structure
func multibranch(start, mid, end int, foldContext context, helix bool) (nucleicAcidStructure, error) {
	var (
		left, right nucleicAcidStructure
		err         error
	)
	if helix {
		left, err = unpairedMinimumFreeEnergyW(start+1, mid, foldContext)
		if err != nil {
			return defaultStructure, fmt.Errorf("multibranch: subsequence (%d, %d, %d): %w", start, end, mid, err)
		}
		right, err = unpairedMinimumFreeEnergyW(mid+1, end-1, foldContext)
		if err != nil {
			return defaultStructure, fmt.Errorf("multibranch: subsequence (%d, %d, %d): %w", start, end, mid, err)
		}
	} else {
		left, err = unpairedMinimumFreeEnergyW(start, mid, foldContext)
		if err != nil {
			return defaultStructure, fmt.Errorf("multibranch: subsequence (%d, %d, %d): %w", start, end, mid, err)
		}
		right, err = unpairedMinimumFreeEnergyW(mid+1, end, foldContext)
		if err != nil {
			return defaultStructure, fmt.Errorf("multibranch: subsequence (%d, %d, %d): %w", start, end, mid, err)
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
		return defaultStructure, fmt.Errorf("multibranch: subsequence (%d, %d, %d): %w", start, mid, end, err)
	}
	err = addBranch(right, &branches, foldContext)
	if err != nil {
		return defaultStructure, fmt.Errorf("multibranch: subsequence (%d, %d, %d): %w", start, mid, end, err)
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
	curSequence := subsequence{start, end}
	for index, currentBranch := range branches {
		leftStart, leftEnd := currentBranch.start, currentBranch.end
		curBranchLeft := branches[abs((index-1)%len(branches))]
		leftOfEnd := curBranchLeft.end
		curBranchRight := branches[abs((index+1)%len(branches))]
		rightStart, rightEnd := curBranchRight.start, curBranchRight.end

		// add foldEnergy from unpaired bp to the right
		// of the helix as though it was a dangling end
		// if there's only one bp, it goes to whichever
		// helix (this or the next) has the more favorable energy
		unpairedLeft := 0
		unpairedRight := 0
		danglingEnergy := 0.0
		if index == len(branches)-1 && !helix {
			// pass
		} else if curBranchRight == curSequence {
			unpairedLeft = leftStart - leftOfEnd - 1
			unpairedRight = rightEnd - leftEnd - 1

			if unpairedLeft != 0 && unpairedRight != 0 {
				danglingEnergy = stack(leftStart-1, leftStart, leftEnd+1, leftEnd, foldContext)
			} else if unpairedRight != 0 {
				danglingEnergy = stack(-1, leftStart, leftEnd+1, leftEnd, foldContext)
				if unpairedRight == 1 {
					danglingEnergy = math.Min(stack(rightStart, -1, rightEnd, rightEnd-1, foldContext), danglingEnergy)
				}
			}
		} else if currentBranch == curSequence {
			unpairedLeft = leftEnd - leftOfEnd - 1
			unpairedRight = rightStart - leftStart - 1

			if unpairedLeft != 0 && unpairedRight != 0 {
				danglingEnergy = stack(leftStart-1, leftStart, leftEnd+1, leftEnd, foldContext)
			} else if unpairedRight != 0 {
				danglingEnergy = stack(leftStart, leftStart+1, leftEnd, -1, foldContext)
				if unpairedRight == 1 {
					danglingEnergy = math.Min(stack(rightStart-1, rightStart, -1, rightEnd, foldContext), danglingEnergy)
				}
			}
		} else {
			unpairedLeft = leftStart - leftOfEnd - 1
			unpairedRight = rightStart - leftEnd - 1

			if unpairedLeft != 0 && unpairedRight != 0 {
				danglingEnergy = stack(leftStart-1, leftStart, leftEnd+1, leftEnd, foldContext)
			} else if unpairedRight != 0 {
				danglingEnergy = stack(-1, leftStart, leftEnd+1, leftEnd, foldContext)
				if unpairedRight == 1 {
					danglingEnergy = math.Min(stack(leftStart-1, leftStart, leftEnd+1, leftEnd, foldContext), danglingEnergy)
				}
			}
		}
		summedEnergy += danglingEnergy
		unpaired += unpairedRight
		if unpairedRight < 0 {
			return defaultStructure, fmt.Errorf("multibranch: subsequence (%d, %d, %d): unpairedRight < 0", start, end, mid)
		}

		if currentBranch != curSequence { // add energy
			w, err := unpairedMinimumFreeEnergyW(leftStart, leftEnd, foldContext)
			if err != nil {
				return defaultStructure, fmt.Errorf("multibranch: subsequence (%d, %d, %d): %w", start, end, mid, err)
			}
			summedEnergy += w.energy
		}
	}

	if unpaired < 0 {
		return defaultStructure, fmt.Errorf("multibranch: subsequence (%d, %d, %d): unpaired < 0", start, end, mid)
	}

	// this is just for readability of the formulas below
	var (
		helicesCount          = foldContext.energies.multibranch.helicesCount
		unpairedCount         = foldContext.energies.multibranch.unpairedCount
		coaxialStackCount     = foldContext.energies.multibranch.coaxialStackCount
		terminalMismatchCount = foldContext.energies.multibranch.terminalMismatchCount
	)

	// penalty for unmatched bp and multi-branch
	multibranchEnergy := helicesCount + unpairedCount*float64(len(branches)) + coaxialStackCount*float64(unpaired)

	if unpaired == 0 {
		multibranchEnergy = helicesCount + terminalMismatchCount
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
//	rightOfStart: The index to the right of start
//		end:  The index of the end of structure on right side
//		leftOfEnd: The index to the left of end
//	 foldContext: The FoldingContext for this sequence
//
// Returns the free energy associated with the internal loop
func internalLoop(start, rightOfStart, end, leftOfEnd int, foldContext context) (float64, error) {
	loopLeftIndex := rightOfStart - start - 1
	loopRightIndex := end - leftOfEnd - 1
	loopLength := loopLeftIndex + loopRightIndex

	if loopLeftIndex < 1 {
		return 0, fmt.Errorf("internal loop: subsequence (%d, %d, %d, %d): missing left part of the loop", start, rightOfStart, end, leftOfEnd)
	}
	if loopRightIndex < 1 {
		return 0, fmt.Errorf("internal loop: subsequence (%d, %d, %d, %d): missing right part of the loop", start, rightOfStart, end, leftOfEnd)
	}

	// single bp mismatch, sum up the two single mismatch pairs
	if loopLeftIndex == 1 && loopRightIndex == 1 {
		mismatchLeftEnergy := stack(start, rightOfStart, end, leftOfEnd, foldContext)
		mismatchRightEnergy := stack(rightOfStart-1, rightOfStart, leftOfEnd+1, leftOfEnd, foldContext)
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
	dG += loopsAsymmetryPenalty * loopAsymmetry

	// apply penalty based on the mismatching pairs on either side of the loop
	pairedMismatchLeftEnergy := pair(foldContext.seq, start, start+1, end, end-1)
	foldEnergy := foldContext.energies.terminalMismatches[pairedMismatchLeftEnergy]
	enthalpyHDifference, entropySDifference = foldEnergy.enthalpyH, foldEnergy.entropyS
	dG += deltaG(enthalpyHDifference, entropySDifference, foldContext.temp)

	pairedMismatchRightEnergy := pair(foldContext.seq, rightOfStart-1, rightOfStart, leftOfEnd+1, leftOfEnd)
	foldEnergy = foldContext.energies.terminalMismatches[pairedMismatchRightEnergy]
	enthalpyHDifference, entropySDifference = foldEnergy.enthalpyH, foldEnergy.entropyS
	dG += deltaG(enthalpyHDifference, entropySDifference, foldContext.temp)

	return dG, nil
}

// stack returns the free energy of a stack.
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
//	rightOfStart: The index to the right of start
//		end: The end index on right side of the pair/stack
//		leftOfEnd: The index to the left of end
//	 foldContext: The FoldingContext for this sequence
//
// Returns the free energy of the nearestNeighbors pairing
func stack(start, rightOfStart, end, leftOfEnd int, foldContext context) float64 {
	// if any(x >= len(seq) for x in [start,rightOfStart, end, leftOfEnd]):
	//    return 0.0
	for _, indices := range []int{start, rightOfStart, end, leftOfEnd} {
		if indices >= len(foldContext.seq) {
			return 0
		}
	}

	paired := pair(foldContext.seq, start, rightOfStart, end, leftOfEnd)
	// if any(x == -1 for x in [start,rightOfStart, end, leftOfEnd]):
	for _, indices := range []int{start, rightOfStart, end, leftOfEnd} {
		if indices == -1 {
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
		dGinc := deltaG(enthalpyHDifference, entropySDifference, foldContext.temp)
		dG += jacobsonStockmayer(hairpinLength, maxLenPreCalulated, dGinc, foldContext.temp)
	}

	// add penalty for a terminal mismatch
	energy, ok := foldContext.energies.terminalMismatches[paired]
	if hairpinLength > 3 && ok {
		enthalpyHDifference, entropySDifference := energy.enthalpyH, energy.entropyS
		dG += deltaG(enthalpyHDifference, entropySDifference, foldContext.temp)
	}

	// add penalty if length 3 and AT closing, formula 8 from SantaLucia, 2004
	if hairpinLength == 3 && (hairpinSeq[0] == 'A' || hairpinSeq[len(hairpinSeq)-1] == 'A') {
		dG += closingATPenalty // convert to entropy
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

// jacobsonStockmayer entropy extrapolation formula is used for bulges,
// hairpins, etc that fall outside the maxLenPreCalulated upper limit for
// pre-calculated free-energies. See SantaLucia and Hicks (2004)
//
// Args:
//
//	queryLen: Length of element without known free energy value
//	knownLen: Length of element with known free energy value (dGx)
//	dGx: The free energy of the element knownLen
//	temp: Temperature in Kelvin
//
// Returns the free energy for a structure of length queryLen
// See SantaLucia and Hicks (2004).
// NOTE: the coefficient 2.44 is based on recent kinetics
// measurements in DNA (Goddard NL, Bonnet G, Krichevsky O, Libchaber A. 2000),
// and thus it is used in preference to the older theoretically derived value
// of 1.75.
func jacobsonStockmayer(queryLen, knownLen int, dGx, temp float64) float64 {
	gasConstant := 1.9872e-3
	return dGx + 2.44*gasConstant*temp*math.Log(float64(queryLen)/float64(knownLen))
}

// pair Returns a stack representation, a key for the nearestNeighbors maps
// Args:
//
//	s: Sequence being folded
//	start: leftmost index
//	rightOfStart: index to right of start
//	end: rightmost index
//	leftOfEnd: index to left of end
//
// Returns string representation of the pair
func pair(s string, start, rightOfStart, end, leftOfEnd int) string {
	seqRunes := []rune(s)
	ret := []rune{'.', '.', '/', '.', '.'}
	if start >= 0 {
		ret[0] = seqRunes[start]
	}
	if rightOfStart >= 0 {
		ret[1] = seqRunes[rightOfStart]
	}
	if end >= 0 {
		ret[3] = seqRunes[end]
	}
	if leftOfEnd >= 0 {
		ret[4] = seqRunes[leftOfEnd]
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
		for _, subseq := range structure.inner {
			rightOfStart, leftOfEnd := subseq.start, subseq.end
			tb := traceback(rightOfStart, leftOfEnd, foldContext)
			if len(tb) > 0 && len(tb[0].inner) > 0 {
				subseq := tb[0].inner[0]
				subStart, subEnd := subseq.start, subseq.end
				summedEnergy += foldContext.unpairedMinimumFreeEnergyW[subStart][subEnd].energy
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
