/*
Package fragment optimally fragments DNA for GoldenGate systems.

Optimal fragmentation is accomplished by using empirical fidelity data derived
by NEB in the paper "Enabling one-pot Golden Gate assemblies of unprecedented
complexity using data-optimized assembly design". We use the BsaI-T4 ligase
data provided in table S1.

Paper link: https://doi.org/10.1371/journal.pone.0238592
Data link: https://doi.org/10.1371/journal.pone.0238592.s001
*/
package fragment

import (
	"fmt"
	"github.com/TimothyStiles/poly/transform"
	"strings"
)

// SetEfficiency gets the estimated fidelity rate of a given set of
// GoldenGate overhangs.
func SetEfficiency(overhangs []string) float64 {
	var efficiency = float64(1.0)
	for _, overhang := range overhangs {
		nCorrect := mismatches[key{overhang, overhang}]
		nTotal := 0
		for _, overhang2 := range overhangs {
			nTotal = nTotal + mismatches[key{overhang, overhang2}]
		}
		if nTotal > 0 {
			efficiency = efficiency * (float64(nCorrect) / float64(nTotal))
		}
	}
	return efficiency
}

// optimizeOverhangIteration takes in a sequence and optimally fragments it.
func optimizeOverhangIteration(sequence string, minFragmentSize int, maxFragmentSize int, existingFragments []string, existingOverhangs []string) ([]string, float64, error) {
	// If the sequence is smaller than maxFragment size, stop iteration.
	if len(sequence) < maxFragmentSize {
		existingFragments = append(existingFragments, sequence)
		return existingFragments, SetEfficiency(existingOverhangs), nil
	}

	// Make sure minFragmentSize > maxFragmentSize
	if minFragmentSize > maxFragmentSize {
		return []string{}, float64(0), fmt.Errorf("minFragmentSize (%d) larger than maxFragmentSize (%d)", minFragmentSize, maxFragmentSize)
	}

	// Minimum lengths (given oligos) for assembly is 8 base pairs
	// https://doi.org/10.1186/1756-0500-3-291
	// For GoldenGate, 2 8bp oligos create 12 base pairs (4bp overhangs on two sides of 4bp),
	// so we check for minimal size of 12 base pairs.
	if minFragmentSize < 12 {
		return []string{}, float64(0), fmt.Errorf("minFragmentSize must be equal to or greater than 12 . Got size of %d", minFragmentSize)
	}

	// If our iteration is approaching the end of the sequence, that means we need to gracefully handle
	// the end so we aren't left with a tiny fragment that cannot be synthesized. For example, if our goal
	// is fragments of 100bp, and we have 110 base pairs left, we want each final fragment to be 55bp, not
	// 100 and 10bp
	if len(sequence) < 2*maxFragmentSize {
		maxAndMinDifference := maxFragmentSize - minFragmentSize
		maxFragmentSizeBuffer := (len(sequence) + maxAndMinDifference) / 2
		if maxFragmentSizeBuffer > maxFragmentSize {
			maxFragmentSizeBuffer = maxFragmentSize
		}
		minFragmentSize = maxFragmentSizeBuffer - maxAndMinDifference
		// Basic check for minimal size
		if minFragmentSize < 12 {
			minFragmentSize = 12
		}
		maxFragmentSize = maxFragmentSizeBuffer // buffer is needed equations above pass.
	}

	// Get all sets of 4 between the min and max FragmentSize
	var bestOverhangEfficiency float64
	var bestOverhangPosition int
	var alreadyExists bool
	for overhangOffset := 0; overhangOffset <= maxFragmentSize-minFragmentSize; overhangOffset++ {
		// We go from max -> min, so we can maximize the size of our fragments
		overhangPosition := maxFragmentSize - overhangOffset
		overhangToTest := sequence[overhangPosition-4 : overhangPosition]

		// Make sure overhang isn't already in set
		alreadyExists = false
		for _, existingOverhang := range existingOverhangs {
			if existingOverhang == overhangToTest || transform.ReverseComplement(existingOverhang) == overhangToTest {
				alreadyExists = true
			}
		}
		if !alreadyExists {
			// Get this overhang set's efficiency
			setEfficiency := SetEfficiency(append(existingOverhangs, overhangToTest))

			// If this overhang is more efficient than any other found so far, set it as the best!
			if setEfficiency > bestOverhangEfficiency {
				bestOverhangEfficiency = setEfficiency
				bestOverhangPosition = overhangPosition
			}
		}
	}
	// Set variables
	if bestOverhangPosition == 0 {
		return []string{}, float64(0), fmt.Errorf("bestOverhangPosition failed by equaling zero")
	}
	existingFragments = append(existingFragments, sequence[:bestOverhangPosition])
	existingOverhangs = append(existingOverhangs, sequence[bestOverhangPosition-4:bestOverhangPosition])
	sequence = sequence[bestOverhangPosition-4:]
	return optimizeOverhangIteration(sequence, minFragmentSize, maxFragmentSize, existingFragments, existingOverhangs)
}

// Fragment fragments a sequence into fragments between the min and max size,
// choosing fragment ends for optimal assembly efficiency. Since fragments will
// be insterted into either a vector or primer binding sites, the first 4 and
// last 4 base pairs are the initial overhang set.
func Fragment(sequence string, minFragmentSize int, maxFragmentSize int) ([]string, float64, error) {
	sequence = strings.ToUpper(sequence)
	return optimizeOverhangIteration(sequence, minFragmentSize, maxFragmentSize, []string{}, []string{sequence[:4], sequence[len(sequence)-4:]})
}
