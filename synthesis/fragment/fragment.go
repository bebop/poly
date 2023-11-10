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
	"strings"

	"github.com/TimothyStiles/poly/checks"
	"github.com/TimothyStiles/poly/transform"
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
			nTotal = nTotal + mismatches[key{overhang, transform.ReverseComplement(overhang2)}]
		}
		if nTotal != nCorrect {
			efficiency = efficiency * (float64(nCorrect) / float64(nTotal))
		}
	}
	return efficiency
}

// NextOverhangs gets a list of possible next overhangs to use for an overhang
// list, along with their efficiencies. This can be used for more optimal
// fragmentation of sequences with potential degeneracy.
func NextOverhangs(currentOverhangs []string) ([]string, []float64) {
	currentOverhangMap := make(map[string]bool)
	for _, overhang := range currentOverhangs {
		currentOverhangMap[overhang] = true
	}

	// These 4 for loops generate all combinations of 4 base pairs
	// checking all 256 4mer combinations for palindromes. Palindromes
	// can cause problems in large combinatorial reactions, so are
	// removed here.
	bases := []rune{'A', 'T', 'G', 'C'}
	var overhangsToTest []string
	for _, base1 := range bases {
		for _, base2 := range bases {
			for _, base3 := range bases {
				for _, base4 := range bases {
					newOverhang := string([]rune{base1, base2, base3, base4})
					_, ok := currentOverhangMap[newOverhang]
					_, okReverse := currentOverhangMap[transform.ReverseComplement(newOverhang)]
					if !ok && !okReverse {
						if !checks.IsPalindromic(newOverhang) {
							overhangsToTest = append(overhangsToTest, newOverhang)
						}
					}
				}
			}
		}
	}

	var efficiencies []float64
	for _, overhang := range overhangsToTest {
		strandEffieciency := SetEfficiency(append(currentOverhangs, overhang))
		complementEfficiency := SetEfficiency(append(currentOverhangs, transform.ReverseComplement(overhang)))
		efficiencies = append(efficiencies, (strandEffieciency+complementEfficiency)/2)
	}
	return overhangsToTest, efficiencies
}

// NextOverhang gets next most efficient overhang to use for a given set of
// overhangs. This is useful for when developing a new set of standard
// overhangs. Note: NextOverhang is biased towards high AT overhangs, but this
// will not affect fidelity at all.
func NextOverhang(currentOverhangs []string) string {
	overhangsToTest, efficiencies := NextOverhangs(currentOverhangs)
	var efficiency float64
	var newOverhang string
	maxEfficiency := float64(0)
	for i, overhang := range overhangsToTest {
		efficiency = efficiencies[i]
		if efficiency > maxEfficiency {
			maxEfficiency = efficiency
			newOverhang = overhang
		}
	}
	return newOverhang
}

// optimizeOverhangIteration takes in a sequence and optimally fragments it.
func optimizeOverhangIteration(sequence string, minFragmentSize int, maxFragmentSize int, existingFragments []string, excludeOverhangs []string, includeOverhangs []string) ([]string, float64, error) {
	// If the sequence is smaller than maxFragment size, stop iteration.
	if len(sequence) < maxFragmentSize {
		existingFragments = append(existingFragments, sequence)
		return existingFragments, SetEfficiency(excludeOverhangs), nil
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
		maxFragmentSize = maxFragmentSizeBuffer // buffer is needed equations above pass.
	}

	// Get all sets of 4 between the min and max FragmentSize
	var bestOverhangEfficiency float64
	var bestOverhangPosition int
	var alreadyExists bool
	var buildAvailable bool
	for overhangOffset := 0; overhangOffset <= maxFragmentSize-minFragmentSize; overhangOffset++ {
		// We go from max -> min, so we can maximize the size of our fragments
		overhangPosition := maxFragmentSize - overhangOffset
		overhangToTest := sequence[overhangPosition-4 : overhangPosition]

		// Make sure overhang isn't already in set
		alreadyExists = false
		for _, excludeOverhang := range excludeOverhangs {
			if excludeOverhang == overhangToTest || transform.ReverseComplement(excludeOverhang) == overhangToTest {
				alreadyExists = true
			}
		}
		// Make sure overhang is in set of includeOverhangs. If includeOverhangs is
		// blank, skip this check.
		buildAvailable = false
		if len(includeOverhangs) == 0 {
			buildAvailable = true
		}
		for _, includeOverhang := range includeOverhangs {
			if includeOverhang == overhangToTest || transform.ReverseComplement(includeOverhang) == overhangToTest {
				buildAvailable = true
			}
		}
		if !alreadyExists && buildAvailable {
			// See if this overhang is a palindrome
			if !checks.IsPalindromic(overhangToTest) {
				// Get this overhang set's efficiency
				setEfficiency := SetEfficiency(append(excludeOverhangs, overhangToTest))

				// If this overhang is more efficient than any other found so far, set it as the best!
				if setEfficiency > bestOverhangEfficiency {
					bestOverhangEfficiency = setEfficiency
					bestOverhangPosition = overhangPosition
				}
			}
		}
	}
	// Set variables
	if bestOverhangPosition == 0 {
		return []string{}, float64(0), fmt.Errorf("bestOverhangPosition failed by equaling zero")
	}
	existingFragments = append(existingFragments, sequence[:bestOverhangPosition])
	excludeOverhangs = append(excludeOverhangs, sequence[bestOverhangPosition-4:bestOverhangPosition])
	sequence = sequence[bestOverhangPosition-4:]
	return optimizeOverhangIteration(sequence, minFragmentSize, maxFragmentSize, existingFragments, excludeOverhangs, includeOverhangs)
}

// Fragment fragments a sequence into fragments between the min and max size,
// choosing fragment ends for optimal assembly efficiency. Since fragments will
// be inserted into either a vector or primer binding sites, the first 4 and
// last 4 base pairs are the initial overhang set.
func Fragment(sequence string, minFragmentSize int, maxFragmentSize int, excludeOverhangs []string) ([]string, float64, error) {
	sequence = strings.ToUpper(sequence)
	return optimizeOverhangIteration(sequence, minFragmentSize, maxFragmentSize, []string{}, append([]string{sequence[:4], sequence[len(sequence)-4:]}, excludeOverhangs...), []string{})
}

// FragmentWithOverhangs fragments a sequence with only a certain overhang set.
// This is useful if you are constraining the set of possible overhangs when
// doing more advanced forms of cloning.
func FragmentWithOverhangs(sequence string, minFragmentSize int, maxFragmentSize int, excludeOverhangs []string, includeOverhangs []string) ([]string, float64, error) {
	sequence = strings.ToUpper(sequence)
	return optimizeOverhangIteration(sequence, minFragmentSize, maxFragmentSize, []string{}, append([]string{sequence[:4], sequence[len(sequence)-4:]}, excludeOverhangs...), includeOverhangs)
}
