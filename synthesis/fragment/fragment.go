/* Package fragment optimally fragments DNA for GoldenGate systems.

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
)

// SetEfficiency gives the predicted fidelity of a list of GoldenGate
// overhangs.
func SetEfficiency(x []string) float64 {
	var efficiency = float64(1.0)
	for _, overhang := range x {
		nCorrect := mismatches[key{overhang, overhang}]
		nTotal := 0
		for _, overhang2 := range x {
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

	// Make sure maxFragmentSize - minFragmentSize < 4
	if (maxFragmentSize - minFragmentSize) < 4 {
		return []string{}, float64(0), fmt.Errorf("maxFragmentSize - minFragmentSize < 4. Got size of %d", maxFragmentSize-minFragmentSize)
	}

	var bestOverhangEfficiency float64
	var bestOverhangPosition int
	var alreadyExists bool
	// If len(sequence < 2*maxFragmentSize), 
	if len(sequence) < 2*maxFragmentSize {
		maxFragmentSizeBuffer := (len(sequence)/2) + (maxFragmentSize - minFragmentSize)
		minFragmentSize = (len(sequence)/2) - (maxFragmentSize - minFragmentSize)
		maxFragmentSize = maxFragmentSizeBuffer // buffer is needed so second equation goes through.
	}
	// Get all sets of 4 between the min and max FragmentSize
	for i := 0; i <= (maxFragmentSize-minFragmentSize)-4; i++ {
		// We go from max -> min, so we can maximize the size of our fragments
		overhangPosition := maxFragmentSize - i
		overhangToTest := sequence[overhangPosition-4 : overhangPosition]

		// Make sure overhang isn't already in set
		alreadyExists = false
		for _, existingOverhang := range existingOverhangs {
			if existingOverhang == overhangToTest {
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
// choosing fragment ends for optimal assembly efficiency. Fragment is designed
// for synthesizing fragments into vectors, the most common time you would need
// optimal fragmentation, so the first 4 base pairs and last 4 base pairs are
// are a part of the fidelity equation.
func Fragment(sequence string, minFragmentSize int, maxFragmentSize int) ([]string, float64, error) {
	sequence = strings.ToUpper(sequence)
	return optimizeOverhangIteration(sequence, minFragmentSize, maxFragmentSize, []string{}, []string{sequence[:4], sequence[len(sequence)-4:]})
}
