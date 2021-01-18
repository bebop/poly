package poly

import (
	"regexp"
	"sort"
	"strings"
)

type CloneSequence struct {
	Sequence string
	Circular bool
}

type Overhang struct {
	Length   int
	Position int
	Forward  bool
}

func GoldenGateCloneBbsI(sequences []CloneSequence) []CloneSequence {
	var fragments []string
	enzymeSite := "GAAGAC"
	enzymeSkip := 2
	enzymeOverhangLen := 4
	simulateAll := false
	for _, seq := range sequences {
		// Check circularity
		var sequence string
		if seq.Circular == true {
			sequence = strings.ToUpper(seq.Sequence + seq.Sequence)
		} else {
			sequence = strings.ToUpper(seq.Sequence)
		}
		// Find and define overhangs
		var overhangs []Overhang
		forwardCuts := regexp.MustCompile(enzymeSite).FindAllStringIndex(sequence, -1)
		reverseCuts := regexp.MustCompile(ReverseComplement(enzymeSite)).FindAllStringIndex(sequence, -1)
		for _, forwardCut := range forwardCuts {
			overhangs = append(overhangs, Overhang{Length: enzymeOverhangLen, Position: forwardCut[1] + enzymeSkip, Forward: true})
		}
		for _, reverseCut := range reverseCuts {
			overhangs = append(overhangs, Overhang{Length: enzymeOverhangLen, Position: reverseCut[0] - enzymeSkip, Forward: false})
		}
		// TODO add checks that make sure overhangs + enzymeSkip is over or below the length of the string

		// Sort overhangs
		sort.SliceStable(overhangs, func(i, j int) bool {
			return overhangs[i].Position < overhangs[j].Position
		})

		// Convert Overhangs into Fragments
		var currentOverhang Overhang
		var nextOverhang Overhang
		for i := 0; i < len(overhangs)-1; i++ {
			currentOverhang = overhangs[i]
			nextOverhang = overhangs[i+1]
			if (seq.Circular == true) && (simulateAll == false) {
				if nextOverhang.Position > len(seq.Sequence) {
					if (currentOverhang.Forward == true) && (nextOverhang.Forward == false) {
						fragments = append(fragments, sequence[currentOverhang.Position:nextOverhang.Position])
						break
					}
				}
			} else {
				if (currentOverhang.Forward == true) && (nextOverhang.Forward == false) {
					fragments = append(fragments, sequence[currentOverhang.Position:nextOverhang.Position])
				}
			}

		}
	}
	// Now that we have Fragments, lets see if we can connect any together to make a circle
	var assembledConstructs []string
	var fragment string
	var appendable bool
	for _, fragmentStarter := range fragments {
		fragment = fragmentStarter
		appendable = true
		if appendable == true {
			appendable = false
			for _, newFragment := range fragments {
				if fragment[len(fragment)-4:] == newFragment[:4] {
					appendable = true
					fragment = fragment + newFragment[4:]
				}
				if fragment[len(fragment)-4:] == fragment[:4] {
					assembledConstructs = append(assembledConstructs, fragment[4:])
					appendable = false
				}
			}
		}
	}
	// For any given cloning reaction, there will be rotations that fulfill the requirement. This removes
	// Rotations. Later, will be replaced with the seqhash algorithm
	var rotatedConstructs []CloneSequence
	var rotatedSeq string
	var withinConstructs bool
	for _, construct := range assembledConstructs {
		rotatedSeq = RotateSequence(construct)
		withinConstructs = true
		for _, rotated := range rotatedConstructs {
			if rotated.Sequence == rotatedSeq {
				withinConstructs = false
			}
		}
		if withinConstructs == true {
			rotatedConstructs = append(rotatedConstructs, CloneSequence{rotatedSeq, true})
		}
	}
	return rotatedConstructs
}
