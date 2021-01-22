package poly

import (
	"errors"
	"regexp"
	"sort"
	"strings"
	"sync"
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

type Fragment struct {
	Sequence        string
	ForwardOverhang string
	ReverseOverhang string
}

type Enzyme struct {
	EnzymeName        string
	RegexpFor         *regexp.Regexp
	RegexpRev         *regexp.Regexp
	EnzymeSkip        int
	EnzymeOverhangLen int
	Directional       bool
}

/******************************************************************************

Base cloning functions begin here.

******************************************************************************/

// https://qvault.io/2019/10/21/golang-constant-maps-slices/
func getBaseRestrictionEnzymes() map[string]Enzyme {
	// Eventually, we want to get the data for this map from ftp://ftp.neb.com/pub/rebase
	enzymeMap := make(map[string]Enzyme)

	// Build default enzymes
	enzymeMap["BsaI"] = Enzyme{"BsaI", regexp.MustCompile("GGTCTC"), regexp.MustCompile("GAGACC"), 1, 4, true}
	enzymeMap["BbsI"] = Enzyme{"BbsI", regexp.MustCompile("GAAGAC"), regexp.MustCompile("GTCTTC"), 2, 4, true}
	enzymeMap["BtgZI"] = Enzyme{"BtgZI", regexp.MustCompile("GCGATG"), regexp.MustCompile("CATCGC"), 10, 4, true}
	enzymeMap["BsmBI"] = Enzyme{"BsmBI", regexp.MustCompile("CGTCTC"), regexp.MustCompile("GAGACG"), 1, 4, true}

	// Return EnzymeMap
	return enzymeMap
}

func RestrictionEnzymeCut(seq CloneSequence, enzymeStr string) ([]Fragment, error) {
	enzymeMap := getBaseRestrictionEnzymes()
	if _, ok := enzymeMap[enzymeStr]; ok == false {
		return []Fragment{}, errors.New("Enzyme " + enzymeStr + " not found in enzymeMap")
	}
	enzyme, _ := enzymeMap[enzymeStr]
	return RestrictionEnzymeCutEnzymeStruct(seq, enzyme), nil
}

func RestrictionEnzymeCutEnzymeStruct(seq CloneSequence, enzyme Enzyme) []Fragment {
	var fragmentSeqs []string
	var sequence string
	if seq.Circular == true {
		sequence = strings.ToUpper(seq.Sequence + seq.Sequence)
	} else {
		sequence = strings.ToUpper(seq.Sequence)
	}
	// Find and define overhangs
	var overhangs []Overhang
	forwardCuts := enzyme.RegexpFor.FindAllStringIndex(sequence, -1)
	reverseCuts := enzyme.RegexpRev.FindAllStringIndex(sequence, -1)
	for _, forwardCut := range forwardCuts {
		overhangs = append(overhangs, Overhang{Length: enzyme.EnzymeOverhangLen, Position: forwardCut[1] + enzyme.EnzymeSkip, Forward: true})
	}
	for _, reverseCut := range reverseCuts {
		overhangs = append(overhangs, Overhang{Length: enzyme.EnzymeOverhangLen, Position: reverseCut[0] - enzyme.EnzymeSkip, Forward: false})
	}
	// TODO add checks that make sure overhangs + enzymeSkip is over or below the length of the string

	// Sort overhangs
	sort.SliceStable(overhangs, func(i, j int) bool {
		return overhangs[i].Position < overhangs[j].Position
	})

	// Convert Overhangs into Fragments
	var currentOverhang Overhang
	var nextOverhang Overhang
	if len(overhangs) == 1 { // Check the case of a single cut
	}
	for i := 0; i < len(overhangs)-1; i++ {
		currentOverhang = overhangs[i]
		nextOverhang = overhangs[i+1]
		switch {
		case (seq.Circular == true) && (enzyme.Directional == true):
			if nextOverhang.Position > len(seq.Sequence) {
				if (currentOverhang.Forward == true) && (nextOverhang.Forward == false) {
					fragmentSeqs = append(fragmentSeqs, sequence[currentOverhang.Position:nextOverhang.Position])
					break
				}
			}
		case (seq.Circular == false) && (enzyme.Directional == true):
			if (currentOverhang.Forward == true) && (nextOverhang.Forward == false) {
				fragmentSeqs = append(fragmentSeqs, sequence[currentOverhang.Position:nextOverhang.Position])
			}
		case enzyme.Directional == false:
			fragmentSeqs = append(fragmentSeqs, sequence[currentOverhang.Position:nextOverhang.Position])
		}

	}

	// Convert fragment sequences into fragments
	var fragments []Fragment
	for _, fragment := range fragmentSeqs {
		fragments = append(fragments, Fragment{Sequence: fragment[enzyme.EnzymeOverhangLen : len(fragment)-enzyme.EnzymeOverhangLen], ForwardOverhang: fragment[:enzyme.EnzymeOverhangLen], ReverseOverhang: fragment[len(fragment)-enzyme.EnzymeOverhangLen:]})
	}
	return fragments
}

func recurseLigate(wg *sync.WaitGroup, c chan string, seedFragment Fragment, fragmentList []Fragment) {
	defer wg.Done()
	if seedFragment.ForwardOverhang == seedFragment.ReverseOverhang {
		c <- seedFragment.ForwardOverhang + seedFragment.Sequence
	} else {
		for _, newFragment := range fragmentList {
			if seedFragment.ReverseOverhang == newFragment.ForwardOverhang {
				newSeed := Fragment{seedFragment.Sequence + seedFragment.ReverseOverhang + newFragment.Sequence, seedFragment.ForwardOverhang, newFragment.ReverseOverhang}
				wg.Add(1)
				go recurseLigate(wg, c, newSeed, fragmentList)
			}
			if (seedFragment.ReverseOverhang == ReverseComplement(newFragment.ReverseOverhang)) && (seedFragment.ReverseOverhang != ReverseComplement(seedFragment.ReverseOverhang)) { // If the second statement isn't there, program will crash on palindromes
				newSeed := Fragment{seedFragment.Sequence + seedFragment.ReverseOverhang + ReverseComplement(newFragment.Sequence), seedFragment.ForwardOverhang, ReverseComplement(newFragment.ForwardOverhang)}
				wg.Add(1)
				go recurseLigate(wg, c, newSeed, fragmentList)
			}
		}
	}
}

func Ligate(fragments []Fragment, maxClones int) []CloneSequence {
	var wg sync.WaitGroup
	c := make(chan string, maxClones) // A buffered channel is needed to prevent blocking.
	wg.Add(len(fragments))
	for _, fragment := range fragments {
		go recurseLigate(&wg, c, fragment, fragments)
	}
	wg.Wait()
	close(c)

	// For any given cloning reaction, there will be rotations that fulfill the requirement. This removes
	// Rotations. Later, will be replaced with the seqhash algorithm
	var rotatedConstructs []CloneSequence
	var rotatedSeq string
	var withinConstructs bool
	for construct := range c {
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

/******************************************************************************

Specific cloning functions begin here.

******************************************************************************/

func GoldenGate(sequences []CloneSequence, enzymeStr string) ([]CloneSequence, error) {
	var fragments []Fragment
	for _, sequence := range sequences {
		newFragments, err := RestrictionEnzymeCut(sequence, enzymeStr)
		if err != nil {
			return []CloneSequence{}, err
		}
		fragments = append(fragments, newFragments...)
	}
	return Ligate(fragments, 1000), nil
}
