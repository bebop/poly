package poly

import (
	"errors"
	"regexp"
	"sort"
	"strings"
	"sync"
)

// CloneSequence is a simple struct that can carry a circular or linear DNA sequence.
type CloneSequence struct {
	Sequence string
	Circular bool
}

// Overhang is a struct that represents the ends of a linearized sequence where Enzymes had cut.
type Overhang struct {
	Length   int
	Position int
	Forward  bool
}

// Fragment is a struct that represents linear DNA sequences with sticky ends.
type Fragment struct {
	Sequence        string
	ForwardOverhang string
	ReverseOverhang string
}

// Enzyme is a struct that represents restriction enzymes.
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

// RestictionEnzymeCut cuts a given sequence with an enzyme represented by the enzyme's name.
func RestrictionEnzymeCut(seq CloneSequence, enzymeStr string) ([]Fragment, error) {
	enzymeMap := getBaseRestrictionEnzymes()
	if _, ok := enzymeMap[enzymeStr]; ok == false {
		return []Fragment{}, errors.New("Enzyme " + enzymeStr + " not found in enzymeMap")
	}
	enzyme, _ := enzymeMap[enzymeStr]
	return RestrictionEnzymeCutEnzymeStruct(seq, enzyme)
}

// RestrictionEnzymeCutEnzymeStruct cuts a given sequence with an enzyme represented by an Enzyme struct.
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

	// If, on a linear sequence, the last overhang's position plus EnzymeSkip is over the length of the sequence, remove that overhang.
	if (seq.Circular == false) && (overhangs[len(overhangs)-1].Position+enzyme.EnzymeSkip > len(sequence)) {
		overhangs = overhangs[len(overhangs)-1]
	}

	// Sort overhangs
	sort.SliceStable(overhangs, func(i, j int) bool {
		return overhangs[i].Position < overhangs[j].Position
	})

	// Convert Overhangs into Fragments
	var currentOverhang Overhang
	var nextOverhang Overhang
	if len(overhangs) == 1 { // Check the case of a single cut
		if seq.Circular == true {
			fragmentSeqs = append(fragmentSeqs, sequence[overhangs[0].Position:]+sequence[:overhangs[0].Position])
		} else {
			fragmentSeqs = append(fragmentSeqs, sequence[overhangs[0].Position:])
			fragmentSeqs = append(fragmentSeqs, sequence[:overhangs[0].Position])
		}
	} else {
		// The following will iterate over the overhangs list to turn them into fragments
		// There are two important variables: if the sequence is circular, and if the enzyme cutting is directional. All Type IIS enzymes
		// are directional, and in normal GoldenGate reactions these fragments would be constantly cut with enzyme as the reaction runs,
		// so are removed from the output sequences. If the enzyme is not directional, all fragments are valid.
		// If the sequence is circular, there is a chance that the nextOverhang's position will be greater than the length of the original sequence.
		// This is ok, and represents a valid cut/fragmentation of a rotation of the sequence. However, everything after will be a repeat fragment
		// of current fragments, so the iteration is terminated.
		for overhangIndex := 0; overhangIndex < len(overhangs)-1; overhangIndex++ {
			currentOverhang = overhangs[overhangIndex]
			nextOverhang = overhangs[overhangIndex+1]
			if enzyme.Directional == true {
				if (currentOverhang.Forward == true) && (nextOverhang.Forward == false) {
					fragmentSeqs = append(fragmentSeqs, sequence[currentOverhang.Position:nextOverhang.Position])
				}
				if nextOverhang.Position > len(seq.Sequence) {
					break
				}
			} else {
				fragmentSeqs = append(fragmentSeqs, sequence[currentOverhang.Position:nextOverhang.Position])
				if nextOverhang.Position > len(seq.Sequence) {
					break
				}
			}
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
	// Recurse ligate simulates all possible ligations of a series of fragments. Each possible combination begins with a "seed" that fragments from the pool can be added to.
	defer wg.Done()
	// If the seed ligates to itself, we can call it done with a successful circularization!
	if seedFragment.ForwardOverhang == seedFragment.ReverseOverhang {
		c <- seedFragment.ForwardOverhang + seedFragment.Sequence
	} else {
		for _, newFragment := range fragmentList {
			// If the seedFragment's reverse overhang is ligates to a fragment's forward overhang, we can ligate those together and seed another ligation reaction
			if seedFragment.ReverseOverhang == newFragment.ForwardOverhang {
				newSeed := Fragment{seedFragment.Sequence + seedFragment.ReverseOverhang + newFragment.Sequence, seedFragment.ForwardOverhang, newFragment.ReverseOverhang}
				wg.Add(1)
				go recurseLigate(wg, c, newSeed, fragmentList)
			}
			// This checks if we can ligate the next fragment in its reverse direction. We have to be careful though - if our seed has a palindrome, it will ligate to itself
			// like [-> <- -> <- -> ...] infinitely. We check for that case here as well.
			if (seedFragment.ReverseOverhang == ReverseComplement(newFragment.ReverseOverhang)) && (seedFragment.ReverseOverhang != ReverseComplement(seedFragment.ReverseOverhang)) { // If the second statement isn't there, program will crash on palindromes
				newSeed := Fragment{seedFragment.Sequence + seedFragment.ReverseOverhang + ReverseComplement(newFragment.Sequence), seedFragment.ForwardOverhang, ReverseComplement(newFragment.ForwardOverhang)}
				wg.Add(1)
				go recurseLigate(wg, c, newSeed, fragmentList)
			}
		}
	}
}

// Ligate simulates ligation of all possible fragment combinations into circular plasmids.
func Ligate(fragments []Fragment, maxClones int) []CloneSequence {
	var wg sync.WaitGroup
	c := make(chan string, maxClones) // A buffered channel is needed to prevent blocking.
	wg.Add(len(fragments))
	for _, fragment := range fragments {
		go recurseLigate(&wg, c, fragment, fragments)
	}
	wg.Wait()
	close(c)

	// Due to how recurseLigate works, any sequence which is used in a complete build can be used to seed
	// a valid rotation of the complete build. This sorts those sequences out using their unique Seqhashes.
	var outputConstructs []CloneSequence
	var seqhashConstruct string
	var withinConstructs bool
	for construct := range c {
		seqhashConstruct, _ = Seqhash(construct)
		withinConstructs = false
		for _, outputConstruct := range outputConstructs {
			if outputConstruct == seqhashConstruct {
				withinConstructs = true
			}
		}
		if withinConstructs == false {
			outputConstructs = append(outputConstructs, CloneSequence{construct, true})
		}
	}
	return outputConstructs
}

/******************************************************************************

Specific cloning functions begin here.

******************************************************************************/

// GoldenGate simulates a GoldenGate cloning reaction.
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
