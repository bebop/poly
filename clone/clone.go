/*
Package clone provides functions for cloning DNA sequences.

Since 1973, the most common way to make recombinant DNA has been restriction
enzyme cloning (though lately, homologous recombination based methods like
Gibson assembly have attracted a lot of use). The cloning functions here allow
for simulation of restriction enzyme cloning.

For a historical review leading up to the discovery:
https://doi.org/10.1073/pnas.1313397110

The idea of restriction enzyme cloning is that you can cut DNA at specific
locations with restriction enzyme and then glue them back together in different
patterns using ligase. The final product is (99.9% of the time) a circular plasmid
that you can transform into a bacterial cell for propagation.

While simulation is simple for simple cases, there are a lot of edge cases to handle, for example:
  - Which input sequences are circular? How do we handle their rotations?
  - Is the enzyme that is cutting directional? How do we handle that directionality?
  - Are there multiple possible outputs of our ligation reaction? For example, ligations may be
    able to create a "library" of plasmids, in which there are millions of valid combinations.
  - How do we handle sequences that get ligated in multiple orientations?

These cloning functions handle all those problems so that they appear simple to the end user.

In particular, there is a focus here on GoldenGate Assembly:
https://en.wikipedia.org/wiki/Golden_Gate_Cloning
https://www.neb.com/applications/cloning-and-synthetic-biology/dna-assembly-and-cloning/golden-gate-assembly

GoldenGate is a particular kind of restriction enzyme cloning reaction that you can do
in a single tube and that is extraordinarily efficient (up to 50 parts) and is popular
for new modular DNA part toolkits. Users can easily simulate GoldenGate assembly reactions
with just their input fragments + the enzyme name.

Unlike many other GoldenGate simulators, we support simulating GoldenGate with
methylated DNA sequences, which are represented as lowercased sequences in user
inputted sequences. Normally, this can be turned off, but can be used in the
special case of recursive GoldenGate reactions.

Let's build some DNA!

# Keoni

PS: We do NOT (yet) handle restriction enzymes which recognize one site but cut
in multiple places (Type IIG enzymes) such as BcgI.
*/
package clone

import (
	"errors"
	"regexp"
	"sort"
	"strings"
	"sync"

	"github.com/TimothyStiles/poly/checks"
	"github.com/TimothyStiles/poly/seqhash"
	"github.com/TimothyStiles/poly/transform"
)

// Part is a simple struct that can carry a circular or linear DNA sequence.
// In the field of synthetic biology, the term "DNA Part" was popularized by
// the iGEM competition http://parts.igem.org/Main_Page , so we use that term
// here.
type Part struct {
	Sequence string
	Circular bool
}

// Overhang is a struct that represents the ends of a linearized sequence where Enzymes had cut.
type Overhang struct {
	Length                        int
	Position                      int
	Forward                       bool
	RecognitionSitePlusSkipLength int
}

// Fragment is a struct that represents linear DNA sequences with sticky ends.
type Fragment struct {
	Sequence        string
	ForwardOverhang string
	ReverseOverhang string
}

// Enzyme is a struct that represents restriction enzymes.
type Enzyme struct {
	Name            string
	RegexpFor       *regexp.Regexp
	RegexpRev       *regexp.Regexp
	Skip            int
	OverhangLen     int
	RecognitionSite string
}

// Eventually, we want to get the data for this map from ftp://ftp.neb.com/pub/rebase
var enzymeMap = map[string]Enzyme{
	"BsaI":  {"BsaI", regexp.MustCompile("GGTCTC"), regexp.MustCompile("GAGACC"), 1, 4, "GGTCTC"},
	"BbsI":  {"BbsI", regexp.MustCompile("GAAGAC"), regexp.MustCompile("GTCTTC"), 2, 4, "GAAGAC"},
	"BtgZI": {"BtgZI", regexp.MustCompile("GCGATG"), regexp.MustCompile("CATCGC"), 10, 4, "GCGATG"},
}

/******************************************************************************

Base cloning functions begin here.

******************************************************************************/

func getBaseRestrictionEnzymes() map[string]Enzyme {
	return enzymeMap
}

// CutWithEnzymeByName cuts a given sequence with an enzyme represented by the
// enzyme's name. It is a convenience wrapper around CutWithEnzyme that
// allows us to specify the enzyme by name. Set methylated flag to true if
// there is lowercase methylated DNA as part of the sequence.
func CutWithEnzymeByName(seq Part, directional bool, enzymeStr string, methylated bool) ([]Fragment, error) {
	enzymeMap := getBaseRestrictionEnzymes()
	if _, ok := enzymeMap[enzymeStr]; !ok {
		return []Fragment{}, errors.New("Enzyme " + enzymeStr + " not found in enzymeMap")
	}
	enzyme := enzymeMap[enzymeStr]
	return CutWithEnzyme(seq, directional, enzyme, methylated), nil
}

// CutWithEnzyme cuts a given sequence with an enzyme represented by an Enzyme struct.
// If there is methylated parts of the target DNA, set the "methylated" flag to
// true and lowercase ONLY methylated DNA.
func CutWithEnzyme(seq Part, directional bool, enzyme Enzyme, methylated bool) []Fragment {
	var fragmentSeqs []string

	// Setup circular sequences
	sequence := seq.Sequence
	if seq.Circular {
		sequence = sequence + sequence
	}

	// If unmethylated, set everything to uppercase so that the enzyme regex
	// works on all the sequence
	if !methylated {
		sequence = strings.ToUpper(sequence)
	}

	// Check for palindromes
	palindromic := checks.IsPalindromic(enzyme.RecognitionSite)

	// Find and define overhangs
	var overhangs []Overhang
	var forwardOverhangs []Overhang
	var reverseOverhangs []Overhang
	forwardCuts := enzyme.RegexpFor.FindAllStringIndex(sequence, -1)
	for _, forwardCut := range forwardCuts {
		forwardOverhangs = append(forwardOverhangs, Overhang{Length: enzyme.OverhangLen, Position: forwardCut[1] + enzyme.Skip, Forward: true, RecognitionSitePlusSkipLength: len(enzyme.RecognitionSite) + enzyme.Skip})
	}
	// Palindromic enzymes won't need reverseCuts
	if !palindromic {
		reverseCuts := enzyme.RegexpRev.FindAllStringIndex(sequence, -1)
		for _, reverseCut := range reverseCuts {
			reverseOverhangs = append(reverseOverhangs, Overhang{Length: enzyme.OverhangLen, Position: reverseCut[0] - enzyme.Skip, Forward: false, RecognitionSitePlusSkipLength: len(enzyme.RecognitionSite) + enzyme.Skip})
		}
	}

	// If, on a linear sequence, the last overhang's position + EnzymeSkip + EnzymeOverhangLen is over the length of the sequence, remove that overhang.
	for _, overhangSet := range [][]Overhang{forwardOverhangs, reverseOverhangs} {
		if len(overhangSet) > 0 {
			if !seq.Circular && (overhangSet[len(overhangSet)-1].Position+enzyme.Skip+enzyme.OverhangLen > len(sequence)) {
				overhangSet = overhangSet[:len(overhangSet)-1]
			}
		}
		overhangs = append(overhangs, overhangSet...)
	}

	// Sort overhangs
	sort.SliceStable(overhangs, func(i, j int) bool {
		return overhangs[i].Position < overhangs[j].Position
	})

	// Convert Overhangs into Fragments
	var fragments []Fragment
	var currentOverhang Overhang
	var nextOverhang Overhang
	// Linear fragments with 1 cut that are no directional will always give a
	// 2 fragments
	if len(overhangs) == 1 && !directional && !seq.Circular { // Check the case of a single cut
		// In the case of a single cut in a linear sequence, we get two fragments with only 1 stick end
		fragmentSeq1 := sequence[overhangs[0].Position+overhangs[0].Length:]
		fragmentSeq2 := sequence[:overhangs[0].Position]
		overhangSeq := sequence[overhangs[0].Position : overhangs[0].Position+overhangs[0].Length]
		fragments = append(fragments, Fragment{fragmentSeq1, overhangSeq, ""})
		fragments = append(fragments, Fragment{fragmentSeq2, "", overhangSeq})
		return fragments
	}

	// Circular fragments with 1 cut will always have 2 overhangs (because of the
	// concat earlier). If we don't require directionality, this will always get
	// cut into a single fragment
	if len(overhangs) == 2 && !directional && seq.Circular {
		// In the case of a single cut in a circular sequence, we get one fragment out with sticky overhangs
		fragmentSeq1 := sequence[overhangs[0].Position+overhangs[0].Length : len(seq.Sequence)]
		fragmentSeq2 := sequence[:overhangs[0].Position]
		fragmentSeq := fragmentSeq1 + fragmentSeq2
		overhangSeq := sequence[overhangs[0].Position : overhangs[0].Position+overhangs[0].Length]
		fragments = append(fragments, Fragment{fragmentSeq, overhangSeq, overhangSeq})
		return fragments
	}

	if len(overhangs) > 1 {
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
			// If we want directional cutting and the enzyme is not palindromic, we
			// can remove fragments that are continuously cut by the enzyme. This is
			// the basis of GoldenGate assembly.
			if directional && !palindromic {
				if currentOverhang.Forward && !nextOverhang.Forward {
					fragmentSeqs = append(fragmentSeqs, sequence[currentOverhang.Position:nextOverhang.Position])
				}
				// We have to subtract RecognitionSitePlusSkipLength in case we have a recognition site on
				// one side of the origin of a circular sequence and the cut site on the other side of the origin
				if nextOverhang.Position-nextOverhang.RecognitionSitePlusSkipLength > len(seq.Sequence) {
					break
				}
			} else {
				fragmentSeqs = append(fragmentSeqs, sequence[currentOverhang.Position:nextOverhang.Position])
				if nextOverhang.Position-nextOverhang.RecognitionSitePlusSkipLength > len(seq.Sequence) {
					break
				}
			}
		}
		// Convert fragment sequences into fragments
		for _, fragment := range fragmentSeqs {
			// Minimum lengths (given oligos) for assembly is 8 base pairs
			// https://doi.org/10.1186/1756-0500-3-291
			if len(fragment) > 8 {
				fragmentSequence := fragment[enzyme.OverhangLen : len(fragment)-enzyme.OverhangLen]
				forwardOverhang := fragment[:enzyme.OverhangLen]
				reverseOverhang := fragment[len(fragment)-enzyme.OverhangLen:]
				fragments = append(fragments, Fragment{Sequence: fragmentSequence, ForwardOverhang: forwardOverhang, ReverseOverhang: reverseOverhang})
			}
		}
	}

	return fragments
}

func recurseLigate(wg *sync.WaitGroup, constructs chan string, infiniteLoopingConstructs chan string, seedFragment Fragment, fragmentList []Fragment, usedFragments []Fragment) {
	// Recurse ligate simulates all possible ligations of a series of fragments. Each possible combination begins with a "seed" that fragments from the pool can be added to.
	defer wg.Done()
	// If the seed ligates to itself, we can call it done with a successful circularization!
	if seedFragment.ForwardOverhang == seedFragment.ReverseOverhang {
		constructs <- seedFragment.ForwardOverhang + seedFragment.Sequence
	} else {
		for _, newFragment := range fragmentList {
			// If the seedFragment's reverse overhang is ligates to a fragment's forward overhang, we can ligate those together and seed another ligation reaction
			var newSeed Fragment
			var fragmentAttached bool
			if seedFragment.ReverseOverhang == newFragment.ForwardOverhang {
				fragmentAttached = true
				newSeed = Fragment{seedFragment.Sequence + seedFragment.ReverseOverhang + newFragment.Sequence, seedFragment.ForwardOverhang, newFragment.ReverseOverhang}
			}
			// This checks if we can ligate the next fragment in its reverse direction. We have to be careful though - if our seed has a palindrome, it will ligate to itself
			// like [-> <- -> <- -> ...] infinitely. We check for that case here as well.
			if (seedFragment.ReverseOverhang == transform.ReverseComplement(newFragment.ReverseOverhang)) && (seedFragment.ReverseOverhang != transform.ReverseComplement(seedFragment.ReverseOverhang)) { // If the second statement isn't there, program will crash on palindromes
				fragmentAttached = true
				newSeed = Fragment{seedFragment.Sequence + seedFragment.ReverseOverhang + transform.ReverseComplement(newFragment.Sequence), seedFragment.ForwardOverhang, transform.ReverseComplement(newFragment.ForwardOverhang)}
			}

			// If fragment is actually attached, move to some checks
			if fragmentAttached {
				// If the newFragment's reverse complement already exists in the used fragment list, we need to cancel the recursion.
				for _, usedFragment := range usedFragments {
					if usedFragment.Sequence == newFragment.Sequence {
						infiniteLoopingConstructs <- usedFragment.ForwardOverhang + usedFragment.Sequence + usedFragment.ReverseOverhang
						return
					}
				}
				wg.Add(1)
				// If everything is clear, append fragment to usedFragments and recurse.
				usedFragments = append(usedFragments, newFragment)
				go recurseLigate(wg, constructs, infiniteLoopingConstructs, newSeed, fragmentList, usedFragments)
			}
		}
	}
}

func getConstructs(c chan string, constructSequences chan []string, circular bool) {
	var constructs []string
	var exists bool
	var existingSeqhashes []string
	for {
		construct, more := <-c
		if more {
			exists = false
			seqhashConstruct, _ := seqhash.Hash(construct, "DNA", circular, true)
			// Check if this construct is unique
			for _, existingSeqhash := range existingSeqhashes {
				if existingSeqhash == seqhashConstruct {
					exists = true
				}
			}
			if !exists {
				constructs = append(constructs, construct)
				existingSeqhashes = append(existingSeqhashes, seqhashConstruct)
			}
		} else {
			constructSequences <- constructs
			close(constructSequences)
			return
		}
	}
}

// CircularLigate simulates ligation of all possible fragment combinations into circular plasmids.
func CircularLigate(fragments []Fragment) ([]string, []string, error) {
	var wg sync.WaitGroup
	var outputConstructs []string
	var outputInfiniteLoopingConstructs []string
	constructs := make(chan string)
	infiniteLoopingConstructs := make(chan string) // sometimes we will get stuck in infinite loops. These are sequences with a recursion break
	constructSequences := make(chan []string)
	infiniteLoopingConstructSequences := make(chan []string)
	for _, fragment := range fragments {
		wg.Add(1)
		go recurseLigate(&wg, constructs, infiniteLoopingConstructs, fragment, fragments, []Fragment{})
	}
	go getConstructs(constructs, constructSequences, true)
	go getConstructs(infiniteLoopingConstructs, infiniteLoopingConstructSequences, false)
	wg.Wait()
	close(constructs)
	close(infiniteLoopingConstructs)
	outputConstructs = <-constructSequences
	outputInfiniteLoopingConstructs = <-infiniteLoopingConstructSequences
	return outputConstructs, outputInfiniteLoopingConstructs, nil
}

/******************************************************************************

Specific cloning functions begin here.

******************************************************************************/

// GoldenGate simulates a GoldenGate cloning reaction. As of right now, we only
// support BsaI, BbsI, BtgZI, and BsmBI. Set methylated flag to true if there
// is lowercase methylated DNA as part of the sequence.
func GoldenGate(sequences []Part, enzymeStr string, methylated bool) ([]string, []string, error) {
	var fragments []Fragment
	for _, sequence := range sequences {
		newFragments, err := CutWithEnzymeByName(sequence, true, enzymeStr, methylated)
		if err != nil {
			return []string{}, []string{}, err
		}
		fragments = append(fragments, newFragments...)
	}
	return CircularLigate(fragments)
}
