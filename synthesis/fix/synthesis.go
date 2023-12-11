/*
Package fix fixes synthetic DNA molecules in preparation for synthesis.

Many synthesis companies have restrictions on the DNA they can synthesize. This
synthesis fixer takes advantage of synonymous codons in protein coding
sequences (CDS) to remove problematic sequences that either users don't want
(like restriction enzymes sites) or that would cause DNA synthesis companies to
reject a synthesis project.

This synthesis fixer is meant to cover the majority of use cases for DNA
fixing. It is not intended to cover all possible use cases, since the majority
of DNA design does not actually have these edge cases.

For most users, using `CdsSimple` will be sufficient to prepare a sequence
for synthesis (you may want to add in restriction enzyme sites to remove).

Cds does not guarantee that all requested features will be removed. If you
have use case that Cds cannot properly fix, please put an issue in the poly
github.
*/
package fix

import (
	"errors"
	"fmt"
	"regexp"
	"sort"
	"strings"
	"sync"

	"github.com/bebop/poly/checks"
	"github.com/bebop/poly/synthesis/codon"
	"github.com/bebop/poly/transform"
)

// DnaSuggestion is a suggestion of a fixer, generated by a
// problematicSequenceFunc. Bias must be `NA`, `GC`, or `AT`, with `NA`
// representing a neutral skew.
type DnaSuggestion struct {
	Start          int    `db:"start"`
	End            int    `db:"end"`
	Bias           string `db:"gcbias"`
	QuantityFixes  int    `db:"quantityfixes"`
	SuggestionType string `db:"suggestiontype"`
}

// Change is a change to a given DNA sequence. A list of changes is given as
// the output of Cds.
type Change struct {
	Position int    `db:"position"`
	Step     int    `db:"step"`
	From     string `db:"codonfrom"`
	To       string `db:"codonto"`
	Reason   string `db:"reason"`
}

// RemoveSequence is a generator for a problematicSequenceFuncs for specific
// sequences.
func RemoveSequence(sequencesToRemove []string, reason string) func(string, chan DnaSuggestion, *sync.WaitGroup) {
	return func(sequence string, c chan DnaSuggestion, waitgroup *sync.WaitGroup) {
		var sequencesToRemoveForReverse []string
		for _, seq := range sequencesToRemove {
			reverseComplementToRemove := transform.ReverseComplement(seq)
			// Palindromes only have to be fixed once, so add a check here
			// for palindromes.
			if reverseComplementToRemove == seq {
				sequencesToRemoveForReverse = []string{seq}
			} else {
				sequencesToRemoveForReverse = []string{seq, reverseComplementToRemove}
			}
			for _, site := range sequencesToRemoveForReverse {
				re := regexp.MustCompile(site)
				locations := re.FindAllStringIndex(sequence, -1)
				for _, location := range locations {
					codonLength := 3
					position := location[0] / codonLength
					c <- DnaSuggestion{position, (location[1] / codonLength) - 1, "NA", 1, reason}
				}
			}
		}
		waitgroup.Done()
	}
}

// RemoveRepeat is a generator to make a problematicSequenceFunc for repeats.
func RemoveRepeat(repeatLen int) func(string, chan DnaSuggestion, *sync.WaitGroup) {
	return func(sequence string, c chan DnaSuggestion, waitgroup *sync.WaitGroup) {
		// Get a kmer list
		kmers := make(map[string]bool)
		for sequencePosition := 0; sequencePosition < len(sequence)-repeatLen; sequencePosition++ {
			_, alreadyFoundForward := kmers[sequence[sequencePosition:sequencePosition+repeatLen]]
			_, alreadyFoundReverse := kmers[transform.ReverseComplement(sequence[sequencePosition:sequencePosition+repeatLen])]
			kmers[sequence[sequencePosition:sequencePosition+repeatLen]] = true
			if alreadyFoundForward || alreadyFoundReverse {
				codonLength := 3
				position := sequencePosition / codonLength
				leftover := sequencePosition % codonLength
				endPosition := (sequencePosition + repeatLen) / codonLength
				switch {
				case leftover == 0:
					c <- DnaSuggestion{position, endPosition, "NA", 1, "Repeat sequence"}
				case leftover != 0:
					c <- DnaSuggestion{position, endPosition - 1, "NA", 1, "Repeat sequence"}
				}
				sequencePosition = sequencePosition + leftover
			}
		}
		waitgroup.Done()
	}
}

// GcContentFixer is a generator to increase or decrease the overall GcContent
// of a CDS. GcContent is defined as the percentage of guanine and cytosine
// base pairs in comparison to adenine and thymine base pairs. Usually, you
// want the range to be somewhere around 50%, with a decent upperBound being
// 80% GC and a decent lowerBound being 20%.
func GcContentFixer(upperBound, lowerBound float64) func(string, chan DnaSuggestion, *sync.WaitGroup) {
	return func(sequence string, c chan DnaSuggestion, waitgroup *sync.WaitGroup) {
		gcContent := checks.GcContent(sequence)
		var numberOfChanges int
		codonLength := 3
		if gcContent > upperBound {
			numberOfChanges = int((gcContent-upperBound)*float64(len(sequence))) + 1
			c <- DnaSuggestion{0, (len(sequence) / codonLength) - 1, "AT", numberOfChanges, "GcContent too high"}
		}
		if gcContent < lowerBound {
			numberOfChanges = int((lowerBound-gcContent)*float64(len(sequence))) + 1
			c <- DnaSuggestion{0, (len(sequence) / codonLength) - 1, "GC", numberOfChanges, "GcContent too low"}
		}
		waitgroup.Done()
	}
}

// getSuggestions gets suggestions from the suggestions channel. This removes
// the need for a magic number.
func getSuggestions(suggestions chan DnaSuggestion, suggestionOutputs chan []DnaSuggestion) {
	var suggestionsList []DnaSuggestion
	for {
		suggestion, more := <-suggestions
		if more {
			suggestionsList = append(suggestionsList, suggestion)
		} else {
			suggestionOutputs <- suggestionsList
			close(suggestionOutputs)
			return
		}
	}
}

// findProblems is a helper function in FixCDS that concurrently runs each
// sequence check and returns a list of all the suggested changes.
func findProblems(sequence string, problematicSequenceFuncs []func(string, chan DnaSuggestion, *sync.WaitGroup)) []DnaSuggestion {
	// Run functions to get suggestions
	suggestions := make(chan DnaSuggestion)
	suggestionOutputs := make(chan []DnaSuggestion)
	var waitgroup sync.WaitGroup
	for _, function := range problematicSequenceFuncs {
		waitgroup.Add(1)
		go function(sequence, suggestions, &waitgroup)
	}
	go getSuggestions(suggestions, suggestionOutputs)
	waitgroup.Wait()
	close(suggestions)

	suggestionsList := <-suggestionOutputs
	return suggestionsList
}

/*
# For developers

FixCDS is the core function of the synthesis fixing package. It takes a CDS
and uses degenerate codons to fix up places with undesirable sequence.

The most important part of synthesis fixing is the functions that go into it
(problematicSequenceFuncs). These provide the range to be fixed as well as the
number of fixes that should be required in that range to fix whatever is
checked for within the function.

FixCDS first builds the following maps:

1. Builds a map of codons at each position, with the last codon in the list
   being taken to rebuild the sequence.
2. Builds a map of potential codons that an input codon can be changed to. For
   example, CAC -> CAT or ATT -> ATA,ATC . There are also codon maps for GC or
   AT nucleotide bias.
3. Builds a map of codon weights. Perhaps the organism really likes CAT codons
   but doesn't code ATA or ATC very often. Given the sequence CACATT, the "CAT"
   change will have a greater relative weight than the ATA or ATC change.

From this map, FixCDS does the following operations:

1. Concurrently runs problematicSequenceFuncs on the sequence to get change
   suggestions.
2. For each output suggestion, get a list of all potential changes to the
   sequence that fix it.
3. Sort each potential change by its codon weight.
4. Append the best change to the positionMap
5. GOTO 1
6. complete

OR, in pseudocode:

1. [x...] = problematicSequenceFuncs()
2. IF len([x...]) == 0; DONE
3. y = positionMap[x][-1]
4. [a,b,c] = potentialChanges[y]
5. [c,a,b] = sort(weights[a], weights[b], weights[c])
6. positionMap[x] = append(positionMap[x], c)
7. GOTO 1

At the end, the user should get a fixed CDS as an output, as well as a list of
changes that were done to the sequence.
*/

// Cds fixes a CDS given the CDS sequence, a codon table, a list of
// functions to solve for, and a number of iterations to attempt fixing.
// Unless you are an advanced user, you should use CdsSimple.
func Cds(sequence string, codontable codon.Table, problematicSequenceFuncs []func(string, chan DnaSuggestion, *sync.WaitGroup)) (string, []Change, error) {
	codonLength := 3
	if len(sequence)%codonLength != 0 {
		return "", []Change{}, errors.New("this sequence isn't a complete CDS, please try to use a CDS without interrupted codons")
	}

	// Setup maps
	// We have a historical map, a relative weight map, and a potential changes map.
	// The historical map gives a history of modifications for a sequence. For
	// each amino acid position of a protein, we have a list of updated codons,
	// starting with the initial codons of the protein. Changes are appended to
	// this history, and at the end the sequence is created by taking each
	// position in the history map and appending the last element in the list
	// to a single sequence string.
	historicalMap := make(map[int][]string)
	weightMap := make(map[string]float64)
	naBiasMap := make(map[string][]string)
	gcBiasMap := make(map[string][]string)
	atBiasMap := make(map[string][]string)

	// Build historical maps and full amino acid weights
	aminoAcidWeightTable := make(map[string]int)
	for _, aminoAcid := range codontable.GetWeightedAminoAcids() {
		var aminoAcidTotal int
		for _, codon := range aminoAcid.Codons {
			// Get the total weights of all the codons for a given amino acid.
			// This will be used later to get a relative weight % for each codon.
			aminoAcidTotal = aminoAcidTotal + codon.Weight

			// This third loop adds in the potential codons that a given codon can
			// switch to. If a bias is required, there are fewer potential changes.
			// Ie, the bias maps.
			codonBias := strings.Count(codon.Triplet, "G") + strings.Count(codon.Triplet, "C")
			for _, toCodon := range aminoAcid.Codons {
				if codon.Triplet != toCodon.Triplet {
					toCodonBias := strings.Count(toCodon.Triplet, "G") + strings.Count(toCodon.Triplet, "C")
					switch {
					case codonBias > toCodonBias:
						atBiasMap[codon.Triplet] = append(atBiasMap[codon.Triplet], toCodon.Triplet)
					case codonBias < toCodonBias:
						gcBiasMap[codon.Triplet] = append(gcBiasMap[codon.Triplet], toCodon.Triplet)
					}
					naBiasMap[codon.Triplet] = append(naBiasMap[codon.Triplet], toCodon.Triplet)
				}
			}
		}
		// If there is an amino acid with no encoding, error with incomplete codon table
		if aminoAcidTotal == 0 {
			return "", []Change{}, errors.New("incomplete codon table")
		}
		aminoAcidWeightTable[aminoAcid.Letter] = aminoAcidTotal
	}

	// Build weight map. The weight map gives the relative normalized weight of
	// any given codon triplet.
	for _, aminoAcid := range codontable.GetWeightedAminoAcids() {
		for _, codon := range aminoAcid.Codons {
			codonWeightRatio := float64(codon.Weight) / float64(aminoAcidWeightTable[aminoAcid.Letter])
			normalizedCodonWeight := 100 * codonWeightRatio
			weightMap[codon.Triplet] = normalizedCodonWeight
		}
	}

	// Build historical map
	position := 0
	for codonPosition := 0; codonPosition < len(sequence); codonPosition = codonPosition + codonLength {
		codon := sequence[codonPosition : codonPosition+codonLength]
		historicalMap[position] = append(historicalMap[position], codon)
		position++
	}

	// For a maximum of 100 iterations, see if we can do better. Usually sequences will be solved within 1-3 rounds,
	// so 100 just effectively acts as the max cap for iterations. Once you get to 100, you pretty much know that
	// we cannot fix the sequence.
	getSequence := func(history map[int][]string) string {
		var sequence string
		for codonPosition := 0; codonPosition < len(history); codonPosition++ {
			codonHistory := history[codonPosition]
			sequence = sequence + codonHistory[len(codonHistory)-1]
		}
		return sequence
	}
	var changes []Change
	var fixIteration int
	for {
		suggestions := findProblems(sequence, problematicSequenceFuncs)
		// If there are no suggestions, break the iteration!
		if len(suggestions) == 0 {
			// Sort changes by fixIteration and position
			sort.Slice(changes, func(i, j int) bool {
				if changes[i].Step == changes[j].Step {
					return changes[i].Position < changes[j].Position
				}
				return changes[i].Step < changes[j].Step
			})
			return sequence, changes, nil
		}
		for _, suggestion := range suggestions { // if you want to add overlaps, add suggestionIndex
			// Check for a valid bias direction
			var validBias bool
			switch suggestion.Bias {
			case "NA", "GC", "AT":
				validBias = true
			}
			if !validBias {
				return sequence, []Change{}, fmt.Errorf("Invalid bias. Expected NA, GC, or AT, got %s", suggestion.Bias)
			}

			// For each suggestion, get a list of potential changes that could fix the problem.
			var potentialChanges []Change
			for positionSelector := suggestion.Start; positionSelector <= suggestion.End && positionSelector < len(historicalMap); positionSelector++ {
				codonList := historicalMap[positionSelector]
				lastCodon := codonList[len(codonList)-1]
				unavailableCodons := make(map[string]bool)
				for _, codonSite := range historicalMap[positionSelector] {
					unavailableCodons[codonSite] = true
				}
				// We will take new potential changes from the respective bias map, given the suggestion bias.
				var biasMap map[string][]string
				switch suggestion.Bias {
				case "NA":
					biasMap = naBiasMap
				case "GC":
					biasMap = gcBiasMap
				case "AT":
					biasMap = atBiasMap
				}
				for _, potentialCodon := range biasMap[lastCodon] {
					if _, ok := unavailableCodons[potentialCodon]; !ok {
						potentialChanges = append(potentialChanges, Change{positionSelector, fixIteration, lastCodon, potentialCodon, suggestion.SuggestionType})
					}
				}
			}

			// Sort potential changes by weight
			sort.SliceStable(potentialChanges, func(i, j int) bool {
				return weightMap[potentialChanges[i].To] > weightMap[potentialChanges[j].To]
			})

			// Remove sorted changes that target the same position.
			var sortedChanges []Change
			usedPositions := make(map[int]bool)
			for _, potentialChange := range potentialChanges {
				if _, ok := usedPositions[potentialChange.Position]; !ok {
					usedPositions[potentialChange.Position] = true
					sortedChanges = append(sortedChanges, potentialChange)
				}
			}

			// Make sure we have enough sorted changes after sorting/removal
			if len(sortedChanges) < suggestion.QuantityFixes {
				return sequence, []Change{}, fmt.Errorf("Too many fixes required. Number of potential fixes: %d , number of required fixes: %d", len(potentialChanges), suggestion.QuantityFixes)
			}
			targetChanges := sortedChanges[:suggestion.QuantityFixes]

			// Update historical map, changes, and sequence
			for _, targetChange := range targetChanges {
				historicalMap[targetChange.Position] = append(historicalMap[targetChange.Position], targetChange.To)
				changes = append(changes, targetChange)
				sequence = getSequence(historicalMap)
			}
		}
		fixIteration++
	}
}

// CdsSimple is FixCds with some defaults for normal usage, including
// removing of homopolymers, removing any repeat larger than 18 base pairs, and
// fixing if a CDS's gc content is above 80% or below 20%
func CdsSimple(sequence string, codontable codon.Table, sequencesToRemove []string) (string, []Change, error) {
	var functions []func(string, chan DnaSuggestion, *sync.WaitGroup)
	// Remove homopolymers
	functions = append(functions, RemoveSequence([]string{"AAAAAAAA", "GGGGGGGG"}, "Homopolymers"))

	// Remove user defined sequences
	functions = append(functions, RemoveSequence(sequencesToRemove, "Removal requested by user"))

	// Remove repeats
	functions = append(functions, RemoveRepeat(18))

	// Ensure normal GC range
	functions = append(functions, GcContentFixer(0.80, 0.20))

	return Cds(sequence, codontable, functions)
}
