/*
Package pcr designs and simulates simple PCR reactions.

PCR, or polymerase chain reaction, is a method developed in 1983 to copy DNA
templates using small fragments of synthesized single-stranded DNA, amplifying
those DNA templates to ~x1,000,000,000 their starting concentration. These small
fragments, referred to as "primers" or "oligos", can be designed on a computer
and then synthesized for amplifying a variety of different templates.

This package allows users to simulate a PCR reaction or design new primers to
amplify a given template. This package assumes perfect annealing to template at
a target temperature, so should only be used for PCR reactions where this is a
reasonable assumption.

If you are trying to simulate amplification out of a large pool, such as an
oligo pool, use the `Simulate` rather than `SimulateSimple` function to detect
if there is concatemerization happening in your multiplex reaction. In most
other cases, use `SimulateSimple`.

IMPORTANT! The targetTm in all functions is specifically for Taq polymerase.
*/
package pcr

import (
	"errors"
	"index/suffixarray"
	"sort"
	"strings"

	"github.com/TimothyStiles/poly/primers"
	"github.com/TimothyStiles/poly/transform"
)

// https://doi.org/10.1089/dna.1994.13.75
var minimalPrimerLength int = 15

// DesignPrimersWithOverhangs designs two primers to amplify a target sequence,
// adding on an overhang to the forward and reverse strand. This overhang can
// contain additional DNA needed for assembly, like Gibson assembly overhangs
// or GoldenGate restriction enzyme sites.
func DesignPrimersWithOverhangs(sequence, forwardOverhang, reverseOverhang string, targetTm float64) (string, string) {
	sequence = strings.ToUpper(sequence)
	forwardPrimer := sequence[0:minimalPrimerLength]
	for additionalNucleotides := 0; primers.MeltingTemp(forwardPrimer) < targetTm; additionalNucleotides++ {
		forwardPrimer = sequence[0 : minimalPrimerLength+additionalNucleotides]
	}
	reversePrimer := transform.ReverseComplement(sequence[len(sequence)-minimalPrimerLength:])
	for additionalNucleotides := 0; primers.MeltingTemp(reversePrimer) < targetTm; additionalNucleotides++ {
		reversePrimer = transform.ReverseComplement(sequence[len(sequence)-(minimalPrimerLength+additionalNucleotides):])
	}

	// Add overhangs to primer
	forwardPrimer = forwardOverhang + forwardPrimer
	reversePrimer = transform.ReverseComplement(reverseOverhang) + reversePrimer

	return forwardPrimer, reversePrimer
}

// DesignPrimers designs two primers to amplify a target sequence and only that
// target sequence (no overhangs).
func DesignPrimers(sequence string, targetTm float64) (string, string) {
	return DesignPrimersWithOverhangs(sequence, "", "", targetTm)
}

// SimulateSimple simulates a PCR reaction. It takes in a list of sequences and
// a list of primers, with support for complex multiplex reactions, produces
// a list of all possible PCR fragments from such a reaction. It does not
// detect concatemerization, which could be useful or very detrimental to
// your reactions. The variable `circular` is for if the target template is
// circular, like a plasmid.
func SimulateSimple(sequences []string, targetTm float64, circular bool, primerList []string) []string {
	// Set all primers to uppercase.
	for primerIndex := range primerList {
		primerList[primerIndex] = strings.ToUpper(primerList[primerIndex])
	}

	var pcrFragments []string
	for _, sequence := range sequences {
		sequence = strings.ToUpper(sequence)
		// Suffix array construction allows function to operate on
		// very large sequences without being worried about exeuction
		// time. For small sequences, it doesn't really matter.
		// https://eli.thegreenplace.net/2016/suffix-arrays-in-the-go-standard-library/
		sequenceIndex := suffixarray.New([]byte(sequence))

		primerLength := len(primerList)

		forwardLocations := make(map[int][]int)
		reverseLocations := make(map[int][]int)
		minimalPrimers := make([]string, primerLength)
		for primerIndex, primer := range primerList {
			var minimalLength int
			for index := minimalPrimerLength; primers.MeltingTemp(primer[len(primer)-index:]) < targetTm; index++ {
				minimalLength = index
				if primer[len(primer)-index:] == primer {
					break
				}
			}
			// Use the minimal binding sites of the primer to find positions in the template
			minimalPrimer := primer[len(primer)-minimalLength:]
			if minimalPrimer != primer {
				minimalPrimers[primerIndex] = minimalPrimer
				// For each primer, we want to look for all possible binding sites in our gene.
				// We then append this to a list of binding sites for that primer.
				for _, forwardLocation := range sequenceIndex.Lookup([]byte(minimalPrimer), -1) {
					forwardLocations[forwardLocation] = append(forwardLocations[forwardLocation], primerIndex)
				}
				for _, reverseLocation := range sequenceIndex.Lookup([]byte(transform.ReverseComplement(minimalPrimer)), -1) {
					reverseLocations[reverseLocation] = append(reverseLocations[reverseLocation], primerIndex)
				}
			}
		}

		var forwardLocationInts []int
		var reverseLocationInts []int
		for key := range forwardLocations {
			forwardLocationInts = append(forwardLocationInts, key)
		}
		for key := range reverseLocations {
			reverseLocationInts = append(reverseLocationInts, key)
		}
		sort.Ints(forwardLocationInts)
		sort.Ints(reverseLocationInts)

		// Next, iterate through the forwardLocations list
		for index, forwardLocation := range forwardLocationInts {
			// First, make sure that this isn't the last element in forwardLocations
			if index+1 != len(forwardLocationInts) {
				// If this isn't the last element in forwardLocations, then we can select the first reverseLocation that is less than the next forwardLocation
				for _, reverseLocation := range reverseLocationInts {
					if (forwardLocation < reverseLocation) && (reverseLocation < forwardLocationInts[index+1]) {
						// If both are true, we have found the sequence we are aiming to PCR! Now, we get all primers from that forwardLocation and then
						// build PCR fragments with each one.
						pcrFragments = append(pcrFragments, generatePcrFragments(sequence, forwardLocation, reverseLocation, forwardLocations[forwardLocation], reverseLocations[reverseLocation], minimalPrimers, primerList)...)
						break
					}
				}
			} else {
				foundFragment := false
				for _, reverseLocation := range reverseLocationInts {
					if forwardLocation < reverseLocation {
						pcrFragments = append(pcrFragments, generatePcrFragments(sequence, forwardLocation, reverseLocation, forwardLocations[forwardLocation], reverseLocations[reverseLocation], minimalPrimers, primerList)...)
						foundFragment = true
					}
				}
				// If the sequence is circular and we haven't found a fragment yet, check the other side of the origin
				if circular && !foundFragment {
					for _, reverseLocation := range reverseLocationInts {
						if forwardLocationInts[0] > reverseLocation {
							// If either one of these are true, create a new pcrFragment and append to pcrFragments
							rotatedSequence := sequence[forwardLocation:] + sequence[:forwardLocation]
							rotatedForwardLocation := 0
							rotatedReverseLocation := len(sequence[forwardLocation:]) + reverseLocation
							pcrFragments = append(pcrFragments, generatePcrFragments(rotatedSequence, rotatedForwardLocation, rotatedReverseLocation, forwardLocations[forwardLocation], reverseLocations[reverseLocation], minimalPrimers, primerList)...)
						}
					}
				}
			}
		}
	}
	return pcrFragments
}

// Simulate simulates a PCR reaction, including concatemerization analysis. It
// takes in a list of sequences and list of primers, produces all possible PCR
// fragments in a given reaction, and then attempts to see if the output
// fragments can amplify themselves. If they can, concatemerization is occurring
// in your reaction, which can lead to confusing results. The variable
// `circular` is for if the target template is circular, like a plasmid.
func Simulate(sequences []string, targetTm float64, circular bool, primerList []string) ([]string, error) {
	initialAmplification := SimulateSimple(sequences, targetTm, circular, primerList)
	subsequentAmplification := SimulateSimple(sequences, targetTm, circular, append(primerList, initialAmplification...))
	if len(initialAmplification) != len(subsequentAmplification) {
		return initialAmplification, errors.New("Concatemerization detected in PCR.")
	}
	return initialAmplification, nil
}

func generatePcrFragments(sequence string, forwardLocation int, reverseLocation int, forwardPrimerIndxs []int, reversePrimerIndxs []int, minimalPrimers []string, primerList []string) []string {
	var pcrFragments []string
	for forwardPrimerIndex := range forwardPrimerIndxs {
		minimalPrimer := minimalPrimers[forwardPrimerIndex]
		fullPrimerForward := primerList[forwardPrimerIndex]
		for _, reversePrimerIndex := range reversePrimerIndxs {
			fullPrimerReverse := transform.ReverseComplement(primerList[reversePrimerIndex])
			pcrFragment := fullPrimerForward[:len(fullPrimerForward)-len(minimalPrimer)] + sequence[forwardLocation:reverseLocation] + fullPrimerReverse
			pcrFragments = append(pcrFragments, pcrFragment)
		}
	}
	return pcrFragments
}
