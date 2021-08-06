package pcr

import (
	"github.com/TimothyStiles/poly/primers"
	"github.com/TimothyStiles/poly/transform"
	"index/suffixarray"
	"sort"
	"strings"
)

// https://doi.org/10.1089/dna.1994.13.75
var minimalPrimerLength int = 15

func DesignPrimersWithOverhangs(sequence, forwardOverhang, reverseOverhang string, targetTm float64) (string, string) {
	sequence = strings.ToUpper(sequence)
	forwardPrimer := sequence[0:minimalPrimerLength]
	for i := 0; primers.MeltingTemp(forwardPrimer) < targetTm; i++ {
		forwardPrimer = sequence[0 : minimalPrimerLength+i]
	}
	reversePrimer := transform.ReverseComplement(sequence[len(sequence)-minimalPrimerLength:])
	for i := 0; primers.MeltingTemp(reversePrimer) < targetTm; i++ {
		reversePrimer = transform.ReverseComplement(sequence[len(sequence)-(minimalPrimerLength+i):])
	}

	// Add overhangs to primer
	forwardPrimer = forwardOverhang + forwardPrimer
	reversePrimer = transform.ReverseComplement(reverseOverhang) + reversePrimer

	return forwardPrimer, reversePrimer
}

func DesignPrimers(sequence string, targetTm float64) (string, string) {
	return DesignPrimersWithOverhangs(sequence, "", "", targetTm)
}

func Simulate(sequence string, targetTm float64, circular bool, primerList []string) []string {
	sequence = strings.ToUpper(sequence)
	sequenceIndex := suffixarray.New([]byte(sequence))

	primerLength := len(primerList)

	forwardLocations := make(map[int][]int)
	reverseLocations := make(map[int][]int)
	minimalPrimers := make([]string, primerLength)
	for primerIdx, primer := range primerList {
		var minimalLength int
		for i := 10; primers.MeltingTemp(primer[len(primer)-i:]) < targetTm; i++ {
			minimalLength = i
			if primer[len(primer)-i:] == primer {
				break
			}
		}
		// Use the minimal binding sites of the primer to find positions in the genome
		if primer[len(primer)-minimalLength:] != primer {
			minimalPrimer := primer[len(primer)-minimalLength:]
			minimalPrimers[primerIdx] = minimalPrimer
			for _, forwardLocation := range sequenceIndex.Lookup([]byte(minimalPrimer), -1) {
				forwardLocations[forwardLocation] = append(forwardLocations[forwardLocation], primerIdx)
			}
			for _, reverseLocation := range sequenceIndex.Lookup([]byte(transform.ReverseComplement(minimalPrimer)), -1) {
				reverseLocations[reverseLocation] = append(reverseLocations[reverseLocation], primerIdx)
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

	var pcrFragment string
	var pcrFragments []string

	// Next, iterate through the forwardLocations list
	for i, forwardLocation := range forwardLocationInts {
		// First, make sure that this isn't the last element in forwardLocations
		if i+1 != len(forwardLocationInts) {
			// If this isn't the last element in forwardLocations, then we can select the first reverseLocation that is less than the next forwardLocation
			for _, reverseLocation := range reverseLocationInts {
				if (forwardLocation < reverseLocation) && (reverseLocation < forwardLocationInts[i+1]) {
					// If both are true, we have found the sequence we are aiming to PCR! Now, we get all primers from that forwardLocation and then
					// build PCR fragments with each one.
					for forwardPrimerIdx := range forwardLocations[forwardLocation] {
						minimalPrimer := minimalPrimers[forwardPrimerIdx]
						fullPrimerForward := primerList[forwardPrimerIdx]
						for _, reversePrimerIdx := range reverseLocations[reverseLocation] {
							fullPrimerReverse := transform.ReverseComplement(primerList[reversePrimerIdx])
							pcrFragment = fullPrimerForward[:len(fullPrimerForward)-len(minimalPrimer)] + sequence[forwardLocation:reverseLocation] + fullPrimerReverse
							pcrFragments = append(pcrFragments, pcrFragment)
						}
					}
					break
				}
			}
		} else {
			foundFragment := false
			for _, reverseLocation := range reverseLocationInts {
				if forwardLocation < reverseLocation {
					for forwardPrimerIdx := range forwardLocations[forwardLocation] {
						minimalPrimer := minimalPrimers[forwardPrimerIdx]
						fullPrimerForward := primerList[forwardPrimerIdx]
						for _, reversePrimerIdx := range reverseLocations[reverseLocation] {
							fullPrimerReverse := transform.ReverseComplement(primerList[reversePrimerIdx])
							pcrFragment = fullPrimerForward[:len(fullPrimerForward)-len(minimalPrimer)] + sequence[forwardLocation:reverseLocation] + fullPrimerReverse
							pcrFragments = append(pcrFragments, pcrFragment)
						}
					}
					foundFragment = true
				}
			}
			// If the sequence is circular and we haven't found a fragment yet, check the other side of the origin
			if circular && !foundFragment {
				for _, reverseLocation := range reverseLocationInts {
					if forwardLocationInts[0] > reverseLocation {
						// If either one of these are true, create a new pcrFragment and append to pcrFragments
						for forwardPrimerIdx := range forwardLocations[forwardLocation] {
							minimalPrimer := minimalPrimers[forwardPrimerIdx]
							fullPrimerForward := primerList[forwardPrimerIdx]
							for _, reversePrimerIdx := range reverseLocations[reverseLocation] {
								fullPrimerReverse := transform.ReverseComplement(primerList[reversePrimerIdx])
								pcrFragment = fullPrimerForward[:len(fullPrimerForward)-len(minimalPrimer)] + sequence[forwardLocation:reverseLocation] + fullPrimerReverse
								pcrFragments = append(pcrFragments, pcrFragment)
							}
						}

					}
				}
			}
		}
	}
	return pcrFragments
}
