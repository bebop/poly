package pcr

import (
	"index/suffixarray"
	"math"
	"sort"
	"strings"

	"github.com/TimothyStiles/poly/checks"
	"github.com/TimothyStiles/poly/primers"
	"github.com/TimothyStiles/poly/transform"
)

// https://doi.org/10.1089/dna.1994.13.75
var minimalPrimerLength int = 15

type PairPrimers struct {
	Foward  string
	Reverse string
}

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

// Generally people that wants a pair of primers for a given sequence doesn't have a target Tm in mind.
// They want to put a sequence as input and receive one or more possible pairs of primers as output.
// So we will try make this by receive a given sequence and following some of the main guidelines for
// design of primers.

// One of the problems of Primer Design Guidelines is they're slightly different from each other.
// One example of these is Addgene Primer Design guidelines says the melting temperature should be between 50-60°C
// Thermo Fisher presents another information, they say the desired melting temperature is actually between 65°C and 75°C
// This another guideline from Premier Biosoft only says 'Primers with melting temperatures in the range of 52-58°C generally produce the best results'
// Sources: https://www.addgene.org/protocols/primer-design/, https://www.thermofisher.com/blog/behindthebench/pcr-primer-design-tips/
// http://www.premierbiosoft.com/tech_notes/PCR_Primer_Design.html

// So to allow users to follow your own guidelines we will make a software that founds every pair of primers
// using some ranges for melting temperature of each primer, GC content and delta melting temperature of a pair
// of primers
//
// So to make this happen I will use the following approach:
//
// 1. I will take a sequence and generate two list of possibles primers using ranges passed by user as arguments
// 2. I will use both lists to try to find pairs of primers that have the correct delta melting temperature between
// a pair of primers
// 3. Return all possible pair of primers

func DesignAllPossiblePrimers(sequence string, minGcContent float64, maxGcContent float64, minMeltingTemp float64, maxMeltingTemp float64, minLength int, maxLength int, maxMeltingTempBetweenPrimers float64) []PairPrimers {
	sequence = strings.ToUpper(sequence)

	checkPrimer := createCheckPrimerFunction(minGcContent, maxGcContent, minMeltingTemp, maxMeltingTemp, minLength, maxLength)

	var fowardPrimers []string
	for i := 0; i <= maxLength; i++ {
		forwardPrimer := sequence[0 : minLength+i]
		if checkPrimer(forwardPrimer) {
			fowardPrimers = append(fowardPrimers, forwardPrimer)
		}
	}

	var reversePrimers []string
	for i := 0; i < maxLength; i++ {
		reversePrimer := transform.ReverseComplement(sequence[len(sequence)-(minLength+i):])
		if checkPrimer(reversePrimer) {
			reversePrimers = append(reversePrimers, reversePrimer)
		}
	}

	pairOfPrimers := getPairOfPrimers(fowardPrimers, reversePrimers, maxMeltingTempBetweenPrimers)
	return pairOfPrimers
}

// Create a function that checks if a prime is following the guideline rules
func createCheckPrimerFunction(minGcContent float64, maxGcContent float64, minMeltingTemp float64, maxMeltingTemp float64, minLength int, maxLength int) func(string) bool {
	return func(primer string) bool {
		meltingTemp := primers.MeltingTemp(primer) / 100
		gcContent := checks.GcContent(primer)
		if len(primer) >= minLength && len(primer) <= maxLength &&
			meltingTemp >= minMeltingTemp && meltingTemp <= maxMeltingTemp &&
			gcContent >= minGcContent && gcContent <= maxGcContent {
			return true
		}
		//fmt.Println(len(primer) >= minLength, len(primer) <= maxLength, meltingTemp >= minMeltingTemp, meltingTemp <= minMeltingTemp, gcContent >= minGcContent, gcContent <= maxGcContent)
		return false
	}

}

func getPairOfPrimers(fowardPrimers []string, reversePrimers []string, maxMeltingTempBetweenPrimers float64) []PairPrimers {
	var pairOfPrimers []PairPrimers

	for _, fowardPrimer := range fowardPrimers {
		for _, reversePrimer := range reversePrimers {
			fowardTm := primers.MeltingTemp(fowardPrimer)
			reverseTm := primers.MeltingTemp(reversePrimer)
			if math.Abs(fowardTm-reverseTm) <= maxMeltingTempBetweenPrimers {
				pairOfPrimers = append(pairOfPrimers, PairPrimers{fowardPrimer, reversePrimer})
			}
		}
	}

	return pairOfPrimers
}
