package dataset

import (
	"fmt"
	"regexp"

	"github.com/TimothyStiles/poly"
	"github.com/TimothyStiles/poly/mfe"
)

/***********************
Functions to compute properties of mRNAs
***********************/

// The default properties to compute values of and add to DatasetWithProperties
// var DefaultProperies [](func(MRNA, map[string]interface{}) (interface{}, string)) = [](func(MRNA, map[string]interface{}) (interface{}, string)){
// 	// same,
// 	lenFivePrimeUTR,
// 	sdSequence,
// }
type Property struct {
	computePropertyFn  func(MRNA) interface{}
	includeInOutputCSV bool
}

// var DefaultProperies map[string]bool = map[string]bool{
// 	"LenFivePrimeUTR": true,
// 	"SDSequence":      true,
// 	"LenSDSequence":   true,
// }
var DefaultProperies []Property = []Property{
	{MRNA.lenFivePrimeUTR, true},
	{MRNA.hasFivePrimeUTR, false},
	{MRNA.lenProteinCodingSequence, false},
	{MRNA.hasProteinCodingSequence, false},

	{MRNA.sdBindingSite, false},
	{MRNA.hasSDBindingSite, false},
	{MRNA.lenSDSequence, true},
	{MRNA.mostFivePrimeSDBindingSite, true},
	{MRNA.mostThreePrimeSDBindingSite, true},

	{MRNA.codingSequenceCutoff, true},
	{MRNA.usedMRNA, true},
	{MRNA.usedMRNAStructure, true},
	{MRNA.dG_mRNA, true},
}

func toString(i interface{}) string {
	return fmt.Sprint(i)
}

func (mrna MRNA) lenFivePrimeUTR() interface{} {
	lenFivePrimeUTR := len(mrna.FivePrimeUTR)
	return lenFivePrimeUTR
}

func (mrna MRNA) hasFivePrimeUTR() interface{} {
	return mrna.properties["lenFivePrimeUTR"] != 0
}

func (mrna MRNA) lenProteinCodingSequence() interface{} {
	lenProteinCodingSequence := len(mrna.ProteinCodingSequence)
	return lenProteinCodingSequence
}

func (mrna MRNA) hasProteinCodingSequence() interface{} {
	return mrna.properties["lenProteinCodingSequence"] != 0
}

type SDBindingSite struct {
	mostFivePrimePairedOnMRNA, mostThreePrimePairedOnMRNA int
	sdSequence                                            string
	lenSDSequence                                         int
}

func (mrna MRNA) sdBindingSite() interface{} {
	fivePrimeUTR := mrna.FivePrimeUTR
	lenFivePrimeUTR := mrna.properties["lenFivePrimeUTR"].(int)
	// fivePrimeUTR = Reverse(fivePrimeUTR)
	fullAntiSDSequence := "ACCUCCUUA"
	fullAntiSDSequence = Reverse(fullAntiSDSequence)
	fullSDSequence := "UAAGGAGGU"
	lenFullSDSequence := len(fullSDSequence)

	// 18 is because SD sequence can end max 17 nucelotides from the start codon
	// and is at max 9 nucleotides long. We remove 2 because GAG is two away from
	// the start of the SD sequence
	smallestSDSequenceStartPos := 5
	smallestSDSequenceEndPos := 8
	smallestSDSequence := fullSDSequence[smallestSDSequenceStartPos-1 : smallestSDSequenceEndPos-1]
	lenSmallestSDSequence := len(smallestSDSequence)

	if lenFivePrimeUTR < lenSmallestSDSequence {
		return SDBindingSite{
			-1, -1, "", -1,
		}
	}

	maxDistFromStartCodonForSDSequece := 17 + lenFullSDSequence - smallestSDSequenceStartPos
	minDistFromStartCodonForSDSequence := 0

	// The sequences of the mRNA to search for the SD sequence in
	sdMRNASearchSequenceStartIdx := lenFivePrimeUTR - maxDistFromStartCodonForSDSequece
	if sdMRNASearchSequenceStartIdx < 0 {
		sdMRNASearchSequenceStartIdx = 0
	}
	sdMRNASearchSequenceEndPos := lenFivePrimeUTR - minDistFromStartCodonForSDSequence
	if sdMRNASearchSequenceEndPos < 1 {
		sdMRNASearchSequenceEndPos = 1
	} else if sdMRNASearchSequenceEndPos > lenFivePrimeUTR {
		sdMRNASearchSequenceEndPos = lenFivePrimeUTR - 1
	}
	sdMRNASearchSequence := fivePrimeUTR[sdMRNASearchSequenceStartIdx:sdMRNASearchSequenceEndPos]

	smallestSDSequenceRegex := regexp.MustCompile(smallestSDSequence)
	matches := smallestSDSequenceRegex.FindAllStringIndex(sdMRNASearchSequence, -1)
	if matches == nil {
		// No SD sequence in 5' untranslated region
		return SDBindingSite{
			-1, -1, "", -1,
		}
	} else {
		// Keep track of best (i.e. longest) SD sequence
		maxLenSDSequence := 0
		// var bestSDSequence string
		var bestSDBindingSite SDBindingSite

		for _, match := range matches {
			smallestSDSequenceStartIdx, _ := sdMRNASearchSequenceStartIdx+match[0], sdMRNASearchSequenceStartIdx+match[1]

			potentialFullSDSequenceStartIdx := smallestSDSequenceStartIdx - (smallestSDSequenceStartPos - 1)
			potentialFullSDSequenceEndIdx := potentialFullSDSequenceStartIdx + lenFullSDSequence
			if potentialFullSDSequenceEndIdx > lenFivePrimeUTR {
				potentialFullSDSequenceEndIdx = lenFivePrimeUTR
			}
			potentialFullSDSequence := fivePrimeUTR[potentialFullSDSequenceStartIdx:potentialFullSDSequenceEndIdx]

			sdSequenceStartIdx, sdSequenceEndIdx := findSDSequence(potentialFullSDSequence, fullAntiSDSequence, smallestSDSequenceStartPos, lenSmallestSDSequence)
			sdSequence := potentialFullSDSequence[sdSequenceStartIdx:sdSequenceEndIdx]
			lenSDSequence := len(sdSequence)
			// lenSDSequence, sdSequence := getFullSDSequence(fivePrimeUTR, smallestSDSequence, smallestSDSequenceStartPos, smallestSDSequenceStartIdx, smallestSDSequenceEndIdx)
			if lenSDSequence >= maxLenSDSequence {
				bestSDBindingSite = SDBindingSite{
					mostFivePrimePairedOnMRNA:  potentialFullSDSequenceStartIdx + sdSequenceStartIdx,
					mostThreePrimePairedOnMRNA: potentialFullSDSequenceStartIdx + sdSequenceEndIdx,
					lenSDSequence:              lenSDSequence,
					sdSequence:                 sdSequence,
				}
				maxLenSDSequence = lenSDSequence
			}
		}
		return bestSDBindingSite
	}
}

func (mrna MRNA) hasSDBindingSite() interface{} {
	sdBindingSite := mrna.properties["sdBindingSite"].(SDBindingSite)
	return sdBindingSite.mostFivePrimePairedOnMRNA != -1
}

func (mrna MRNA) lenSDSequence() interface{} {
	// panic(len(mrna.properties["SDSequence"].(string)))
	if mrna.properties["hasSDBindingSite"].(bool) {
		sdBindingSite := mrna.properties["sdBindingSite"].(SDBindingSite)
		return sdBindingSite.lenSDSequence
	} else {
		return 0
	}
}

func (mrna MRNA) mostFivePrimeSDBindingSite() interface{} {
	// panic(len(mrna.properties["SDSequence"].(string)))
	if mrna.properties["hasSDBindingSite"].(bool) {
		sdBindingSite := mrna.properties["sdBindingSite"].(SDBindingSite)
		return sdBindingSite.mostFivePrimePairedOnMRNA
	} else {
		return 0
	}
}

func (mrna MRNA) mostThreePrimeSDBindingSite() interface{} {
	// panic(len(mrna.properties["SDSequence"].(string)))
	if mrna.properties["hasSDBindingSite"].(bool) {
		sdBindingSite := mrna.properties["sdBindingSite"].(SDBindingSite)
		return sdBindingSite.mostThreePrimePairedOnMRNA
	} else {
		return 0
	}
}

var doesBasePair map[byte]map[byte]bool = map[byte]map[byte]bool{
	'A': {'U': true},
	'U': {'A': true, 'G': true},
	'C': {'G': true},
	'G': {'U': true, 'C': true},
}

func findSDSequence(potentialSDSequence, fullAntiSDSequence string, smallestSDSequenceStartPos, lenSmallestSDSequence int) (int, int) {
	sdSequenceStartIdx, sdSequenceEndIdx := smallestSDSequenceStartPos-1, smallestSDSequenceStartPos+lenSmallestSDSequence-1
	for i := sdSequenceStartIdx - 1; i >= 0; i-- {
		currSDSequenceChar := potentialSDSequence[i]
		currAntiSDSequenceChar := fullAntiSDSequence[i]
		if doesBasePair[currSDSequenceChar][currAntiSDSequenceChar] {
			sdSequenceStartIdx--
		} else {
			break
		}
	}

	for i := sdSequenceEndIdx + 1; i < len(potentialSDSequence); i++ {
		currSDSequenceChar := potentialSDSequence[i]
		currAntiSDSequenceChar := fullAntiSDSequence[i]
		if doesBasePair[currSDSequenceChar][currAntiSDSequenceChar] {
			sdSequenceEndIdx++
		} else {
			break
		}
	}
	return sdSequenceStartIdx, sdSequenceEndIdx
}

func (mrna MRNA) codingSequenceCutoff() interface{} {
	codingSequenceCutoff := 35
	if len(mrna.FivePrimeUTR) == 0 {
		codingSequenceCutoff = 76
	}

	var proteinCodingSequence string
	if mrna.properties["hasProteinCodingSequence"].(bool) {
		proteinCodingSequence = mrna.ProteinCodingSequence
	} else {
		proteinCodingSequence = mrna.Sequence
	}

	if codingSequenceCutoff > len(proteinCodingSequence) {
		codingSequenceCutoff = len(proteinCodingSequence)
	}
	return codingSequenceCutoff
}

func (mrna MRNA) usedMRNA() interface{} {
	var proteinCodingSequence string
	if mrna.properties["hasProteinCodingSequence"].(bool) {
		proteinCodingSequence = mrna.ProteinCodingSequence
	} else {
		proteinCodingSequence = mrna.Sequence
	}

	codingSequenceCutoff := mrna.properties["codingSequenceCutoff"].(int)
	return mrna.FivePrimeUTR + proteinCodingSequence[:codingSequenceCutoff]
}

func (mrna MRNA) usedMRNAStructure() interface{} {
	mRNA := mrna.properties["usedMRNA"].(string)
	if mRNA == "" {
		return ""
	} else {
		structure, _ := poly.LinearFold(mRNA)
		return structure
	}
}

func (mrna MRNA) dG_mRNA() interface{} {
	mRNA := mrna.properties["usedMRNA"].(string)
	mRNAStructure := mrna.properties["usedMRNAStructure"].(string)
	if mRNA == "" {
		return 0.0
	} else {
		dG_mRNA, _, err := mfe.MinimumFreeEnergy(mRNA, mRNAStructure, mfe.DefaultTemperature)
		if err != nil {
			panic(err)
		}
		return dG_mRNA
	}
}

/********************
Helper functions used to compute properties
*********************/

func Reverse(s string) string {
	runes := []rune(s)
	for i, j := 0, len(runes)-1; i < j; i, j = i+1, j-1 {
		runes[i], runes[j] = runes[j], runes[i]
	}
	return string(runes)
}

/****
Properties from original RBS calculator

What about dG_mRNA?
- Need to add ability to read in parameters into MinimumFreeEnergy
- Add ability for specifying no dangling ends

How do you add dG_mRNA_rRNA?


dG_standby
	1. Get list of modules (hairpins with single stranded up and down regions)
	2. For each module, make note of D, P, H and calculate best As
*****/
