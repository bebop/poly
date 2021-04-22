package poly

import (
	"github.com/jmoiron/sqlx"
	"github.com/juliangruber/go-intersect"
	"math/rand"
	"regexp"
	"strings"
	"sync"
	"time"

	weightedRand "github.com/mroth/weightedrand"

	"encoding/json"
	"io/ioutil"
)

/******************************************************************************
Oct, 15, 2020

File is structured as so:

	Structs:
		CodonTable - holds all information mapping codons <-> aminoacids during transformations.
		AnimoAcid - holds amino acide related info for CodonTable struct
		Codon - holds codon related info for AminoAcid struct

	Big functions that everything else is related to:

		Translate - given a nucleic sequence string and codon table it translates sequences
								to UPPERCASE amino acid sequences.

		Optimize - given an amino acid sequence string and codon table it translates
							 sequences to UPPERCASE nucleic acid sequences.

This file contains almost everything you need to do standard codon optimization.

Biological context: certain cells favor certain codons and will reject or under
express sequences that don't use a similar ratio of codons.
This is called codon bias: https://en.wikipedia.org/wiki/Codon_usage_bias

Furthermore, different ribosomes in different organisms will interpret codons differently.
What may be a start codon for one ribosome may be a stop in the other.
Heck, apparently nucleomorphs contain 4 different kinds of ribosomes.
https://en.wikipedia.org/wiki/Nucleomorph <- Thanks Keoni for mentioning this example!

Anywho, most of this file and CodonTable's struct methods are meant to help overcome
this codon bias. There's a default CodonTable generator near the bottom of this file
with a whole section on how it works and why it's gotta be that way.

Like most codebases, best usage examples come from tests. You can check out
TestTranslate and TestOptimize in transformations_test.go for pretty solid
examples of what you can do with this code. Will add more to docs site before merge
into prime.

TTFN,
Tim

******************************************************************************/

// Codon holds information for a codon triplet in a struct
type Codon struct {
	Triplet string `json:"triplet"`
	Weight  int    `json:"weight"` // needs to be set to 1 for random chooser
}

// AminoAcid holds information for an amino acid and related codons in a struct
type AminoAcid struct {
	Letter string  `json:"letter"`
	Codons []Codon `json:"codons"`
}

// CodonTable holds information for a codon table.
type CodonTable struct {
	StartCodons []string    `json:"start_codons"`
	StopCodons  []string    `json:"stop_codons"`
	AminoAcids  []AminoAcid `json:"amino_acids"`
}

// Translate translates a codon sequence to an amino acid sequence
func Translate(sequence string, codonTable CodonTable) string {

	var aminoAcids strings.Builder
	var currentCodon strings.Builder
	translationTable := codonTable.generateTranslationTable()

	for _, letter := range sequence {

		// add current nucleotide to currentCodon
		currentCodon.WriteRune(letter)

		// if current nucleotide is the third in a codon translate to aminoAcid write to aminoAcids and reset currentCodon.
		if currentCodon.Len() == 3 {
			aminoAcids.WriteString(translationTable[currentCodon.String()])

			// reset codon string builder for next codon.
			currentCodon.Reset()
		}
	}
	return aminoAcids.String()
}

// Optimize takes an amino acid sequence and CodonTable and returns an optimized codon sequence
func Optimize(aminoAcids string, codonTable CodonTable) string {

	// weightedRand library insisted setting seed like this. Not sure what environmental side effects exist.
	rand.Seed(time.Now().UTC().UnixNano())

	var codons strings.Builder
	codonChooser := codonTable.chooser()

	for _, aminoAcid := range aminoAcids {
		aminoAcidString := string(aminoAcid)
		codons.WriteString(codonChooser[aminoAcidString].Pick().(string))
	}
	return codons.String()
}

// GetOptimizationTable is a Sequence method that takes a CodonTable and weights it to be used to optimize inserts.
func (sequence Sequence) GetOptimizationTable(codonTable CodonTable) CodonTable {
	sequenceString := getCodingRegions(sequence)
	return codonTable.OptimizeTable(sequenceString)
}

// OptimizeTable weights each codon in a codon table according to input string codon frequency.
// This function actually mutates the CodonTable struct itself.
func (codonTable CodonTable) OptimizeTable(sequence string) CodonTable {

	sequence = strings.ToUpper(sequence)
	codonFrequencyMap := getCodonFrequency(sequence)

	for aminoAcidIndex, aminoAcid := range codonTable.AminoAcids {
		// apply weights to codonTable
		for codonIndex, codon := range aminoAcid.Codons {
			codonTable.AminoAcids[aminoAcidIndex].Codons[codonIndex].Weight = codonFrequencyMap[codon.Triplet]
		}

	}
	return codonTable
}

// getCodonFrequency takes a DNA sequence and returns a hashmap of its codons and their frequencies.
func getCodonFrequency(sequence string) map[string]int {

	codonFrequencyHashMap := map[string]int{}
	var currentCodon strings.Builder

	for _, letter := range sequence {

		// add current nucleotide to currentCodon
		currentCodon.WriteRune(letter)

		// if current nucleotide is the third in a codon add to hashmap
		if currentCodon.Len() == 3 {
			// if codon is already initalized in map increment
			if _, ok := codonFrequencyHashMap[currentCodon.String()]; ok {
				codonString := currentCodon.String()
				codonFrequencyHashMap[codonString]++
				// if codon is not already initalized in map initialize with value at 1
			} else {
				codonString := currentCodon.String()
				codonFrequencyHashMap[codonString] = 1
			}
			// reset codon string builder for next codon.
			currentCodon.Reset()
		}
	}
	return codonFrequencyHashMap
}

// chooser is a CodonTable method to convert a codon table to a chooser
func (codonTable CodonTable) chooser() map[string]weightedRand.Chooser {

	// This maps codon tables structure to weightRand.NewChooser structure
	codonChooser := make(map[string]weightedRand.Chooser)

	// iterate over every amino acid in the codonTable
	for _, aminoAcid := range codonTable.AminoAcids {

		// create a list of codon choices for this specific amino acid
		codonChoices := make([]weightedRand.Choice, len(aminoAcid.Codons))

		// Get sum of codon occurences for particular amino acid
		codonOccurenceSum := 0
		for _, codon := range aminoAcid.Codons {
			codonOccurenceSum += codon.Weight
		}

		// Threshold codons that occur less than 10% for coding a particular amino acid
		for _, codon := range aminoAcid.Codons {
			codonPercentage := float64(codon.Weight) / float64(codonOccurenceSum)

			if codonPercentage > 0.10 {
				// for every codon related to current amino acid append its Triplet and Weight to codonChoices after thresholding
				codonChoices = append(codonChoices, weightedRand.Choice{Item: codon.Triplet, Weight: uint(codon.Weight)})
			}
		}

		// add this chooser set to the codonChooser map under the name of the aminoAcid it represents.
		codonChooser[aminoAcid.Letter] = weightedRand.NewChooser(codonChoices...)
	}
	return codonChooser
}

// Generate map of codons -> amino acid
func (codonTable CodonTable) generateTranslationTable() map[string]string {
	var translationMap = make(map[string]string)
	for _, aminoAcid := range codonTable.AminoAcids {
		for _, codon := range aminoAcid.Codons {
			translationMap[codon.Triplet] = aminoAcid.Letter
		}
	}
	return translationMap
}

// helper function to pull coding regions out of an Sequence
func getCodingRegions(sequence Sequence) string {
	// pick out the each coding region in the Sequence and add it to the sequence Builder
	var sequenceBuilder strings.Builder

	for _, feature := range sequence.Features {
		if feature.Type == "CDS" {
			sequenceBuilder.WriteString(feature.GetSequence())
		}
	}

	return sequenceBuilder.String()
}

/******************************************************************************
Oct, 15, 2020

Codon table generation stuff begins here.

Alright, I know it's ugly down below this comment block but there ain't much
we can do until we can experimentally derive our own codon tables.

The story is this. Different organisms use different codons to represent
different things.

The NCBI publishes this weird data format for developers to use for generating
codon tables and mapping codons to amino acids for different organisms.

All this stuff is experimentally derived and I'm not sure how it's done really.
I won't really have a chance to find out for a while but there's some future
work where I may want to do experiments like this and you'll see more about it.

There are two tables. I got annoyed since the original only went by number so
I made one that went by name too. Looking back on it this is probably useless
but maybe someone will find a use for it so here it stays.

Happy hacking,
Tim

P.S

Maybe we should publish our JSON representations? Let me know if want you to
organize that.

******************************************************************************/

// Function to generate default codon tables from NCBI https://www.ncbi.nlm.nih.gov/Taxonomy/Utils/wprintgc.cgi
func generateCodonTable(aminoAcids, starts string) CodonTable {
	base1 := "TTTTTTTTTTTTTTTTCCCCCCCCCCCCCCCCAAAAAAAAAAAAAAAAGGGGGGGGGGGGGGGG"
	base2 := "TTTTCCCCAAAAGGGGTTTTCCCCAAAAGGGGTTTTCCCCAAAAGGGGTTTTCCCCAAAAGGGG"
	base3 := "TCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAG"
	// Add triplets to an amino acid -> triplet map, and if a possible start codon, add to start codon list
	var aminoAcidMap = make(map[rune][]Codon)
	var startCodons []string
	var stopCodons []string
	for i, aminoAcid := range aminoAcids {
		if _, ok := aminoAcidMap[aminoAcid]; ok == false {
			aminoAcidMap[aminoAcid] = []Codon{}
		}
		triplet := string([]byte{base1[i], base2[i], base3[i]})
		aminoAcidMap[aminoAcid] = append(aminoAcidMap[aminoAcid], Codon{triplet, 1})
		if starts[i] == 77 { // M rune
			startCodons = append(startCodons, triplet)
		}
		if starts[i] == 42 { // * rune
			stopCodons = append(stopCodons, triplet)
		}
	}
	// Convert amino acid -> triplet map to an amino acid list
	var aminoAcidSlice []AminoAcid
	for k, v := range aminoAcidMap {
		aminoAcidSlice = append(aminoAcidSlice, AminoAcid{string(k), v})
	}
	return CodonTable{startCodons, stopCodons, aminoAcidSlice}
}

// GetCodonTable takes the index of desired NCBI codon table and returns it.
func GetCodonTable(index int) CodonTable {
	return defaultCodonTablesByNumber[index]
}

// defaultCodonTablesByNumber stores all codon tables published by NCBI https://www.ncbi.nlm.nih.gov/Taxonomy/Utils/wprintgc.cgi using numbered indeces.
var defaultCodonTablesByNumber = map[int]CodonTable{
	1:  generateCodonTable("FFLLSSSSYY**CC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG", "---M------**--*----M---------------M----------------------------"),
	2:  generateCodonTable("FFLLSSSSYY**CCWWLLLLPPPPHHQQRRRRIIMMTTTTNNKKSS**VVVVAAAADDEEGGGG", "----------**--------------------MMMM----------**---M------------"),
	3:  generateCodonTable("FFLLSSSSYY**CCWWTTTTPPPPHHQQRRRRIIMMTTTTNNKKSSRRVVVVAAAADDEEGGGG", "----------**----------------------MM---------------M------------"),
	4:  generateCodonTable("FFLLSSSSYY**CCWWLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG", "--MM------**-------M------------MMMM---------------M------------"),
	5:  generateCodonTable("FFLLSSSSYY**CCWWLLLLPPPPHHQQRRRRIIMMTTTTNNKKSSSSVVVVAAAADDEEGGGG", "---M------**--------------------MMMM---------------M------------"),
	6:  generateCodonTable("FFLLSSSSYYQQCC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG", "--------------*--------------------M----------------------------"),
	9:  generateCodonTable("FFLLSSSSYY**CCWWLLLLPPPPHHQQRRRRIIIMTTTTNNNKSSSSVVVVAAAADDEEGGGG", "----------**-----------------------M---------------M------------"),
	10: generateCodonTable("FFLLSSSSYY**CCCWLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG", "----------**-----------------------M----------------------------"),
	11: generateCodonTable("FFLLSSSSYY**CC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG", "---M------**--*----M------------MMMM---------------M------------"),
	12: generateCodonTable("FFLLSSSSYY**CC*WLLLSPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG", "----------**--*----M---------------M----------------------------"),
	13: generateCodonTable("FFLLSSSSYY**CCWWLLLLPPPPHHQQRRRRIIMMTTTTNNKKSSGGVVVVAAAADDEEGGGG", "---M------**----------------------MM---------------M------------"),
	14: generateCodonTable("FFLLSSSSYYY*CCWWLLLLPPPPHHQQRRRRIIIMTTTTNNNKSSSSVVVVAAAADDEEGGGG", "-----------*-----------------------M----------------------------"),
	16: generateCodonTable("FFLLSSSSYY*LCC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG", "----------*---*--------------------M----------------------------"),
	21: generateCodonTable("FFLLSSSSYY**CCWWLLLLPPPPHHQQRRRRIIMMTTTTNNNKSSSSVVVVAAAADDEEGGGG", "----------**-----------------------M---------------M------------"),
	22: generateCodonTable("FFLLSS*SYY*LCC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG", "------*---*---*--------------------M----------------------------"),
	23: generateCodonTable("FF*LSSSSYY**CC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG", "--*-------**--*-----------------M--M---------------M------------"),
	24: generateCodonTable("FFLLSSSSYY**CCWWLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSSKVVVVAAAADDEEGGGG", "---M------**-------M---------------M---------------M------------"),
	25: generateCodonTable("FFLLSSSSYY**CCGWLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG", "---M------**-----------------------M---------------M------------"),
	26: generateCodonTable("FFLLSSSSYY**CC*WLLLAPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG", "----------**--*----M---------------M----------------------------"),
	27: generateCodonTable("FFLLSSSSYYQQCCWWLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG", "--------------*--------------------M----------------------------"),
	28: generateCodonTable("FFLLSSSSYYQQCCWWLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG", "----------**--*--------------------M----------------------------"),
	29: generateCodonTable("FFLLSSSSYYYYCC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG", "--------------*--------------------M----------------------------"),
	30: generateCodonTable("FFLLSSSSYYEECC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG", "--------------*--------------------M----------------------------"),
	31: generateCodonTable("FFLLSSSSYYEECCWWLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG", "----------**-----------------------M----------------------------"),
	33: generateCodonTable("FFLLSSSSYYY*CCWWLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSSKVVVVAAAADDEEGGGG", "---M-------*-------M---------------M---------------M------------")}

// DefaultCodonTablesByName stores all codon tables published by NCBI https://www.ncbi.nlm.nih.gov/Taxonomy/Utils/wprintgc.cgi using named indeces.
var defaultCodonTablesByName = map[string]CodonTable{
	"Standard":                         generateCodonTable("FFLLSSSSYY**CC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG", "---M------**--*----M---------------M----------------------------"), // 1
	"VertebrateMitochondrial":          generateCodonTable("FFLLSSSSYY**CCWWLLLLPPPPHHQQRRRRIIMMTTTTNNKKSS**VVVVAAAADDEEGGGG", "----------**--------------------MMMM----------**---M------------"), // 2
	"YeastMitochondrial":               generateCodonTable("FFLLSSSSYY**CCWWTTTTPPPPHHQQRRRRIIMMTTTTNNKKSSRRVVVVAAAADDEEGGGG", "----------**----------------------MM---------------M------------"), // 3
	"MoldMitochondrial":                generateCodonTable("FFLLSSSSYY**CCWWLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG", "--MM------**-------M------------MMMM---------------M------------"), // 4
	"ProtozoanMitochondrial":           generateCodonTable("FFLLSSSSYY**CCWWLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG", "--MM------**-------M------------MMMM---------------M------------"), // 4
	"CoelenterateMitochondrial":        generateCodonTable("FFLLSSSSYY**CCWWLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG", "--MM------**-------M------------MMMM---------------M------------"), // 4
	"Mycoplasma":                       generateCodonTable("FFLLSSSSYY**CCWWLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG", "--MM------**-------M------------MMMM---------------M------------"), // 4
	"Spiroplasma":                      generateCodonTable("FFLLSSSSYY**CCWWLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG", "--MM------**-------M------------MMMM---------------M------------"), // 4
	"InvertebrateMitochondrial":        generateCodonTable("FFLLSSSSYY**CCWWLLLLPPPPHHQQRRRRIIMMTTTTNNKKSSSSVVVVAAAADDEEGGGG", "---M------**--------------------MMMM---------------M------------"), // 5
	"CiliateNuclear":                   generateCodonTable("FFLLSSSSYYQQCC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG", "--------------*--------------------M----------------------------"), // 6
	"DasycladaceanNuclear":             generateCodonTable("FFLLSSSSYYQQCC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG", "--------------*--------------------M----------------------------"), // 6
	"HexamitaNuclear":                  generateCodonTable("FFLLSSSSYYQQCC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG", "--------------*--------------------M----------------------------"), // 6
	"EchinodermMitochondrial":          generateCodonTable("FFLLSSSSYY**CCWWLLLLPPPPHHQQRRRRIIIMTTTTNNNKSSSSVVVVAAAADDEEGGGG", "----------**-----------------------M---------------M------------"), // 9
	"FlatwormMitochondrial":            generateCodonTable("FFLLSSSSYY**CCWWLLLLPPPPHHQQRRRRIIIMTTTTNNNKSSSSVVVVAAAADDEEGGGG", "----------**-----------------------M---------------M------------"), // 9
	"EuplotidNuclear":                  generateCodonTable("FFLLSSSSYY**CCCWLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG", "----------**-----------------------M----------------------------"), // 10
	"BacterialPlastid":                 generateCodonTable("FFLLSSSSYY**CC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG", "---M------**--*----M------------MMMM---------------M------------"), // 11
	"ArchaelPlastid":                   generateCodonTable("FFLLSSSSYY**CC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG", "---M------**--*----M------------MMMM---------------M------------"), // 11
	"PlantPlastid":                     generateCodonTable("FFLLSSSSYY**CC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG", "---M------**--*----M------------MMMM---------------M------------"), // 11
	"AlternativeYeastNuclear":          generateCodonTable("FFLLSSSSYY**CC*WLLLSPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG", "----------**--*----M---------------M----------------------------"), // 12
	"AscidianMitochondrial":            generateCodonTable("FFLLSSSSYY**CCWWLLLLPPPPHHQQRRRRIIMMTTTTNNKKSSGGVVVVAAAADDEEGGGG", "---M------**----------------------MM---------------M------------"), // 13
	"AlternativeFlatwormMitochondrial": generateCodonTable("FFLLSSSSYYY*CCWWLLLLPPPPHHQQRRRRIIIMTTTTNNNKSSSSVVVVAAAADDEEGGGG", "-----------*-----------------------M----------------------------"), // 14
	"ChlorophyceanMitochondrial":       generateCodonTable("FFLLSSSSYY*LCC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG", "----------*---*--------------------M----------------------------"), // 16
	"TrematodeMitochondrial":           generateCodonTable("FFLLSSSSYY**CCWWLLLLPPPPHHQQRRRRIIMMTTTTNNNKSSSSVVVVAAAADDEEGGGG", "----------**-----------------------M---------------M------------"), // 21
	"ScenedesmusObliquusMitochondrial": generateCodonTable("FFLLSS*SYY*LCC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG", "------*---*---*--------------------M----------------------------"), // 22
	"ThraustochytriumMitochondrial":    generateCodonTable("FF*LSSSSYY**CC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG", "--*-------**--*-----------------M--M---------------M------------"), // 23
	"RhabdopleuridaeMitochondrial":     generateCodonTable("FFLLSSSSYY**CCWWLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSSKVVVVAAAADDEEGGGG", "---M------**-------M---------------M---------------M------------"), // 24
	"CandidateDivisionSR1":             generateCodonTable("FFLLSSSSYY**CCGWLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG", "---M------**-----------------------M---------------M------------"), // 25
	"Gracilibacteria":                  generateCodonTable("FFLLSSSSYY**CCGWLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG", "---M------**-----------------------M---------------M------------"), // 25
	"PachysolenTannophilusNuclear":     generateCodonTable("FFLLSSSSYY**CC*WLLLAPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG", "----------**--*----M---------------M----------------------------"), // 26
	"KaryorelictNuclear":               generateCodonTable("FFLLSSSSYYQQCCWWLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG", "--------------*--------------------M----------------------------"), // 27
	"CondylostomaNuclear":              generateCodonTable("FFLLSSSSYYQQCCWWLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG", "----------**--*--------------------M----------------------------"), // 28
	"MesodiniumNuclear":                generateCodonTable("FFLLSSSSYYYYCC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG", "--------------*--------------------M----------------------------"), // 29
	"PeritrichNuclear":                 generateCodonTable("FFLLSSSSYYEECC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG", "--------------*--------------------M----------------------------"), // 30
	"BlastocrithidiaNuclear":           generateCodonTable("FFLLSSSSYYEECCWWLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG", "----------**-----------------------M----------------------------"), // 31
	"CephalodiscidaeMitochondrial":     generateCodonTable("FFLLSSSSYYY*CCWWLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSSKVVVVAAAADDEEGGGG", "---M-------*-------M---------------M---------------M------------")} // 33

/******************************************************************************
Nov, 20, 2020

Codon table JSON stuff begins here.

In 2007, a Japanese group made a codon usage database based off of the current
iteration of GenBank at the time. It became pretty much the #1 place to get
codon tables for the last 13 years (https://www.kazusa.or.jp/codon/).

Here is an example of how they are usually stored:

(fields: [triplet] [frequency: per thousand] ([number]))
```
UUU 17.6(714298)  UCU 15.2(618711)  UAU 12.2(495699)  UGU 10.6(430311)
UUC 20.3(824692)  UCC 17.7(718892)  UAC 15.3(622407)  UGC 12.6(513028)
UUA  7.7(311881)  UCA 12.2(496448)  UAA  1.0( 40285)  UGA  1.6( 63237)
UUG 12.9(525688)  UCG  4.4(179419)  UAG  0.8( 32109)  UGG 13.2(535595)

CUU 13.2(536515)  CCU 17.5(713233)  CAU 10.9(441711)  CGU  4.5(184609)
CUC 19.6(796638)  CCC 19.8(804620)  CAC 15.1(613713)  CGC 10.4(423516)
CUA  7.2(290751)  CCA 16.9(688038)  CAA 12.3(501911)  CGA  6.2(250760)
CUG 39.6(1611801)  CCG  6.9(281570)  CAG 34.2(1391973)  CGG 11.4(464485)

AUU 16.0(650473)  ACU 13.1(533609)  AAU 17.0(689701)  AGU 12.1(493429)
AUC 20.8(846466)  ACC 18.9(768147)  AAC 19.1(776603)  AGC 19.5(791383)
AUA  7.5(304565)  ACA 15.1(614523)  AAA 24.4(993621)  AGA 12.2(494682)
AUG 22.0(896005)  ACG  6.1(246105)  AAG 31.9(1295568)  AGG 12.0(486463)

GUU 11.0(448607)  GCU 18.4(750096)  GAU 21.8(885429)  GGU 10.8(437126)
GUC 14.5(588138)  GCC 27.7(1127679)  GAC 25.1(1020595)  GGC 22.2(903565)
GUA  7.1(287712)  GCA 15.8(643471)  GAA 29.0(1177632)  GGA 16.5(669873)
GUG 28.1(1143534)  GCG  7.4(299495)  GAG 39.6(1609975)  GGG 16.5(669768)
```

You'll notice a couple of things here - Namely, this format isn't very
amenable to scripting (non-standardized IO format), and that the table
data (What is the amino acid to codon pairing?) has to be stored else-
where.

The database hasn't been updated for 13 years. The format isn't nice
for automation or for bulk analysis. We can do better. Poly's codonTable
format is a basic JSON file that can be used in different programs, and
since we have a nice GenBank parser, we can continuously run codon table
analysis on the most up-to-date GenBank files.

I just bought codontable.com and codontables.com, so hopefully we can
set up a webservice that is better than https://www.kazusa.or.jp/codon/

Also, I need codon tables in JSON for the codon optimizer app. That is
the real reason I'm doing this now.

Jolly Good!
Keoni

PS: On it @Tim, publishing those JSON representations. Woohoo!
******************************************************************************/

// ParseCodonJSON parses a CodonTable JSON file.
func ParseCodonJSON(file []byte) CodonTable {
	var codontable CodonTable
	_ = json.Unmarshal([]byte(file), &codontable)
	return codontable
}

// ReadCodonJSON reads a CodonTable JSON file.
func ReadCodonJSON(path string) CodonTable {
	file, err := ioutil.ReadFile(path)
	if err != nil {
		// return 0, fmt.Errorf("Failed to open file %s for unpack: %s", gzFilePath, err)
	}
	codontable := ParseCodonJSON(file)
	return codontable
}

// WriteCodonJSON writes a CodonTable struct out to JSON.
func WriteCodonJSON(codontable CodonTable, path string) {
	file, _ := json.MarshalIndent(codontable, "", " ")
	_ = ioutil.WriteFile(path, file, 0644)
}

/******************************************************************************
Dec, 9, 2020

Synthesis fixing stuff begins here.

Codon optimization is not good enough for synthetic DNA - there are certain
elements that synthesis doesn't really work for. For example, long homopolymers
will nearly always need removal in order to synthesize a gene efficiently.

Poly codon optimization + synthesis fixing is meant to be as efficient yet
accurate as possible.

# How do we do it?

From experience using Twist, there are a couple of major elements to fix:

1. Global GC content
2. ~50bp GC hotspots
3. Homopolymers
4. K-mer repeats
5. Various banned restriction enzyme sites

To implement this functionality, I imagine each major element will have a function
that takes in a sequence and returns 2 integers and either "GC", "AT", or "NA". The
2 integers are start and end positions of problematic sequences, and the "GC"/"AT"
are for if the sequence needs to be biased in a particular direction (we call this
the change type).

If there are two elements that have an overlap (for example, a high "GC" region
with a restriction enzyme site), fixing the outer element is tabled until the next
round of fixes. However, the nearest outer element change type will be inherited for
biased stochastic selection of a new codon. With the example of high "GC" and an
internal restriction enzyme site, the internal restriction enzyme will be changed
with an AT bias if possible.

Once all change sites are identified, changes are implemented across the entire
sequence. If any given change causes a problem which is *smaller* than the change,
that change is reverted, the entire context of the problematic region is noted, and
a change at a different position is attempted, until there is a change that fixes
the original problem and does not cause a problem which is *smaller* than that change.

The procedure above is repeated up the chain until all problems are resolved. If a problem
cannot be resolved, this is noted in the output, but it is not an error - "It is possible to
commit no mistakes and still lose. That is not a weakness. That is life."

[Removed notes about backtracking. It's hard and I don't really want to implement it right now]

Although there is lots of recursion and rechecking, I roughly estimate only 1~2 rounds will be
necessary for ~95% of sequences, which should be rather fast.

# Anything else?

This functionality should be fundamentally valuable to companies and labs. We will also
be able to improve the functionality of these fixes as we collect more data, or even
have company-by-company recommendations on things to fix. Hopefully by working in the
open, we can improve codon optimization and usage in a transparent way. One day,
this will enable us to be far better at expressing proteins at certain levels within
cells. That is the dream - working in the open will enable us to be far better at
engineering life.

Information wants to be free,

Keoni Gandall







This is Keoni coming back after a few months. Let me re-explain what this code will do.

1 First, we generate an in-memory SQLite DB. This database lets us do relational queries
  on important parts of our sequence that would be difficult to write in Golang
2 Next, we use the inputed functions to look for problematic regions (problematicRegion).
  Each problematicRegion has a list of suggested fixes (suggestedFixes) that, if fixed up,
  should resolve the problematicRegion.
3 Once we have the list of problematicRegions + suggestedFixes, we sort the suggestedFixes by
  the estimated amount of changes needed to fix that region.
3a The suggestedFix with the smallest quantity of changes necessary to fix the problem is then
   fixed in sequence space.
3b The specific codon that is changed to fix a problem is first sorted by if it has been used already.
   For example, if one codon has already been switched around, the system will avoid changing it
   again unless absolutely necessary.
3c Codons are then sorted to move towards optimal codons. For example, if one codon can be changed from
   a fairly rare codon to the most common codon for that amino acid, it will be chosen to be changed over
   a codon who is already fairly optimal.
4. All problematicRegions *without overlaps* are repeated with [GOTO 3] in order of quantity of changes
5. After all independent problematicRegions have been checked, repeat search with [GOTO 2]
6. If we end up with a sequence without any problematicRegions, we return the sequence!



******************************************************************************/

type DnaSuggestion struct {
	Start          int    `db:"start"`
	End            int    `db:"end"`
	Bias           string `db:"bias"`
	QuantityFixes  int    `db:"quantityfixes"`
	SuggestionType string `db:"suggestiontype"`
}

func FindBsaI(sequence string, c chan DnaSuggestion, wg *sync.WaitGroup) {
	defer wg.Done()

	re := regexp.MustCompile(`GGTCTC`)
	locs := re.FindAllStringIndex(sequence, -1)
	for _, loc := range locs {
		position := loc[0] / 3
		leftover := loc[0] % 3
		switch {
		case leftover == 0:
			c <- DnaSuggestion{position, (loc[1] / 3), "NA", 1, "BsaI removal"}
		case leftover != 0:
			c <- DnaSuggestion{position - 1, (loc[1] / 3) - 1, "NA", 1, "BsaI removal"}
		}
	}
}

func findProblems(sequence string, problematicSequenceFuncs []func(string, chan DnaSuggestion, *sync.WaitGroup)) []DnaSuggestion {
	// Run functions to get suggestions
	suggestions := make(chan DnaSuggestion, 100)
	var wg sync.WaitGroup
	for _, f := range problematicSequenceFuncs {
		wg.Add(1)
		go f(sequence, suggestions, &wg)
	}
	wg.Wait()

	var suggestionsList []DnaSuggestion
	for suggestion := range suggestions {
		suggestionsList = append(suggestionsList, suggestion)
	}
	return suggestionsList
}

func FixCds(sequence string, protectedPositions []int, codontable CodonTable, problematicSequenceFuncs []func(string, chan DnaSuggestion, *sync.WaitGroup)) (string, error) {
	var db *sqlx.DB
	db = sqlx.MustConnect("sqlite3", ":memory:")
	createMemoryDbSql := `
	CREATE TABLE codon (
		codon TEXT PRIMARY KEY,
		aa TEXT
	);

	CREATE TABLE seq (
		pos INT PRIMARY KEY,
		protected BOOL DEFAULT false
	);

	CREATE TABLE history (
		pos INT REFERENCES seq(pos),
		codon TEXT NOT NULL REFERNECES codon(codon),
		step INT
	);

	-- Weights are set on a per position basis for codon harmonization at a later point
	CREATE TABLE weights (
		pos INT REFERENCES seq(pos),
		codon TEXT NOT NULL REFERENCES codon(codon),
		weight INT
	);

	CREATE TABLE codonbias (
		fromcodon TEXT REFERENCES codon(codon),
		tocodon TEXT REFERENCES codon(codon),
		gcbias string
	);

	CREATE TABLE suggestedFix (
		id INT PRIMARY KEY,
		step INT,
		start INT REFERENCES seq(pos),
		end INT REFERENCES seq(pos),
		gcbias string,
		quantityfixes INT,
		suggestiontype TEXT
	);
`
	db.MustExec(createMemoryDbSql)
	// Insert codons
	weightTable := make(map[string]int)
	codonInsert := `INSERT INTO codon(codon, aa) VALUES (?, ?)`
	for _, aminoAcid := range codontable.AminoAcids {
		for _, codon := range aminoAcid.Codons {
			db.MustExec(codonInsert, codon.Triplet, aminoAcid.Letter)
			weightTable[codon.Triplet] = codon.Weight

			codonBias := strings.Count(codon.Triplet, "G") + strings.Count(codon.Triplet, "C")
			for _, toCodon := range aminoAcid.Codons {
				toCodonBias := strings.Count(toCodon.Triplet, "G") + strings.Count(toCodon.Triplet, "C")
				switch {
				case codonBias == toCodonBias:
					db.MustExec(`INSERT INTO codonbias(fromcodon, tocodon, gcbias) VALUES (?, ?, ?)`, codon.Triplet, toCodon.Triplet, "NA")
				case codonBias > toCodonBias:
					db.MustExec(`INSERT INTO codonbias(fromcodon, tocodon, gcbias) VALUES (?, ?, ?)`, codon.Triplet, toCodon.Triplet, "AT")
				case codonBias < toCodonBias:
					db.MustExec(`INSERT INTO codonbias(fromcodon, tocodon, gcbias) VALUES (?, ?, ?)`, codon.Triplet, toCodon.Triplet, "GC")
				}
			}
		}
	}

	// Insert seq and history
	pos := 1
	var protected bool
	for i := 0; i < len(sequence)-3; i = i + 3 {
		codon := sequence[i : i+3]
		protected = false
		for _, protectedPosition := range protectedPositions {
			if protectedPosition == i {
				protected = true
			}
		}
		db.MustExec(`INSERT INTO seq(pos, protected) VALUES (?, ?)`, pos, protected)
		db.MustExec(`INSERT INTO history(pos, codon, step) VALUES (?, ?, 0)`, pos, codon)
		db.MustExec(`INSERT INTO weights(pos, codon, weight) VALUES (?,?,?)`, pos, codon, weightTable[codon])
	}

	// For a maximum of 1000 iterations, see if we can do better
	for i := 1; i < 1000; i++ {
		suggestions := findProblems(sequence, problematicSequenceFuncs)
		// If there are no suggestions, break the iteration!
		if len(suggestions) == 0 {
			break
		}
		for suggestionIndex, suggestion := range suggestions {
			// First, let's insert the suggestions that we found using our problematicSequenceFuncs
			db.MustExec(`INSERT INTO suggestedfix(step, start, end, gcbias, quantityfixes) VALUES (?, ?, ?, ?, ?)`, i, suggestion.Start, suggestion.End, suggestion.Bias, suggestion.QuantityFixes)

			// Second, lets look for if there are any overlapping suggestions. This is equivalent to a spot where we can "kill two birds with one stone"
			// TODO: make this into a function so you can get overlapping overlaps (for when 3 problems are overlapped on top of each other)
			for secondSuggestionIndex, secondSuggestion := range suggestions {
				if suggestionIndex != secondSuggestionIndex {
					// makeRange is a helper function that will allow us to generate lists for intersection
					makeRange := func(min, max int) []int {
						a := make([]int, max-min+1)
						for i := range a {
							a[i] = min + i
						}
						return a
					}

					// The overlapping regions would need to have the same bias applied to them. So we check first if they have no preference or the same bias
					if (suggestion.Bias == "NA") || (secondSuggestion.Bias == "NA") || (suggestion.Bias == secondSuggestion.Bias) {
						// Each suggestion will have a start and end. If we take the intersection of the numbers between the start and the end of both suggestions,
						// we should be able to find a region where an edit would affect both suggested changes.
						overlap := intersect.Sorted(makeRange(suggestion.Start, suggestion.End), makeRange(secondSuggestion.Start, secondSuggestion.End)).([]int)
						if len(overlap) != 0 {
							// Each overlap that is successfully found will be inserted as a suggestedfix. First, we need to find its bias
							var overlapBias string
							if suggestion.Bias == "NA" {
								overlapBias = secondSuggestion.Bias
							} else {
								overlapBias = suggestion.Bias
							}

							// Each overlap will have a number of fixes that need to be hit to fix the outer objects
							var overlapQuantityFixes int
							if suggestion.QuantityFixes >= secondSuggestion.QuantityFixes {
								overlapQuantityFixes = suggestion.QuantityFixes
							} else {
								overlapQuantityFixes = secondSuggestion.QuantityFixes
							}

							// There can't be more desired fixes than there are positions
							if overlapQuantityFixes > len(overlap) {
								overlapQuantityFixes = len(overlap)
							}

							// Lets check if there is already an overlap at this position during this step
							var id int
							err := db.Get(&id, `SELECT id FROM suggestedfix WHERE step = ? AND start = ? AND end = ? AND suggestiontype = "overlap"`, i, overlap[0], overlap[len(overlap)-1])
							// If err is not nil, then we couldn't scan a single instance of the above query
							if err != nil {
								db.MustExec(`INSERT INTO suggestedfix(step, start, end, gcbias, quantityfixes, suggestiontype) VALUES (?, ?, ?, ?, ?, "overlap")`, i, overlap[0], overlap[len(overlap)-1], overlapBias, overlapQuantityFixes)
							}

						}
					}
				}
			}
		}

		// The following statements are the magic sauce that makes this all worthwhile.
		// Parameters: step, gcbias, start, end, quantityfix
		sqlFix := `INSERT INTO history (codon, pos, step) SELECT t.codon, t.pos, ? FROM (
				SELECT w.codon, w.pos 
				FROM seq AS s 
					JOIN history AS h ON h.pos = s.pos 
					JOIN weights AS w ON w.pos = s.pos 
					JOIN codon AS c ON h.codon = c.codon 
					JOIN codonbias AS cb ON cb.fromcodon = c.codon 
				WHERE cb.gbias = ?
					AND s.pos >= ?
					AND s.pos <= ?
					AND w.codon != h.codon 
					AND s.protected = false
				ORDER BY w.weight
			) AS t
			GROUP BY t.pos
			LIMIT ?;`

		// During this step, first see if there are any overlapping suggestedfixes
		var overlapSuggestions []DnaSuggestion
		db.Select(&overlapSuggestions, `SELECT * FROM suggestedfix WHERE suggestiontype = "overlap" AND step = ?`, i)
		// If there are overlapping suggestedfixes, fix those and iterate along
		if len(overlapSuggestions) > 0 {
			for _, overlapSuggestion := range overlapSuggestions {
				db.MustExec(sqlFix, i, overlapSuggestion.Bias, overlapSuggestion.Start, overlapSuggestion.End, overlapSuggestion.QuantityFixes)
			}
		} else { // If there aren't any overlapping suggestedfixes, fix the the remaining independent suggestedfixes
			var independentSuggestions []DnaSuggestion
			db.Select(&independentSuggestions, `SELECT * FROM suggestedfix WHERE step = ?`, i)
			for _, independentSuggestion := range independentSuggestions {
				db.MustExec(sqlFix, i, independentSuggestion.Bias, independentSuggestion.Start, independentSuggestion.End, independentSuggestion.QuantityFixes)
			}
		}
	}
	var codons []string
	db.Select(&codons, `SELECT codon FROM history GROUP BY pos ORDER BY (pos, step)`)
	finalSeq := strings.Join(codons, "")
	return finalSeq, nil
}
