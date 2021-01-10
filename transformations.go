package poly

import (
	"math/rand"
	"strings"
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

******************************************************************************/

type DnaFix struct {
	Start: int,
	End: int,
	Bias: string,
	InternalFix: []DnaFix
}

type Triplet struct {
	Letter: string
	Position: int
	Choices: []weightedRand.Choice
}

func FindBsaI(searchSeq, c chan DnaFix, wg *sync.WaitGroup) {
	defer wg.Done()

	// This should be done with a regex to find all occurences, but this is just a test
	// function so.... yea.
	i := searchSeq(searchSeq, "GGTCTC")
	c <- DnaFix{i, i+6, "NA"}
}

func FixCodons(sequence string, codontable CodonTable, problematicSequenceFuncs []func(string, chan DnaFix, *sync.WaitGroup)) string, error{

	// Check that sequence is triplets
	// {}

	// Uppercase the sequence
	// {}

	// Translate sequence
	translation := Translate(sequence, codontable)

	// Generate a list of triplets that represent the entire sequence
	var triplets []Triplet
	for i, aa := range translation {
		for _, aminoAcid := range codontable.AminoAcids {
			if aminoAcid.Letter == aa {

				// Get sum of codon occurences for a particular amino acid
				codonOccurenceSum := 0
				for _, codon in aminoAcid.Codons {
					codonOccurenceSum += codon.Weight
				}

				// Threshold codons that occur less than 10% for coding a particular amino acid
				for _, codon := range aminoAcids {
					codonPercentage := float64(codon.Weight) / float64(codonOccurenceSum)

					var codonChoices []weightedRand.Choice
					if codonPercentage > 0.10 {
						// for every codon related to current amino acid append its Triplet and Weight to codonChoices after thresholding
						// except for the codon currently at the translation location
						if codon.Triplet != sequence[i*3:(i*3)+3] {
							codonChoices = append(codonChoices, weightedRand.Choice{Item: codon.Triplet, Weight: uint(codon.Weight)})
						}
					}
				}
				triplets = append(triplets, Triplet{aa, i, codonChoices})
			}
		}
	}

	// Build initial fix list
	var fixes chan DnaFix
	var wg sync.WaitGroup
	for _, f := range problematicSequenceFuncs {
		wg.Add(1)
		go f(sequence, fixes, wg)
	}
	wg.Wait()

	iterationCount := 0
	// Loop until all fixes are fixed or there has been 1000 iterations
	for (len(fixes) != 0) || (iterationCount < 1000) {
		for _, fix := range fixes {
			// bestFixScore keeps score of the best externalFix (ie, the smallest one that
			// still covers the entire sequence)
			var bestFixScore int
			var bestExternalFix DnaFix
			fittingExternal := false
			for _, externalFix := range fixes {
				if (fix.Start > externalFix.Start) && (fix.End < externalFix.End) {
					externalFixScore := (fix.Start - externalFix.Start) + (externalFix.End - fix.End)
					if externalFixScore < bestFixScore {
						bestExternalFix = externalFix
						fittingExternal = true
					}
				}
			}
			// If there is a best externalFix, add the current fix to its fix list
			if fittingExternal == true {
				externalFix.InternalFix = append(externalFix.InternalFix, fix)
			}
		}
		// Only go actually fix things without any internal fixes
		for _, fix := range fixes {
			if fix.InternalFix == []DnaFix {
				// It is ok to round down since wobbles are usually on right side
				tripletStart = fix.Start / 3
				tripletEnd = fix.End / 3
				var fixableTriplets []Triplet
				for i := tripletStart; i < tripletEnd; i++ {
					fixableTriplets = append(fixableTriplets, triplets[i])
				}

				// This entire switch is code to Bias towards GC or AT content
				tripletChoosers := make(map[int][]weightedRand.Choice)
				switch {
				case fix.Bias == "GC":
					for i, triplet := fixableTriplets {
						var finalChoices []weightedRand.Choice
						for _, choice := range triplet.Choices {
							if (strings.Count(choice.Item, "G") + strings.Count(choice.Item, "C")) > (strings.Count(sequence[Triplet.Position * 3:(Triplet.Position* 3) + 3], "G") + strings.Count(sequence[Triplet.Position * 3:(Triplet.Position * 3) + 3], "C")) {
								finalChoices = append(finalChoices, choice)
							}
						}
						tripletChoosers[i] = finalChoices
					}
				case fix.Bias == "AT":
					for i, triplet := fixableTriplets {
                                                var finalChoices []weightedRand.Choice
                                                for _, choice := range triplet.Choices {
                                                        if (strings.Count(choice.Item, "A") + strings.Count(choice.Item, "T")) > (strings.Count(sequence[Triplet.Position * 3:(Triplet.Position * 3) + 3], "A") + strings.Count(sequence[Triplet.Position * 3:(Triplet.Position * 3) + 3], "T")) {
                                                                finalChoices = append(finalChoices, choice)
                                                        }
                                                }
						tripletChoosers[i] = finalChoices
                                        }
				case fix.Bias == "NA":
					for i, triplet := fixableTriplets {
						var finalChoices []weightedRand.Choice
						for _, choice := range triplet.Choices {
							finalChoices = append(finalChoices, choice)
						}
					}
					tripletChoosers[i] = finalChoices
				}
				// Now that we have our tripletChoosers to yoink from, we have to choose one to use
				targetPosition = int(len(fixableTriplets) / 2)
				var middleOutTripletList []Triplet
				var frontTriplets []Triplet
				var endTriplets []Triplet
				for i, t := range fixableTriplets {
					switch {
					case i < targetPosition:
						frontTriplets = append(frontTriplets, t)
					case i > targetPosition:
						endTriplets = append(endTriplets, t)
					case i == targetPosition:
						middleOutTripletList = append(middleOutTripletList, t)
					}
				}

				// Now sort em in middle out.
				// For example, we have list [1, 2, 3, 4]
				// we will reorder to list [2, 3, 1, 4]
				if len(frontTriplets) < len(endTriplets) {
					oddTriplets := endTriplets
					evenTriplets := frontTriplets
				} else {
					oddTriplets := frontTriplets
					evenTriplets := endTriplets
				}
				for i := 1; i < len(frontTriplets) + len(endTriplets); i++ {
					if(n%2==0) {
						middleOutTripletList = append(middleOutTripletList, endTriplets[(i/2)-1])
					} else {
						middleOutTripletList = append(middleOutTripletList, frontTriplets[((i+1)/2)-1])
					}
				}

				// For each Triplet in middleOutTripletList, check if there is any option to change to. If there is,
				// change it, and then note that in the original Triplet list. 
				for i, tr := range middleOutTripletList {
					if len(tripletChoosers[tr.Position]) == 0 {
						if i+1 == len(middleOutTripletList) {
							return sequence, errors.New("FAILED - NO TRIPLET TO CHOOSE FROM")
						}
					} else {
						change := weightedRand.Chooser(tripletChoosers[tr.Position]...).Pick().(string)
						changePosition := tr.Position
						break
					}
				}

				// Now we have the change we want to make. Let's make the change to the sequence itself, and then
				//



					}
				}
			}
		}
	}
}



type DnaFix2 struct {
	Positions: []int
	Bias: string

}

func FixSequence2(sequence string, codontable CodonTable, problematicSequenceFuncs []func(string, chan DnaFix *sync.WaitGroup)) (string, error) {
	aaSequence = Translation(sequence, codontable)
	// Initialize a SQLite database representing the sequence
	db = database, _ := sql.Open("sqlite3", ":memory:")

	// Initialize Tables
	statement, _ := db.Prepare("CREATE TABLE aminoacid(aminoacid TEXT PRIMARY KEY)")
	statement.Exec()
	statement, _ = db.Prepare("CREATE TABLE codonchoices(id INTEGER PRIMARY KEY, tripletchoice TEXT, weight INT, aminoacid TEXT, FOREIGN KEY(aminoacid) REFERENCES aminoacid(aminoacid))")
	statement.Exec()
	statement, _ = db.Prepare("CREATE TABLE position (position INTEGER PRIMARY KEY, aminoacid TEXT, FOREIGN KEY(aminoacid) REFERENCES aminoacid(aminoacid))")
	statement.Exec()
	statement, _ = db.Prepare("CREATE TABLE codon (id INTEGER PRIMARY KEY, position INTEGER, triplet TEXT, level INTEGER, FOREIGN KEY(position) REFERENCES position(position))")
	statement.Exec()

	// Fill aminoacid with codonTable data
	statement, _ = db.Prepare("INSERT INTO aminoacid(aminoacid) VALUES (?)")
	for _, aa := range codontable.AminoAcids {
		statement.Exec(aa.Letter)
	}

	// Fill codonchoices with codonTable data
	statement, _ = db.Prepare("INSERT INTO codonchoices(tripletchoice, weight, aminoacid) VALUES (?,?,?)")
	for _, aa := range codontable.AminoAcids {
		for _, codon := range aa.Codons {
			statement.Exec(codon.Triplet, codon.Weight, aa.Letter)
		}
	}

	// Fill position with translation data
	statement, _ = db.Prepare("INSERT INTO position (position, aminoacid) VALUES (?,?)")
	for i, aa := range aaSequence {
		statement.Exec(i, aa)
	}

	// Fill codon with sequence specific data
	statement, _ = db.Prepare("INSERT INTO codon (position, triplet, level) VALUES (?,?,?)")
	for i, aa := range aaSequence {
		statement.Exec(i, sequence[i*3:(i*3)+3], 0)
	}

	func getErrors() {
		// Get full sequence to pass to problematicSequenceFuncs
	}
	return "hi", nil

}



