package poly

import (
	"math/rand"
	"strings"
	"time"

	weightedRand "github.com/mroth/weightedrand"
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
	Triplet string
	Weight  int // needs to be set to 1 for random chooser
}

// AminoAcid holds information for an amino acid and related codons in a struct in a struct
type AminoAcid struct {
	Letter string
	Codons []Codon
}

// CodonTable holds information for a codon table.
type CodonTable struct {
	StartCodons []string
	StopCodons  []string
	AminoAcids  []AminoAcid
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

// CreateWeights weights each codon in a codon table according to input string codon frequency
// this function actually mutates the CodonTable struct itself and doesn't technically need to
// return but is useful if you want to create a copy while mutating... Not sure if deep copy is
// appropriate
func (codonTable CodonTable) CreateWeights(sequence string) CodonTable {

	sequence = strings.ToUpper(sequence)
	codonFrequencyMap := getCodonFrequency(sequence)

	for aminoAcidIndex, aminoAcid := range codonTable.AminoAcids {

		// Get sum of codon occurences for particular amino acid
		codonOccurenceSum := 0
		for _, codon := range aminoAcid.Codons {
			triplet := codonFrequencyMap[codon.Triplet]
			codonOccurenceSum += triplet
		}

		// Threshold codons that occur less than 10% for coding a particular amino acid
		for _, codon := range aminoAcid.Codons {
			codonPercentage := float64(codonFrequencyMap[codon.Triplet]) / float64(codonOccurenceSum)

			if codonPercentage < 0.10 {
				codonFrequencyMap[codon.Triplet] = 0
			}

		}

		// apply thresholded weights to codonTable
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

		// for every codon related to current amino acid append its Triplet and Weight to codonChoices
		for _, codon := range aminoAcid.Codons {
			codonChoices = append(codonChoices, weightedRand.Choice{Item: codon.Triplet, Weight: uint(codon.Weight)})
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
