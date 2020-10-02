package main

import (
	"math/rand"
	"strings"
	"time"

	"github.com/mroth/weightedrand"
)

// CodonTranslationMap contains a one to one relation between codons and their corresponding amino acid.
// * represents stop codon
var CodonTranslationMap = map[string]string{
	"AAA": "K",
	"AAC": "N",
	"AAG": "K",
	"AAT": "N",
	"ACA": "T",
	"ACC": "T",
	"ACG": "T",
	"ACT": "T",
	"AGA": "R",
	"AGC": "S",
	"AGG": "R",
	"AGT": "S",
	"ATA": "I",
	"ATC": "I",
	"ATG": "M",
	"ATT": "I",
	"CAA": "Q",
	"CAC": "H",
	"CAG": "Q",
	"CAT": "H",
	"CCA": "P",
	"CCC": "P",
	"CCG": "P",
	"CCT": "P",
	"CGA": "R",
	"CGC": "R",
	"CGG": "R",
	"CGT": "R",
	"CTA": "L",
	"CTC": "L",
	"CTG": "L",
	"CTT": "L",
	"GAA": "E",
	"GAC": "D",
	"GAG": "E",
	"GAT": "D",
	"GCA": "A",
	"GCC": "A",
	"GCG": "A",
	"GCT": "A",
	"GGA": "G",
	"GGC": "G",
	"GGG": "G",
	"GGT": "G",
	"GTA": "V",
	"GTC": "V",
	"GTG": "V",
	"GTT": "V",
	"TAA": "*",
	"TAC": "Y",
	"TAG": "*",
	"TAT": "Y",
	"TCA": "S",
	"TCC": "S",
	"TCG": "S",
	"TCT": "S",
	"TGA": "*",
	"TGC": "C",
	"TGG": "W",
	"TGT": "C",
	"TTA": "L",
	"TTC": "F",
	"TTG": "L",
	"TTT": "F",
}

// AminoAcidMap contains a one to many mapping of AminoAcid single letter notations to their corresponding codons.
// * represents the stop codon
var AminoAcidMap = map[string]AminoAcid{
	"*": {Letter: "*", Codons: []string{"TAA", "TAG", "TGA"}}, //stop
	"A": {Letter: "A", Codons: []string{"GCA", "GCC", "GCG", "GCT"}},
	"C": {Letter: "C", Codons: []string{"TGC", "TGT"}},
	"D": {Letter: "D", Codons: []string{"GAC", "GAT"}},
	"E": {Letter: "E", Codons: []string{"GAA", "GAG"}},
	"F": {Letter: "F", Codons: []string{"TTC", "TTT"}},
	"G": {Letter: "G", Codons: []string{"GGA", "GGC", "GGG", "GGT"}},
	"H": {Letter: "H", Codons: []string{"CAC", "CAT"}},
	"I": {Letter: "I", Codons: []string{"ATA", "ATC", "ATT"}},
	"K": {Letter: "K", Codons: []string{"AAA", "AAG"}},
	"L": {Letter: "L", Codons: []string{"CTA", "CTC", "CTG", "CTT", "TTA", "TTG"}},
	"M": {Letter: "M", Codons: []string{"ATG"}}, // start
	"N": {Letter: "N", Codons: []string{"AAC", "AAT"}},
	"P": {Letter: "P", Codons: []string{"CCA", "CCC", "CCG", "CCT"}},
	"Q": {Letter: "Q", Codons: []string{"CAA", "CAG"}},
	"R": {Letter: "R", Codons: []string{"AGA", "AGG", "CGA", "CGC", "CGG", "CGT"}},
	"S": {Letter: "S", Codons: []string{"AGC", "AGT", "TCA", "TCC", "TCG", "TCT"}},
	"T": {Letter: "T", Codons: []string{"ACA", "ACC", "ACG", "ACT"}},
	"V": {Letter: "V", Codons: []string{"GTA", "GTC", "GTG", "GTT"}},
	"W": {Letter: "W", Codons: []string{"TGG"}},
	"Y": {Letter: "Y", Codons: []string{"TAC", "TAT"}},
}

// GetCodonFrequency takes a DNA sequence and returns a hashmap of its codons and their frequencies.
func GetCodonFrequency(sequence string) map[string]int {

	codonFrequencyHashMap := map[string]int{}
	var currentCodon strings.Builder

	for _, letter := range sequence {

		// add current nucleotide to currentCodon
		currentCodon.WriteRune(letter)

		// if current nucleotide is the third in a codon add to hashmap
		if currentCodon.Len() == 3 {
			// if codon is already initalized in map increment
			if _, ok := codonFrequencyHashMap[currentCodon.String()]; ok {
				codonFrequencyHashMap[currentCodon.String()]++
				// if codon is not already initalized in map initialize with value at 1
			} else {
				codonFrequencyHashMap[currentCodon.String()] = 1
			}
			// reset codon string builder for next codon.
			currentCodon.Reset()
		}
	}
	return codonFrequencyHashMap
}

func getCodonWeightMap(codonFrequencyHashMap map[string]int) map[string]weightedrand.Chooser {

	var codonWeightMap map[string]weightedrand.Chooser

	for _, aminoAcid := range AminoAcidMap {

		codonChoices := make([]weightedrand.Choice, len(aminoAcid.Codons))

		// get initial sum of aminoAcid in table
		var aminoAcidSum int
		for _, codon := range aminoAcid.Codons {
			aminoAcidSum += codonFrequencyHashMap[codon]
		}

		// construct codon choices and omit codons that occur less than 10%
		for _, codon := range aminoAcid.Codons {
			if codonFrequencyHashMap[codon]/aminoAcidSum >= 10.0 {
				codonChoices = append(codonChoices, weightedrand.Choice{Item: codon, Weight: uint(codonFrequencyHashMap[codon])})
			}
		}

		codonWeightMap[aminoAcid.Letter] = weightedrand.NewChooser(codonChoices...)

	}
	return codonWeightMap
}

// Translate translates a codon sequence to an amino acid sequence
func Translate(sequence string) string {
	var aminoAcids strings.Builder
	var currentCodon strings.Builder

	for _, letter := range sequence {

		// add current nucleotide to currentCodon
		currentCodon.WriteRune(letter)

		// if current nucleotide is the third in a codon translate to aminoAcid write to aminoAcids and reset currentCodon.
		if currentCodon.Len() == 3 {
			aminoAcids.WriteString(CodonTranslationMap[currentCodon.String()])

			// reset codon string builder for next codon.
			currentCodon.Reset()
		}
	}
	return aminoAcids.String()
}

// ReverseTranslate takes an amino acid string and a table of codons and returns a codon optimized nucleic acid string
func ReverseTranslate(aminoAcids string, randSeed int64, optimizationTable map[string]weightedrand.Chooser) string {

	var nucleotides strings.Builder

	// set seed makes testing of "randomness" possible for regression tests?
	rand.Seed(randSeed)

	// for each amino acid pick a random weighted codon and add to string builder
	for _, aminoAcid := range aminoAcids {
		nucleotides.WriteString(optimizationTable[string(aminoAcid)].Pick().(string))
	}

	return nucleotides.String()
}

func getOrganismOptimizationTable(annotatedSequence AnnotatedSequence) map[string]weightedrand.Chooser {

	// agglomerate all coding sequences into one big string buffer
	var sequenceBuffer strings.Builder
	for _, feature := range annotatedSequence.Features {
		if feature.Type == "CDS" {
			sequenceBuffer.WriteString(feature.getSequence())
		}
	}

	// get codon frequency hashmap to pass to construct weight map and construct weights then
	codonFrequencyHashMap := GetCodonFrequency(sequenceBuffer.String())
	optimizationTable := getCodonWeightMap(codonFrequencyHashMap)

	return optimizationTable
}

// CodonOptimizeForOrganism takes a sequence string and an AnnotatedSequence and return a codon optimized sequence for that AnnotatedSequence
func CodonOptimizeForOrganism(sequence string, annotatedSequence AnnotatedSequence) string {
	aminoAcids := Translate(sequence)
	optimizationTable := getOrganismOptimizationTable(annotatedSequence)
	randSeed := time.Now().UTC().UnixNano()

	optimizedSequence := ReverseTranslate(aminoAcids, randSeed, optimizationTable)

	return optimizedSequence
}

// Codon holds information for a codon triplet in a struct
type Codon struct {
	Triplet    string
	Occurrence int
}

// AminoAcid holds information for an amino acid and related codons in a struct in a struct
type AminoAcid struct {
	Letter string
	Codons []string
}

// // Function to generate default codon tables from NCBI https://www.ncbi.nlm.nih.gov/Taxonomy/Utils/wprintgc.cgi
// func generateCodonTable(aminoAcids, starts string) CodonTable {
// 	base1 := "TTTTTTTTTTTTTTTTCCCCCCCCCCCCCCCCAAAAAAAAAAAAAAAAGGGGGGGGGGGGGGGG"
// 	base2 := "TTTTCCCCAAAAGGGGTTTTCCCCAAAAGGGGTTTTCCCCAAAAGGGGTTTTCCCCAAAAGGGG"
// 	base3 := "TCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAG"
// 	// Add triplets to an amino acid -> triplet map, and if a possible start codon, add to start codon list
// 	var aminoAcidMap = make(map[rune][]Codon)
// 	var startCodons []string
// 	var stopCodons []string
// 	for i, aminoAcid := range aminoAcids {
// 		if _, ok := aminoAcidMap[aminoAcid]; ok == false {
// 			aminoAcidMap[aminoAcid] = []Codon{}
// 		}
// 		triplet := string([]byte{base1[i], base2[i], base3[i]})
// 		aminoAcidMap[aminoAcid] = append(aminoAcidMap[aminoAcid], Codon{triplet, 0})
// 		if starts[i] == 77 { // M rune
// 			startCodons = append(startCodons, triplet)
// 		}
// 		if starts[i] == 42 { // * rune
// 			stopCodons = append(stopCodons, triplet)
// 		}
// 	}
// 	// Convert amino acid -> triplet map to an amino acid list
// 	var aminoAcidSlice []AminoAcid
// 	for k, v := range aminoAcidMap {
// 		aminoAcidSlice = append(aminoAcidSlice, AminoAcid{string(k), v})
// 	}
// 	return CodonTable{startCodons, stopCodons, aminoAcidSlice}
// }

// // DefaultCodonTablesByNumber stores all codon tables published by NCBI https://www.ncbi.nlm.nih.gov/Taxonomy/Utils/wprintgc.cgi using numbered indeces.
// var DefaultCodonTablesByNumber = map[int]CodonTable{
// 	1:  generateCodonTable("FFLLSSSSYY**CC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG", "---M------**--*----M---------------M----------------------------"),
// 	2:  generateCodonTable("FFLLSSSSYY**CCWWLLLLPPPPHHQQRRRRIIMMTTTTNNKKSS**VVVVAAAADDEEGGGG", "----------**--------------------MMMM----------**---M------------"),
// 	3:  generateCodonTable("FFLLSSSSYY**CCWWTTTTPPPPHHQQRRRRIIMMTTTTNNKKSSRRVVVVAAAADDEEGGGG", "----------**----------------------MM---------------M------------"),
// 	4:  generateCodonTable("FFLLSSSSYY**CCWWLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG", "--MM------**-------M------------MMMM---------------M------------"),
// 	5:  generateCodonTable("FFLLSSSSYY**CCWWLLLLPPPPHHQQRRRRIIMMTTTTNNKKSSSSVVVVAAAADDEEGGGG", "---M------**--------------------MMMM---------------M------------"),
// 	6:  generateCodonTable("FFLLSSSSYYQQCC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG", "--------------*--------------------M----------------------------"),
// 	9:  generateCodonTable("FFLLSSSSYY**CCWWLLLLPPPPHHQQRRRRIIIMTTTTNNNKSSSSVVVVAAAADDEEGGGG", "----------**-----------------------M---------------M------------"),
// 	10: generateCodonTable("FFLLSSSSYY**CCCWLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG", "----------**-----------------------M----------------------------"),
// 	11: generateCodonTable("FFLLSSSSYY**CC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG", "---M------**--*----M------------MMMM---------------M------------"),
// 	12: generateCodonTable("FFLLSSSSYY**CC*WLLLSPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG", "----------**--*----M---------------M----------------------------"),
// 	13: generateCodonTable("FFLLSSSSYY**CCWWLLLLPPPPHHQQRRRRIIMMTTTTNNKKSSGGVVVVAAAADDEEGGGG", "---M------**----------------------MM---------------M------------"),
// 	14: generateCodonTable("FFLLSSSSYYY*CCWWLLLLPPPPHHQQRRRRIIIMTTTTNNNKSSSSVVVVAAAADDEEGGGG", "-----------*-----------------------M----------------------------"),
// 	16: generateCodonTable("FFLLSSSSYY*LCC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG", "----------*---*--------------------M----------------------------"),
// 	21: generateCodonTable("FFLLSSSSYY**CCWWLLLLPPPPHHQQRRRRIIMMTTTTNNNKSSSSVVVVAAAADDEEGGGG", "----------**-----------------------M---------------M------------"),
// 	22: generateCodonTable("FFLLSS*SYY*LCC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG", "------*---*---*--------------------M----------------------------"),
// 	23: generateCodonTable("FF*LSSSSYY**CC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG", "--*-------**--*-----------------M--M---------------M------------"),
// 	24: generateCodonTable("FFLLSSSSYY**CCWWLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSSKVVVVAAAADDEEGGGG", "---M------**-------M---------------M---------------M------------"),
// 	25: generateCodonTable("FFLLSSSSYY**CCGWLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG", "---M------**-----------------------M---------------M------------"),
// 	26: generateCodonTable("FFLLSSSSYY**CC*WLLLAPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG", "----------**--*----M---------------M----------------------------"),
// 	27: generateCodonTable("FFLLSSSSYYQQCCWWLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG", "--------------*--------------------M----------------------------"),
// 	28: generateCodonTable("FFLLSSSSYYQQCCWWLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG", "----------**--*--------------------M----------------------------"),
// 	29: generateCodonTable("FFLLSSSSYYYYCC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG", "--------------*--------------------M----------------------------"),
// 	30: generateCodonTable("FFLLSSSSYYEECC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG", "--------------*--------------------M----------------------------"),
// 	31: generateCodonTable("FFLLSSSSYYEECCWWLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG", "----------**-----------------------M----------------------------"),
// 	33: generateCodonTable("FFLLSSSSYYY*CCWWLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSSKVVVVAAAADDEEGGGG", "---M-------*-------M---------------M---------------M------------")}

// // DefaultCodonTablesByName stores all codon tables published by NCBI https://www.ncbi.nlm.nih.gov/Taxonomy/Utils/wprintgc.cgi using named indeces.
// var DefaultCodonTablesByName = map[string]CodonTable{
// 	"Standard":                         generateCodonTable("FFLLSSSSYY**CC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG", "---M------**--*----M---------------M----------------------------"), // 1
// 	"VertebrateMitochondrial":          generateCodonTable("FFLLSSSSYY**CCWWLLLLPPPPHHQQRRRRIIMMTTTTNNKKSS**VVVVAAAADDEEGGGG", "----------**--------------------MMMM----------**---M------------"), // 2
// 	"YeastMitochondrial":               generateCodonTable("FFLLSSSSYY**CCWWTTTTPPPPHHQQRRRRIIMMTTTTNNKKSSRRVVVVAAAADDEEGGGG", "----------**----------------------MM---------------M------------"), // 3
// 	"MoldMitochondrial":                generateCodonTable("FFLLSSSSYY**CCWWLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG", "--MM------**-------M------------MMMM---------------M------------"), // 4
// 	"ProtozoanMitochondrial":           generateCodonTable("FFLLSSSSYY**CCWWLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG", "--MM------**-------M------------MMMM---------------M------------"), // 4
// 	"CoelenterateMitochondrial":        generateCodonTable("FFLLSSSSYY**CCWWLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG", "--MM------**-------M------------MMMM---------------M------------"), // 4
// 	"Mycoplasma":                       generateCodonTable("FFLLSSSSYY**CCWWLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG", "--MM------**-------M------------MMMM---------------M------------"), // 4
// 	"Spiroplasma":                      generateCodonTable("FFLLSSSSYY**CCWWLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG", "--MM------**-------M------------MMMM---------------M------------"), // 4
// 	"InvertebrateMitochondrial":        generateCodonTable("FFLLSSSSYY**CCWWLLLLPPPPHHQQRRRRIIMMTTTTNNKKSSSSVVVVAAAADDEEGGGG", "---M------**--------------------MMMM---------------M------------"), // 5
// 	"CiliateNuclear":                   generateCodonTable("FFLLSSSSYYQQCC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG", "--------------*--------------------M----------------------------"), // 6
// 	"DasycladaceanNuclear":             generateCodonTable("FFLLSSSSYYQQCC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG", "--------------*--------------------M----------------------------"), // 6
// 	"HexamitaNuclear":                  generateCodonTable("FFLLSSSSYYQQCC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG", "--------------*--------------------M----------------------------"), // 6
// 	"EchinodermMitochondrial":          generateCodonTable("FFLLSSSSYY**CCWWLLLLPPPPHHQQRRRRIIIMTTTTNNNKSSSSVVVVAAAADDEEGGGG", "----------**-----------------------M---------------M------------"), // 9
// 	"FlatwormMitochondrial":            generateCodonTable("FFLLSSSSYY**CCWWLLLLPPPPHHQQRRRRIIIMTTTTNNNKSSSSVVVVAAAADDEEGGGG", "----------**-----------------------M---------------M------------"), // 9
// 	"EuplotidNuclear":                  generateCodonTable("FFLLSSSSYY**CCCWLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG", "----------**-----------------------M----------------------------"), // 10
// 	"BacterialPlastid":                 generateCodonTable("FFLLSSSSYY**CC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG", "---M------**--*----M------------MMMM---------------M------------"), // 11
// 	"ArchaelPlastid":                   generateCodonTable("FFLLSSSSYY**CC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG", "---M------**--*----M------------MMMM---------------M------------"), // 11
// 	"PlantPlastid":                     generateCodonTable("FFLLSSSSYY**CC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG", "---M------**--*----M------------MMMM---------------M------------"), // 11
// 	"AlternativeYeastNuclear":          generateCodonTable("FFLLSSSSYY**CC*WLLLSPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG", "----------**--*----M---------------M----------------------------"), // 12
// 	"AscidianMitochondrial":            generateCodonTable("FFLLSSSSYY**CCWWLLLLPPPPHHQQRRRRIIMMTTTTNNKKSSGGVVVVAAAADDEEGGGG", "---M------**----------------------MM---------------M------------"), // 13
// 	"AlternativeFlatwormMitochondrial": generateCodonTable("FFLLSSSSYYY*CCWWLLLLPPPPHHQQRRRRIIIMTTTTNNNKSSSSVVVVAAAADDEEGGGG", "-----------*-----------------------M----------------------------"), // 14
// 	"ChlorophyceanMitochondrial":       generateCodonTable("FFLLSSSSYY*LCC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG", "----------*---*--------------------M----------------------------"), // 16
// 	"TrematodeMitochondrial":           generateCodonTable("FFLLSSSSYY**CCWWLLLLPPPPHHQQRRRRIIMMTTTTNNNKSSSSVVVVAAAADDEEGGGG", "----------**-----------------------M---------------M------------"), // 21
// 	"ScenedesmusObliquusMitochondrial": generateCodonTable("FFLLSS*SYY*LCC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG", "------*---*---*--------------------M----------------------------"), // 22
// 	"ThraustochytriumMitochondrial":    generateCodonTable("FF*LSSSSYY**CC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG", "--*-------**--*-----------------M--M---------------M------------"), // 23
// 	"RhabdopleuridaeMitochondrial":     generateCodonTable("FFLLSSSSYY**CCWWLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSSKVVVVAAAADDEEGGGG", "---M------**-------M---------------M---------------M------------"), // 24
// 	"CandidateDivisionSR1":             generateCodonTable("FFLLSSSSYY**CCGWLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG", "---M------**-----------------------M---------------M------------"), // 25
// 	"Gracilibacteria":                  generateCodonTable("FFLLSSSSYY**CCGWLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG", "---M------**-----------------------M---------------M------------"), // 25
// 	"PachysolenTannophilusNuclear":     generateCodonTable("FFLLSSSSYY**CC*WLLLAPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG", "----------**--*----M---------------M----------------------------"), // 26
// 	"KaryorelictNuclear":               generateCodonTable("FFLLSSSSYYQQCCWWLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG", "--------------*--------------------M----------------------------"), // 27
// 	"CondylostomaNuclear":              generateCodonTable("FFLLSSSSYYQQCCWWLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG", "----------**--*--------------------M----------------------------"), // 28
// 	"MesodiniumNuclear":                generateCodonTable("FFLLSSSSYYYYCC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG", "--------------*--------------------M----------------------------"), // 29
// 	"PeritrichNuclear":                 generateCodonTable("FFLLSSSSYYEECC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG", "--------------*--------------------M----------------------------"), // 30
// 	"BlastocrithidiaNuclear":           generateCodonTable("FFLLSSSSYYEECCWWLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG", "----------**-----------------------M----------------------------"), // 31
// 	"CephalodiscidaeMitochondrial":     generateCodonTable("FFLLSSSSYYY*CCWWLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSSKVVVVAAAADDEEGGGG", "---M-------*-------M---------------M---------------M------------")} // 33

// func main() {
// 	rand.Seed(time.Now().UTC().UnixNano())
// 	codonTable := CodonTable{[]string{}, []string{}, []AminoAcid{{"M", []Codon{{"ATG", 1}}}, {"G", []Codon{{"GGA", 1}}}, {"*", []Codon{{"TAA", 1}, {"TGA", 1}}}}}
// 	translation := Translate("ATGGGCTGA", DefaultCodonTablesByName["PeritrichNuclear"])

// 	fmt.Println(Optimize(translation, codonTable))
// }
