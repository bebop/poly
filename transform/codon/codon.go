package codon

import (
	"errors"

	"math/rand"
	"strings"
	"time"

	"github.com/TimothyStiles/poly"
	weightedRand "github.com/mroth/weightedrand"

	"encoding/json"
	"io/ioutil"
)

/******************************************************************************
Oct, 15, 2020

File is structured as so:

	Structs:
		Table - holds all information mapping codons <-> aminoacids during transformations.
		AnimoAcid - holds amino acide related info for Table struct
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

Anywho, most of this file and Table's struct methods are meant to help overcome
this codon bias. There's a default Table generator near the bottom of this file
with a whole section on how it works and why it's gotta be that way.

Like most codebases, best usage examples come from tests. You can check out
TestTranslate and TestOptimize in transformations_test.go for pretty solid
examples of what you can do with this code.

TTFN,
Tim

******************************************************************************/

var errEmtpyCodonTable error = errors.New("empty codon table")
var errEmtpyAminoAcidString error = errors.New("empty amino acid string")
var errEmtpySequenceString error = errors.New("empty sequence string")

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

// Table holds information for a codon table.
type Table struct {
	StartCodons []string    `json:"start_codons"`
	StopCodons  []string    `json:"stop_codons"`
	AminoAcids  []AminoAcid `json:"amino_acids"`
}

// Translate translates a codon sequence to an amino acid sequence
func Translate(sequence string, codonTable Table) (string, error) {
	if len(codonTable.StartCodons) == 0 && len(codonTable.StopCodons) == 0 && len(codonTable.AminoAcids) == 0 {
		return "", errEmtpyCodonTable
	}
	if len(sequence) == 0 {
		return "", errEmtpySequenceString
	}

	var aminoAcids strings.Builder
	var currentCodon strings.Builder
	translationTable := codonTable.generateTranslationTable()

	for _, letter := range sequence {

		// add current nucleotide to currentCodon
		currentCodon.WriteRune(letter)

		// if current nucleotide is the third in a codon translate to aminoAcid write to aminoAcids and reset currentCodon.
		if currentCodon.Len() == 3 {
			aminoAcids.WriteString(translationTable[strings.ToUpper(currentCodon.String())])

			// reset codon string builder for next codon.
			currentCodon.Reset()
		}
	}
	return aminoAcids.String(), nil
}

// Optimize takes an amino acid sequence and Table and returns an optimized codon sequence
func Optimize(aminoAcids string, codonTable Table) (string, error) {
	if len(codonTable.StartCodons) == 0 && len(codonTable.StopCodons) == 0 && len(codonTable.AminoAcids) == 0 {
		return "", errEmtpyCodonTable
	}
	if len(aminoAcids) == 0 {
		return "", errEmtpyAminoAcidString
	}

	// weightedRand library insisted setting seed like this. Not sure what environmental side effects exist.
	rand.Seed(time.Now().UTC().UnixNano())

	var codons strings.Builder
	codonChooser := codonTable.chooser()

	for _, aminoAcid := range aminoAcids {
		aminoAcidString := string(aminoAcid)
		codons.WriteString(codonChooser[aminoAcidString].Pick().(string))
	}
	return codons.String(), nil
}

// OptimizeTable weights each codon in a codon table according to input string codon frequency.
// This function actually mutates the Table struct itself.
func (codonTable Table) OptimizeTable(sequence string) Table {

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

// GetCodingRegions is a helper function to pull coding regions out of an Sequence as input for optimizing codon tables.
func GetCodingRegions(sequence poly.Sequence) string {
	// pick out the each coding region in the Sequence and add it to the sequence Builder
	var sequenceBuilder strings.Builder

	for _, feature := range sequence.Features {
		if feature.Type == "CDS" {
			sequenceBuilder.WriteString(feature.GetSequence())
		}
	}

	return sequenceBuilder.String()
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

// chooser is a Table method to convert a codon table to a chooser
func (codonTable Table) chooser() map[string]weightedRand.Chooser {

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
func (codonTable Table) generateTranslationTable() map[string]string {
	var translationMap = make(map[string]string)
	for _, aminoAcid := range codonTable.AminoAcids {
		for _, codon := range aminoAcid.Codons {
			translationMap[codon.Triplet] = aminoAcid.Letter
		}
	}
	return translationMap
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
I made one that went by name too. Looking back on it this was useless so I removed
it.

Happy hacking,
Tim

******************************************************************************/

// Function to generate default codon tables from NCBI https://www.ncbi.nlm.nih.gov/Taxonomy/Utils/wprintgc.cgi
func generateCodonTable(aminoAcids, starts string) Table {
	base1 := "TTTTTTTTTTTTTTTTCCCCCCCCCCCCCCCCAAAAAAAAAAAAAAAAGGGGGGGGGGGGGGGG"
	base2 := "TTTTCCCCAAAAGGGGTTTTCCCCAAAAGGGGTTTTCCCCAAAAGGGGTTTTCCCCAAAAGGGG"
	base3 := "TCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAG"
	// Add triplets to an amino acid -> triplet map, and if a possible start codon, add to start codon list
	var aminoAcidMap = make(map[rune][]Codon)
	var startCodons []string
	var stopCodons []string
	for i, aminoAcid := range aminoAcids {
		if _, ok := aminoAcidMap[aminoAcid]; !ok {
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
	return Table{startCodons, stopCodons, aminoAcidSlice}
}

// GetCodonTable takes the index of desired NCBI codon table and returns it.
func GetCodonTable(index int) Table {
	return defaultCodonTablesByNumber[index]
}

// defaultCodonTablesByNumber stores all codon tables published by NCBI https://www.ncbi.nlm.nih.gov/Taxonomy/Utils/wprintgc.cgi using numbered indeces.
var defaultCodonTablesByNumber = map[int]Table{
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

******************************************************************************/

// ParseCodonJSON parses a Table JSON file.
func ParseCodonJSON(file []byte) Table {
	var codontable Table
	_ = json.Unmarshal([]byte(file), &codontable)
	return codontable
}

// ReadCodonJSON reads a Table JSON file.
func ReadCodonJSON(path string) Table {
	file, _ := ioutil.ReadFile(path)
	codontable := ParseCodonJSON(file)
	return codontable
}

// WriteCodonJSON writes a Table struct out to JSON.
func WriteCodonJSON(codontable Table, path string) {
	file, _ := json.MarshalIndent(codontable, "", " ")
	_ = ioutil.WriteFile(path, file, 0644)
}

/******************************************************************************
Dec, 17, 2020

Compromise + Add codon table stuff begins here


== Compromise tables ==
Basically, I want to codon optimize a protein for two or more organisms.
In order to do that, I need to be able to generate a codon table that
is a compromise between the codon tables between two different organisms.

The method is fairly simple: standardize codon counts so the weights are
equal between both organisms, then add them together. In addition, have
a variable percentage for removing rare codons (this makes compromise
tables lossy).

Simple code, but very powerful if it can be used to encode genes for
multiple organisms.

== Add tables ==
Some organisms have multiple chromosomes. We need to add em all up
to get an accurate codon table (different than compromise tables,
since these are all already balanced).

Godspeed,

Keoni
******************************************************************************/

// CompromiseCodonTable takes 2 CodonTables and makes a new Table
// that is an equal compromise between the two tables.
func CompromiseCodonTable(firstCodonTable Table, secondCodonTable Table, cutOff float64) (Table, error) {
	// Initialize output Table, c
	var c Table
	// Check if cutOff is too high or low (this is converted to a percent)
	if cutOff < 0 {
		return c, errors.New("Cut off too low. Cannot be less than 0 or greater than 1")
	}
	if cutOff > 1 {
		return c, errors.New("Cut off too high. Cannot be greater than 1")
	}

	// Take start and stop strings from first table
	// and use them as start + stops in final Table
	c.StartCodons = firstCodonTable.StartCodons
	c.StopCodons = firstCodonTable.StopCodons

	// Initialize the finalAminoAcid list for the output Table
	var finalAminoAcids []AminoAcid

	// Loop over all AminoAcids represented in the first Table
	for _, firstAa := range firstCodonTable.AminoAcids {
		var firstTriplets []string
		var firstWeights []int
		var firstTotal int

		var secondWeights []int
		var secondTotal int
		// For each amino acid in firstCodonTable, get list of all codons, and append triplets
		// and weights to a list
		for _, firstCodon := range firstAa.Codons {
			firstTriplets = append(firstTriplets, firstCodon.Triplet)
			firstWeights = append(firstWeights, firstCodon.Weight)
			firstTotal = firstTotal + firstCodon.Weight
			for _, secondAa := range secondCodonTable.AminoAcids {
				if secondAa.Letter == firstAa.Letter {
					for _, secondCodon := range secondAa.Codons {
						// For each codon from firstCodonTable, get the
						// corresponding triplet and weight from secondCodonTable
						if secondCodon.Triplet == firstCodon.Triplet {
							secondWeights = append(secondWeights, secondCodon.Weight)
							secondTotal = secondTotal + secondCodon.Weight
						}
					}
				}
			}
		}

		var finalTriplets []string
		var finalWeights []int
		cutOffWeight := int(10000 * cutOff)

		// For each of the Triplets in the amino acid, output a triplet weight
		// for the first and second triplet, which is the percentage of Triplets
		// coding for that amino acid multiplied by 10,000
		for i, firstTriplet := range firstTriplets {
			finalTriplets = append(finalTriplets, firstTriplet)
			firstTripletWeight := int((float64(firstWeights[i]) / float64(firstTotal)) * 10000)
			secondTripletWeight := int((float64(secondWeights[i]) / float64(secondTotal)) * 10000)
			// If the triplet is less than the cutoff weight in either the first or second table,
			// set its weight to zero. Otherwise, append the average of the first and second weight
			// to final weights
			if (firstTripletWeight < cutOffWeight) || (secondTripletWeight < cutOffWeight) {
				finalWeights = append(finalWeights, 0)
			} else {
				finalWeights = append(finalWeights, int((float64(firstTripletWeight)+float64(secondTripletWeight))/2))
			}
		}
		// From those final weights and final triplets, build a list of Codons
		var finalCodons []Codon
		for i, finalTriplet := range finalTriplets {
			finalCodons = append(finalCodons, Codon{finalTriplet, finalWeights[i]})
		}

		// Append list of Codons to finalAminoAcids
		finalAminoAcids = append(finalAminoAcids, AminoAcid{firstAa.Letter, finalCodons})
	}
	c.AminoAcids = finalAminoAcids
	return c, nil
}

// AddCodonTable takes 2 CodonTables and adds them together to create
// a new Table.
func AddCodonTable(firstCodonTable Table, secondCodonTable Table) Table {
	var c Table

	// Take start and stop strings from first table
	c.StartCodons = firstCodonTable.StartCodons
	c.StopCodons = firstCodonTable.StopCodons

	// Add up codons
	var finalAminoAcids []AminoAcid
	for _, firstAa := range firstCodonTable.AminoAcids {
		var finalCodons []Codon
		for _, firstCodon := range firstAa.Codons {
			for _, secondAa := range secondCodonTable.AminoAcids {
				for _, secondCodon := range secondAa.Codons {
					if firstCodon.Triplet == secondCodon.Triplet {
						finalCodons = append(finalCodons, Codon{firstCodon.Triplet, firstCodon.Weight + secondCodon.Weight})
					}
				}
			}
		}
		finalAminoAcids = append(finalAminoAcids, AminoAcid{firstAa.Letter, finalCodons})
	}
	c.AminoAcids = finalAminoAcids
	return c
}
