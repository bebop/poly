package main

import (
	"math/rand"
	"time"

	weightedRand "github.com/mroth/weightedrand"
)

// Codon holds information for a codon triplet in a struct
type Codon struct {
	Triplet    string
	Occurrence int
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

// Generate map of amino acid -> codon chooser
func (codonTable CodonTable) generateOptimizationTable() map[string]weightedRand.Chooser {
	rand.Seed(time.Now().UTC().UnixNano())
	var optimizationMap = make(map[string]weightedRand.Chooser)
	for _, aminoAcid := range codonTable.AminoAcids {
		// Get list of triplets and their weights
		codonChoices := make([]weightedRand.Choice, len(aminoAcid.Codons))
		totals := 0
		for _, codon := range aminoAcid.Codons {
			codonChoices = append(codonChoices, weightedRand.Choice{Item: codon.Triplet, Weight: uint(codon.Occurrence)})
			totals += codon.Occurrence
		}
		optimizationMap[aminoAcid.Letter] = weightedRand.NewChooser(codonChoices...)
	}
	return optimizationMap
}

// Optimize takes an amino acid sequence and CodonTable and returns an optimized codon sequence
func Optimize(aminoAcids string, codonTable CodonTable) string {
	var codons string
	optimizationTable := codonTable.generateOptimizationTable()
	for _, aminoAcid := range aminoAcids {
		codons += optimizationTable[string(aminoAcid)].Pick().(string)
	}
	return codons
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

// Translate translates a codon sequence to an amino acid sequence
func Translate(nucleotides string, codonTable CodonTable) string {
	var aminoAcids string
	translationTable := codonTable.generateTranslationTable()
	length := len(nucleotides) / 3 // Assumes input sequences are divisible by 3
	for i := 0; length > i; i++ {
		aminoAcids += translationTable[nucleotides[i*3:(i+1)*3]]
	}
	return aminoAcids
}

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
		aminoAcidMap[aminoAcid] = append(aminoAcidMap[aminoAcid], Codon{triplet, 0})
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

// DefaultCodonTables stores all codon tables published by NCBI https://www.ncbi.nlm.nih.gov/Taxonomy/Utils/wprintgc.cgi
var DefaultCodonTables = map[int]CodonTable{
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

// func main() {
// 	rand.Seed(time.Now().UTC().UnixNano())
// 	codonTable := CodonTable{[]string{}, []string{}, []AminoAcid{{"M", []Codon{{"ATG", 1}}}, {"G", []Codon{{"GGA", 1}}}, {"*", []Codon{{"TAA", 1}, {"TGA", 1}}}}}
// 	translation := Translate("ATGGGCTGA", DefaultCodonTables[1])

// 	fmt.Println(Optimize(translation, codonTable))
// }
