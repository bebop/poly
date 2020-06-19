package main

import (
	"math/rand"
	"time"

	wr "github.com/mroth/weightedrand"
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
func (codonTable CodonTable) optimizationTable() map[string]wr.Chooser {
	rand.Seed(time.Now().UTC().UnixNano())
	var optMap = make(map[string]wr.Chooser)
	for _, aa := range codonTable.AminoAcids {
		// Get list of triplets and their weights
		codonChoices := make([]wr.Choice, len(aa.Codons))
		totals := 0
		for _, codon := range aa.Codons {
			codonChoices = append(codonChoices, wr.Choice{Item: codon.Triplet, Weight: uint(codon.Occurrence)})
			totals += codon.Occurrence
		}
		optMap[aa.Letter] = wr.NewChooser(codonChoices...)
	}
	return optMap
}

// Optimize takes an amino acid sequence and CodonTable and returns an optimized codon sequence
func Optimize(aminoAcids string, codonTable CodonTable) string {
	var codons string
	optTable := codonTable.optimizationTable()
	for _, aa := range aminoAcids {
		codons += optTable[string(aa)].Pick().(string)
	}
	return codons
}

// Generate map of codons -> amino acid
func (codonTable CodonTable) translationTable() map[string]string {
	var transMap = make(map[string]string)
	for _, aa := range codonTable.AminoAcids {
		for _, codon := range aa.Codons {
			transMap[codon.Triplet] = aa.Letter
		}
	}
	return transMap
}

// Translate translates a codon sequence to an amino acid sequence
func Translate(nucleotides string, codonTable CodonTable) string {
	var aminoAcids string
	transTable := codonTable.translationTable()
	l := len(nucleotides) / 3 // Assumes input sequences are divisble by 3
	for i := 0; l > i; i++ {
		aminoAcids += transTable[nucleotides[i*3:(i+1)*3]]
	}
	return aminoAcids
}

// Function to generate default codon tables from NCBI https://www.ncbi.nlm.nih.gov/Taxonomy/Utils/wprintgc.cgi
func genTable(aas, starts string) CodonTable {
	base1 := "TTTTTTTTTTTTTTTTCCCCCCCCCCCCCCCCAAAAAAAAAAAAAAAAGGGGGGGGGGGGGGGG"
	base2 := "TTTTCCCCAAAAGGGGTTTTCCCCAAAAGGGGTTTTCCCCAAAAGGGGTTTTCCCCAAAAGGGG"
	base3 := "TCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAG"
	// Add triplets to an amino acid -> triplet map, and if a possible start codon, add to start codon list
	var aaMap = make(map[rune][]Codon)
	var startCodons []string
	var stopCodons []string
	for i, aa := range aas {
		if _, ok := aaMap[aa]; ok == false {
			aaMap[aa] = []Codon{}
		}
		triplet := string([]byte{base1[i], base2[i], base3[i]})
		aaMap[aa] = append(aaMap[aa], Codon{triplet, 0})
		if starts[i] == 77 { // M rune
			startCodons = append(startCodons, triplet)
		}
		if starts[i] == 42 { // * rune
			stopCodons = append(stopCodons, triplet)
		}
	}
	// Convert amino acid -> triplet map to an amino acid list
	var aaList []AminoAcid
	for k, v := range aaMap {
		aaList = append(aaList, AminoAcid{string(k), v})
	}
	return CodonTable{startCodons, stopCodons, aaList}
}

// DefaultCodonTables stores all codon tables published by NCBI https://www.ncbi.nlm.nih.gov/Taxonomy/Utils/wprintgc.cgi
var DefaultCodonTables = map[int]CodonTable{
	1:  genTable("FFLLSSSSYY**CC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG", "---M------**--*----M---------------M----------------------------"),
	2:  genTable("FFLLSSSSYY**CCWWLLLLPPPPHHQQRRRRIIMMTTTTNNKKSS**VVVVAAAADDEEGGGG", "----------**--------------------MMMM----------**---M------------"),
	3:  genTable("FFLLSSSSYY**CCWWTTTTPPPPHHQQRRRRIIMMTTTTNNKKSSRRVVVVAAAADDEEGGGG", "----------**----------------------MM---------------M------------"),
	4:  genTable("FFLLSSSSYY**CCWWLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG", "--MM------**-------M------------MMMM---------------M------------"),
	5:  genTable("FFLLSSSSYY**CCWWLLLLPPPPHHQQRRRRIIMMTTTTNNKKSSSSVVVVAAAADDEEGGGG", "---M------**--------------------MMMM---------------M------------"),
	6:  genTable("FFLLSSSSYYQQCC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG", "--------------*--------------------M----------------------------"),
	9:  genTable("FFLLSSSSYY**CCWWLLLLPPPPHHQQRRRRIIIMTTTTNNNKSSSSVVVVAAAADDEEGGGG", "----------**-----------------------M---------------M------------"),
	10: genTable("FFLLSSSSYY**CCCWLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG", "----------**-----------------------M----------------------------"),
	11: genTable("FFLLSSSSYY**CC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG", "---M------**--*----M------------MMMM---------------M------------"),
	12: genTable("FFLLSSSSYY**CC*WLLLSPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG", "----------**--*----M---------------M----------------------------"),
	13: genTable("FFLLSSSSYY**CCWWLLLLPPPPHHQQRRRRIIMMTTTTNNKKSSGGVVVVAAAADDEEGGGG", "---M------**----------------------MM---------------M------------"),
	14: genTable("FFLLSSSSYYY*CCWWLLLLPPPPHHQQRRRRIIIMTTTTNNNKSSSSVVVVAAAADDEEGGGG", "-----------*-----------------------M----------------------------"),
	16: genTable("FFLLSSSSYY*LCC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG", "----------*---*--------------------M----------------------------"),
	21: genTable("FFLLSSSSYY**CCWWLLLLPPPPHHQQRRRRIIMMTTTTNNNKSSSSVVVVAAAADDEEGGGG", "----------**-----------------------M---------------M------------"),
	22: genTable("FFLLSS*SYY*LCC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG", "------*---*---*--------------------M----------------------------"),
	23: genTable("FF*LSSSSYY**CC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG", "--*-------**--*-----------------M--M---------------M------------"),
	24: genTable("FFLLSSSSYY**CCWWLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSSKVVVVAAAADDEEGGGG", "---M------**-------M---------------M---------------M------------"),
	25: genTable("FFLLSSSSYY**CCGWLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG", "---M------**-----------------------M---------------M------------"),
	26: genTable("FFLLSSSSYY**CC*WLLLAPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG", "----------**--*----M---------------M----------------------------"),
	27: genTable("FFLLSSSSYYQQCCWWLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG", "--------------*--------------------M----------------------------"),
	28: genTable("FFLLSSSSYYQQCCWWLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG", "----------**--*--------------------M----------------------------"),
	29: genTable("FFLLSSSSYYYYCC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG", "--------------*--------------------M----------------------------"),
	30: genTable("FFLLSSSSYYEECC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG", "--------------*--------------------M----------------------------"),
	31: genTable("FFLLSSSSYYEECCWWLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG", "----------**-----------------------M----------------------------"),
	33: genTable("FFLLSSSSYYY*CCWWLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSSKVVVVAAAADDEEGGGG", "---M-------*-------M---------------M---------------M------------")}

//func main() {
//	rand.Seed(time.Now().UTC().UnixNano())
//	codonTable := CodonTable{[]string{},[]string{},[]aminoAcid{{"M",[]Codon{{"ATG", 1}}},{"G",[]Codon{{"GGA", 1}}},{"*",[]Codon{{"TAA",1},{"TGA",1}}}}}
//	trans := translate("ATGGGCTGA", defaultCodonTables[1])
//	fmt.Println(optimize(trans,codonTable))
//	fmt.Println(optimize(trans,codonTable))
//	fmt.Println(optimize(trans,codonTable))
//	fmt.Println(optimize(trans,codonTable))
//}
