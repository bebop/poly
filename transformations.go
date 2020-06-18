package main

import (
	wr "github.com/mroth/weightedrand"
	"math/rand"
	"time"
)

type codon struct {
	triplet    string
	occurrence int
}

type aminoAcid struct {
	letter string
	codons []codon
}

type codonTable struct {
	startCodons []string
	stopCodons []string
	aminoAcid   []aminoAcid
}

// Generate map of amino acid -> codon chooser
func (ct codonTable) optimizationTable() map[string]wr.Chooser {
	rand.Seed(time.Now().UTC().UnixNano())
	var optMap = make(map[string]wr.Chooser)
	for _, aa := range ct.aminoAcid {
		// Get list of triplets and their weights
		codonChoices := make([]wr.Choice, len(aa.codons))
		totals := 0
		for _, codon := range aa.codons {
			codonChoices = append(codonChoices, wr.Choice{Item: codon.triplet, Weight: uint(codon.occurrence)})
			totals += codon.occurrence
		}
		optMap[aa.letter] = wr.NewChooser(codonChoices...)
	}
	return optMap
}

//Convert from amino acid sequence -> codon sequence
func Optimize(aminoAcids string, ct codonTable) string {
	var codons string
	optTable := ct.optimizationTable()
	for _, aa := range aminoAcids {
		codons += optTable[string(aa)].Pick().(string)
	}
	return codons
}

// Generate map of codons -> amino acid
func (ct codonTable) translationTable() map[string]string {
	var transMap = make(map[string]string)
	for _, aa := range ct.aminoAcid {
		for _, codon := range aa.codons {
			transMap[codon.triplet] = aa.letter
		}
	}
	return transMap
}

// Convert from codon sequence -> amino acid sequence
func Translate(nucleotides string, ct codonTable) string {
	var aminoAcids string
	transTable := ct.translationTable()
	l := len(nucleotides)/3 // Assumes input sequences are divisble by 3
	for i := 0; l>i; i++ {
		aminoAcids += transTable[nucleotides[i*3:(i+1)*3]]
	}
	return aminoAcids
}

// Function to generate default codon tables from NCBI https://www.ncbi.nlm.nih.gov/Taxonomy/Utils/wprintgc.cgi
func genTable(aas, starts string) codonTable {
	base1 := "TTTTTTTTTTTTTTTTCCCCCCCCCCCCCCCCAAAAAAAAAAAAAAAAGGGGGGGGGGGGGGGG"
	base2 := "TTTTCCCCAAAAGGGGTTTTCCCCAAAAGGGGTTTTCCCCAAAAGGGGTTTTCCCCAAAAGGGG"
	base3 := "TCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAGTCAG"
	// Add triplets to an amino acid -> triplet map, and if a possible start codon, add to start codon list
	var aaMap = make(map[rune][]codon)
	var startCodons []string
	var stopCodons []string
	for i, aa := range aas {
		if _, ok := aaMap[aa]; ok == false {
			aaMap[aa] = []codon{}
		}
		triplet := string([]byte{base1[i],base2[i],base3[i]})
		aaMap[aa] = append(aaMap[aa], codon{triplet, 0})
		if starts[i] == 77 { // M rune
			startCodons = append(startCodons, triplet)
		}
		if starts[i] == 42 { // * rune
			stopCodons = append(stopCodons, triplet)
		}
	}
	// Convert amino acid -> triplet map to an amino acid list
	var aaList []aminoAcid
	for k,v := range aaMap {
		aaList = append(aaList, aminoAcid{string(k),v})
	}
	return codonTable{startCodons,stopCodons,aaList}
}

// Public variable for all codon tables published by NCBI https://www.ncbi.nlm.nih.gov/Taxonomy/Utils/wprintgc.cgi
var DefaultCodonTables = map[int]codonTable{
	1: genTable("FFLLSSSSYY**CC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG","---M------**--*----M---------------M----------------------------"),
	2: genTable("FFLLSSSSYY**CCWWLLLLPPPPHHQQRRRRIIMMTTTTNNKKSS**VVVVAAAADDEEGGGG","----------**--------------------MMMM----------**---M------------"),
	3: genTable("FFLLSSSSYY**CCWWTTTTPPPPHHQQRRRRIIMMTTTTNNKKSSRRVVVVAAAADDEEGGGG","----------**----------------------MM---------------M------------"),
	4: genTable("FFLLSSSSYY**CCWWLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG","--MM------**-------M------------MMMM---------------M------------"),
	5: genTable("FFLLSSSSYY**CCWWLLLLPPPPHHQQRRRRIIMMTTTTNNKKSSSSVVVVAAAADDEEGGGG","---M------**--------------------MMMM---------------M------------"),
	6: genTable("FFLLSSSSYYQQCC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG","--------------*--------------------M----------------------------"),
	9: genTable("FFLLSSSSYY**CCWWLLLLPPPPHHQQRRRRIIIMTTTTNNNKSSSSVVVVAAAADDEEGGGG","----------**-----------------------M---------------M------------"),
	10: genTable("FFLLSSSSYY**CCCWLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG","----------**-----------------------M----------------------------"),
	11: genTable("FFLLSSSSYY**CC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG","---M------**--*----M------------MMMM---------------M------------"),
	12: genTable("FFLLSSSSYY**CC*WLLLSPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG","----------**--*----M---------------M----------------------------"),
	13: genTable("FFLLSSSSYY**CCWWLLLLPPPPHHQQRRRRIIMMTTTTNNKKSSGGVVVVAAAADDEEGGGG","---M------**----------------------MM---------------M------------"),
	14: genTable("FFLLSSSSYYY*CCWWLLLLPPPPHHQQRRRRIIIMTTTTNNNKSSSSVVVVAAAADDEEGGGG","-----------*-----------------------M----------------------------"),
	16: genTable("FFLLSSSSYY*LCC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG","----------*---*--------------------M----------------------------"),
	21: genTable("FFLLSSSSYY**CCWWLLLLPPPPHHQQRRRRIIMMTTTTNNNKSSSSVVVVAAAADDEEGGGG","----------**-----------------------M---------------M------------"),
	22: genTable("FFLLSS*SYY*LCC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG","------*---*---*--------------------M----------------------------"),
	23: genTable("FF*LSSSSYY**CC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG","--*-------**--*-----------------M--M---------------M------------"),
	24: genTable("FFLLSSSSYY**CCWWLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSSKVVVVAAAADDEEGGGG","---M------**-------M---------------M---------------M------------"),
	25: genTable("FFLLSSSSYY**CCGWLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG","---M------**-----------------------M---------------M------------"),
	26: genTable("FFLLSSSSYY**CC*WLLLAPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG","----------**--*----M---------------M----------------------------"),
	27: genTable("FFLLSSSSYYQQCCWWLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG","--------------*--------------------M----------------------------"),
	28: genTable("FFLLSSSSYYQQCCWWLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG","----------**--*--------------------M----------------------------"),
	29: genTable("FFLLSSSSYYYYCC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG","--------------*--------------------M----------------------------"),
	30: genTable("FFLLSSSSYYEECC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG","--------------*--------------------M----------------------------"),
	31: genTable("FFLLSSSSYYEECCWWLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG","----------**-----------------------M----------------------------"),
	33: genTable("FFLLSSSSYYY*CCWWLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSSKVVVVAAAADDEEGGGG","---M-------*-------M---------------M---------------M------------")}

//func main() {
//	rand.Seed(time.Now().UTC().UnixNano())
//	ct := codonTable{[]string{},[]string{},[]aminoAcid{{"M",[]codon{{"ATG", 1}}},{"G",[]codon{{"GGA", 1}}},{"*",[]codon{{"TAA",1},{"TGA",1}}}}}
//	trans := translate("ATGGGCTGA", defaultCodonTables[1])
//	fmt.Println(optimize(trans,ct))
//	fmt.Println(optimize(trans,ct))
//	fmt.Println(optimize(trans,ct))
//	fmt.Println(optimize(trans,ct))
//}
