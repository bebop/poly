package poly

import (
	"bytes"
	"errors"
	"math/rand"
	"strings"
)

// complementBaseRuneMap provides 1:1 mapping between bases and their complements
var complementBaseRuneMap = map[rune]rune{
	65:  84,  // A -> T
	66:  86,  // B -> V
	67:  71,  // C -> G
	68:  72,  // D -> H
	71:  67,  // G -> C
	72:  68,  // H -> D
	75:  77,  // K -> M
	77:  75,  // M -> K
	78:  78,  // N -> N
	82:  89,  // R -> Y
	83:  83,  // S -> S
	84:  65,  // T -> A
	85:  65,  // U -> A
	86:  66,  // V -> B
	87:  87,  // W -> W
	89:  82,  // Y -> R
	97:  116, // a -> t
	98:  118, // b -> v
	99:  103, // c -> g
	100: 104, // d -> h
	103: 99,  // g -> a
	104: 100, // h -> d
	107: 109, // k -> m
	109: 107, // m -> k
	110: 110, // n -> n
	114: 121, // r -> y
	115: 115, // s -> s
	116: 97,  // t -> a
	117: 97,  // u -> a
	118: 98,  // v -> b
	119: 119, // w -> w
	121: 114, // y -> r
}

// GetSequence is a method to get the full sequence of an annotated sequence
func (sequence Sequence) GetSequence() string {
	return sequence.Sequence
}

// GetSequence is a method wrapper to get a Feature's sequence. Mutates with Sequence.
func (feature Feature) GetSequence() string {
	return getFeatureSequence(feature, feature.SequenceLocation)
}

// ReverseComplement takes the reverse complement of a sequence
func ReverseComplement(sequence string) string {
	complementString := strings.Map(ComplementBase, sequence)
	n := len(complementString)
	newString := make([]rune, n)
	for _, base := range complementString {
		n--
		newString[n] = base
	}
	return string(newString)
}

// ComplementBase accepts a base pair and returns its complement base pair
func ComplementBase(basePair rune) rune {
	return complementBaseRuneMap[basePair]
}

//RandomProteinSequence returns a random protein sequence as a string that have size length, starts with aminoacid M (Methionine) and finishes with * (stop codon). The random generator uses the seed provided as parameter.
func RandomProteinSequence(length int, seed int64) (string, error) {
	//The length needs to be greater than two because the random protein sequenced returned always contain a start and stop codon. You could see more about this stuff here: https://en.wikipedia.org/wiki/Genetic_code#Start_and_stop_codons
	if length <= 2 {
		err := errors.New("The length needs to be greater than two because the random protein sequenced returned always contain a start and stop codon. Please select a higher length in RandomProteinSequence function")
		return "", err
	}

	// https://en.wikipedia.org/wiki/Amino_acid#Table_of_standard_amino_acid_abbreviations_and_properties
	var aminoAcidsAlphabet = []rune("ACDEFGHIJLMNPQRSTVWY")
	rand.Seed(seed)

	randomSequence := make([]rune, length)

	for peptide := range randomSequence {
		if peptide == 0 {
			//M is the standard abbreviation for the Methionine aminoacid. A protein sequence start with M because the start codon is translated to Methionine
			randomSequence[peptide] = 'M'
		} else if peptide == length-1 {
			//* is the standard abbreviation for the stop codon. That's a signal for the ribosome to stop the translation and because of that a protein sequence is finished with *
			randomSequence[peptide] = '*'
		} else {
			randomIndex := rand.Intn(len(aminoAcidsAlphabet))
			randomSequence[peptide] = aminoAcidsAlphabet[randomIndex]
		}
	}

	return string(randomSequence), nil
}

// AllVariantsIUPAC takes a string as input
// and returns all iupac variants as output
func AllVariantsIUPAC(seq string) []string {
	seqVariantList := [][]rune{}
	seqVariants := []string{}

	iupac := map[rune][]rune{ // rune map of all iupac nucleotide variants
		'G': []rune{'G'},
		'A': []rune{'A'},
		'T': []rune{'T'},
		'C': []rune{'C'},
		'R': []rune{'G', 'A'},
		'Y': []rune{'T', 'C'},
		'M': []rune{'A', 'C'},
		'K': []rune{'G', 'T'},
		'S': []rune{'G', 'C'},
		'W': []rune{'A', 'T'},
		'H': []rune{'A', 'C', 'T'},
		'B': []rune{'G', 'T', 'C'},
		'V': []rune{'G', 'C', 'A'},
		'D': []rune{'G', 'A', 'T'},
		'N': []rune{'G', 'A', 'T', 'C'},
	}

	for _, s := range seq {
		seqVariantList = append(seqVariantList, iupac[s])
	}

	cartesianProducts := cartRune(seqVariantList...)
	for _, product := range cartesianProducts {
		seqVariants = append(seqVariants, string(product))
	}
	return seqVariants
}

func cartRune(inList ...[]rune) [][]rune {
	// An iteratitive approach to calculate Cartesian product of two or more lists
	// Adapted from https://rosettacode.org/wiki/Cartesian_product_of_two_or_more_lists
	// supposedly "minimizes allocations and computes and fills the result sequentially"

	var possibleVariants int = 1 // a counter used to determine the possible number of variants
	for _, inList := range inList {
		possibleVariants *= len(inList)
	}
	if possibleVariants == 0 {
		return nil // in the future this could be part of an error return?
	}
	allVariants := make([][]rune, possibleVariants)              // this is the 2D slice where all variants will be stored
	variantHolders := make([]rune, possibleVariants*len(inList)) // this is an empty slice with a length totaling the size of all input characters
	variantChoices := make([]int, len(inList))                   // these will be all the possible variants
	start := 0
	for variant := range allVariants {
		end := start + len(inList) // define end point
		variantHolder := variantHolders[start:end]

		allVariants[variant] = variantHolder

		start = end // start at end point

		for variantChoicesIndex, variantChoice := range variantChoices {
			variantHolder[variantChoicesIndex] = inList[variantChoicesIndex][variantChoice]
		}
		for variantChoicesIndex := len(variantChoices) - 1; variantChoicesIndex >= 0; variantChoicesIndex-- {
			variantChoices[variantChoicesIndex]++
			if variantChoices[variantChoicesIndex] < len(inList[variantChoicesIndex]) {
				break
			}
			variantChoices[variantChoicesIndex] = 0
		}
	}
	return allVariants
}

// getFeatureSequence takes a feature and location object and returns a sequence string.
func getFeatureSequence(feature Feature, location Location) string {
	var sequenceBuffer bytes.Buffer
	var sequenceString string
	parentSequence := feature.ParentSequence.Sequence

	if len(location.SubLocations) == 0 {
		sequenceBuffer.WriteString(parentSequence[location.Start:location.End])
	} else {

		for _, subLocation := range location.SubLocations {
			sequenceBuffer.WriteString(getFeatureSequence(feature, subLocation))
		}
	}

	// reverse complements resulting string if needed.
	if location.Complement {
		sequenceString = ReverseComplement(sequenceBuffer.String())
	} else {
		sequenceString = sequenceBuffer.String()
	}

	return sequenceString
}
