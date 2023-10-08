/*
Package random provides functions to generate random DNA, RNA, and protein sequences.
*/
package random

import (
	"errors"
	"math/rand"

	"github.com/TimothyStiles/poly/transform/variants"
)

var iupacToBases = variants.IUPAC2Bases()
var iupacToAAs = variants.IUPAC2AAs()

// ProteinSequence returns a random protein sequence string of a given length and seed.
// All returned sequences start M (Methionine) and end with * (stop codon).
func ProteinSequence(length int, seed int64) (string, error) {
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

func RandomRune(runes []rune) rune {
	randomIndex := rand.Intn(len(runes))
	return runes[randomIndex]
}

func ProteinSequenceFromPattern(pattern string, seed int64) string {
	randomSequence := make([]rune, len(pattern))
	rand.Seed(seed)

	for res, peptide := range pattern {
		if options, found := iupacToAAs[peptide]; found {
			randomSequence[res] = RandomRune(options)
		} else {
			randomSequence[res] = peptide
		}
	}
	return string(randomSequence)
}

// DNASequence returns a random DNA sequence string of a given length and seed.
func DNASequence(length int, seed int64) (string, error) {
	return randomNucelotideSequence(length, seed, []rune("ACTG")), nil
}

func DNASequenceFromPattern(pattern string, seed int64) string {
	randomSequence := make([]rune, len(pattern))
	rand.Seed(seed)

	for i, base := range pattern {
		if options, found := iupacToBases[base]; found {
			randomSequence[i] = RandomRune(options)
		} else {
			randomSequence[i] = base
		}
	}
	return string(randomSequence)
}

// RNASequence returns a random DNA sequence string of a given length and seed.
func RNASequence(length int, seed int64) (string, error) {
	return randomNucelotideSequence(length, seed, []rune("ACUG")), nil
}

func randomNucelotideSequence(length int, seed int64, alphabet []rune) string {
	alphabetLength := len(alphabet)
	rand.Seed(seed)

	randomSequence := make([]rune, length)
	for basepair := range randomSequence {
		randomIndex := rand.Intn(alphabetLength)
		randomSequence[basepair] = alphabet[randomIndex]
	}

	return string(randomSequence)
}
