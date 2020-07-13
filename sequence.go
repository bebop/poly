package main

import (
	"strconv"
	"strings"
)

// getSequence is a method to get a Feature's sequence. Mutates with AnnotatedSequence.
func (f Feature) getSequence(as AnnotatedSequence) string {
	return getFeatureSequence(f.Location, as.Sequence.Sequence)
}

// ReverseComplement takes the reverse complement of a sequence
func ReverseComplement(sequence string) string {
	complementString := strings.Map(complementBase, sequence)
	n := len(complementString)
	newString := make([]rune, n)
	for _, base := range complementString {
		n--
		newString[n] = base
	}
	return string(newString)
}

// complementBase accepts a base pair and returns its complement base pair
func complementBase(basePair rune) rune {
	complementMap := make(map[rune]rune)
	bases := "ATUGCYRSWKMBDHVN"
	complementBases := "TAACGRYSWMKVHDBN"
	lowerComplementBases := strings.ToLower(complementBases)

	for i, base := range bases {
		complementMap[base] = rune(complementBases[i])
	}

	for i, base := range strings.ToLower(bases) {
		complementMap[base] = rune(lowerComplementBases[i])
	}

	return complementMap[basePair]
}

// getFeatureSequence is a recursive function to get the sequence of a feature, given that feature's Location and parent full sequence
func getFeatureSequence(s string, fullSeq string) string {
	if !(strings.ContainsAny(s, "(")) { // Case checks for simple expression of x..x
		startEndList := strings.Split(s, "..")
		start, _ := strconv.Atoi(startEndList[0])
		end, _ := strconv.Atoi(startEndList[1])
		return fullSeq[start-1 : end] // Sequences are indexed at 1, but inclusive of last element
	} else {
		firstOuterParathesis := strings.Index(s, "(")
		expression := s[firstOuterParathesis+1 : strings.LastIndex(s, ")")]
		switch command := s[0:firstOuterParathesis]; command {
		case "join":
			// This case checks for join(complement(x..x),complement(x..x)), or any more complicated derivatives
			if strings.ContainsAny(expression, "(") {
				firstInnerParathesis := strings.Index(expression, "(")
				parathesisCount := 1
				comma := 0
				for i := 1; parathesisCount > 0; i++ { // "(" is at 0, so we start at 1
					comma = i
					switch expression[firstInnerParathesis+i] {
					case []byte("(")[0]:
						parathesisCount++
					case []byte(")")[0]:
						parathesisCount--
					}
				}
				return getFeatureSequence(expression[:firstInnerParathesis+comma+1], fullSeq) + getFeatureSequence(expression[2+firstInnerParathesis+comma:], fullSeq)
			} else { // This is the default join(x..x,x..x)
				joinedSequence := ""
				for _, numberRange := range strings.Split(expression, ",") {
					joinedSequence += getFeatureSequence(numberRange, fullSeq)
				}
				return joinedSequence
			}

		case "complement":
			return ReverseComplement(getFeatureSequence(expression, fullSeq))
		}
	}
	return "failed"
}
