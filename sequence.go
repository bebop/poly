package main

import (
	"strconv"
	"strings"
)

// getSequence is a method to get a Feature's sequence. Mutates with AnnotatedSequence.
func (feature Feature) getSequence() string {
	return getFeatureSequence(feature.Location, feature.ParentAnnotatedSequence.Sequence.Sequence)
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

// ComplementBaseRuneMap provides 1:1 mapping between bases and their complements
var ComplementBaseRuneMap = map[rune]rune{
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

// complementBase accepts a base pair and returns its complement base pair
func complementBase(basePair rune) rune {
	return ComplementBaseRuneMap[basePair]
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
