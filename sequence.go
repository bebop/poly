package main

import (
	"bytes"
	"strconv"
	"strings"
)

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

// getSequence is a method to get a Feature's sequence. Mutates with AnnotatedSequence.
func (feature Feature) getSequence() string {
	return getFeatureSequence(feature)
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
	return ComplementBaseRuneMap[basePair]
}

func parseGbkLocation(locationString string) Location {
	var location Location
	if !(strings.ContainsAny(locationString, "(")) { // Case checks for simple expression of x..x
		startEndSplit := strings.Split(locationString, "..")
		start, _ := strconv.Atoi(startEndSplit[0])
		end, _ := strconv.Atoi(startEndSplit[1])
		location = Location{Start: start - 1, End: end}
	} else {
		firstOuterParentheses := strings.Index(locationString, "(")
		expression := locationString[firstOuterParentheses+1 : strings.LastIndex(locationString, ")")]
		switch command := locationString[0:firstOuterParentheses]; command {
		case "join":
			// This case checks for join(complement(x..x),complement(x..x)), or any more complicated derivatives
			if strings.ContainsAny(expression, "(") {
				firstInnerParentheses := strings.Index(expression, "(")
				ParenthesesCount := 1
				comma := 0
				for i := 1; ParenthesesCount > 0; i++ { // "(" is at 0, so we start at 1
					comma = i
					switch expression[firstInnerParentheses+i] {
					case []byte("(")[0]:
						ParenthesesCount++
					case []byte(")")[0]:
						ParenthesesCount--
					}
				}
				location.SubLocations = append(location.SubLocations, parseGbkLocation(expression[:firstInnerParentheses+comma+1]), parseGbkLocation(expression[2+firstInnerParentheses+comma:]))
			} else { // This is the default join(x..x,x..x)
				for _, numberRange := range strings.Split(expression, ",") {
					location.SubLocations = append(location.SubLocations, parseGbkLocation(numberRange))
				}
			}

		case "complement":
			subLocation := parseGbkLocation(expression)
			subLocation.Complement = true
			location.SubLocations = append(location.SubLocations, subLocation)
		}
	}
	return location
}

func getFeatureSequenceHelper(feature Feature, location Location) string {
	var sequenceBuffer bytes.Buffer
	var sequenceString string
	parentSequence := feature.ParentAnnotatedSequence.Sequence.Sequence
	if len(location.SubLocations) == 0 {
		sequenceBuffer.WriteString(parentSequence[location.Start:location.End])
	} else {
		for _, subLocation := range location.SubLocations {
			sequenceBuffer.WriteString(getFeatureSequenceHelper(feature, subLocation))
		}
	}

	if location.Complement {
		sequenceString = ReverseComplement(sequenceBuffer.String())
	} else {
		sequenceString = sequenceBuffer.String()
	}

	return sequenceString
}

func getFeatureSequence(feature Feature) string {
	location := parseGbkLocation(feature.GbkLocationString)

	return getFeatureSequenceHelper(feature, location)
}
