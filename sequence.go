package poly

import (
	"bytes"
	"encoding/hex"
	"encoding/json"
	"strings"

	"lukechampine.com/blake3"
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

// AddFeature is the canonical way to add a Feature into a Sequence struct. Appending a Feature struct directly to Sequence.Feature's will break .GetSequence() method.
func (sequence *Sequence) AddFeature(feature Feature) {
	feature.parentSequence = sequence
	sequence.Features = append(sequence.Features, feature)
}

// RemoveFeature filters a given Feature from a sequences Features list.
func (sequence *Sequence) RemoveFeature(feature Feature) {

	var filteredFeatures []Feature
	for _, sequenceFeature := range sequence.Features {
		if !Equals(feature, sequenceFeature) {
			filteredFeatures = append(filteredFeatures, sequenceFeature)
		}
	}
	sequence.Features = filteredFeatures
}

// Length is a helper method to get Length of sequence.
func (sequence Sequence) Length() int {
	return len(sequence.Sequence)
}

type Equaler interface {
	Equals(*Feature) bool
}

// Equals compares two Features to determine if they're equal while ignoring unexported fields.
func Equals(feature1, feature2 Feature) bool {
	feature1Bytes, _ := json.Marshal(feature1)
	feature2Bytes, _ := json.Marshal(feature2)

	hashBytes1 := blake3.Sum256(feature1Bytes)
	hashBytes2 := blake3.Sum256(feature2Bytes)

	hashString1 := hex.EncodeToString(hashBytes1[:])
	hashString2 := hex.EncodeToString(hashBytes2[:])

	// this was seperated for debugging.
	hashEqual := hashString1 == hashString2

	return hashEqual

}

// InsertSequence inserts a subsequence into a sequence. This will ruin location coordinates unless features are updated.
func (sequence *Sequence) insertSequence(start int, insertSequence string) {

	var newSequence strings.Builder

	for index, basepair := range sequence.Sequence {
		if index == start {
			for _, insertBasepair := range insertSequence {
				newSequence.WriteRune(insertBasepair)
			}
		}
		newSequence.WriteRune(basepair)
	}

	sequence.Sequence = newSequence.String()
}

// deleteRange sequence method deletes a substring within a sequence. This will ruin location coordinates unless updated.
func (sequence *Sequence) deleteRange(start int, stop int) {

	// TODO handle corner case where stop is less than start

	var newSequence strings.Builder
	for index, basepair := range sequence.Sequence {
		// if index < start || index > stop {
		if !(index >= start && index < stop) {

			newSequence.WriteRune(basepair)
		}
	}

	sequence.Sequence = newSequence.String()
}

// GetRange is a method to get a subsequence over a certain range.
func (sequence Sequence) GetRange(start int, end int) string {
	var sequenceString strings.Builder

	for index, basepair := range sequence.Sequence {
		if index >= start && index < end {
			sequenceString.WriteRune(basepair)
		}
	}

	return sequenceString.String()
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

// getFeatureSequence takes a feature and location object and returns a sequence string.
func getFeatureSequence(feature Feature, location Location) string {
	var sequenceBuffer bytes.Buffer
	var sequenceString string
	parentSequence := feature.parentSequence.Sequence

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
