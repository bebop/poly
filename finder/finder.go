package finder

import (
	"regexp"

	"github.com/Open-Science-Global/poly"
	"github.com/Open-Science-Global/poly/transform"
)

type Match struct {
	Start   int
	End     int
	Message string
}

func Find(sequence string, finderFunctions []func(string) []Match) []Match {
	var matches []Match
	for _, finderFunction := range finderFunctions {
		matches = append(matches, finderFunction(sequence)...)
	}
	return matches
}

func AddMatchesToSequence(matches []Match, sequence poly.Sequence) poly.Sequence {
	for _, match := range matches {
		var feature poly.Feature
		attributes := make(map[string]string)
		attributes["label"] = match.Message
		feature.Attributes = attributes
		feature.Type = "Match"
		feature.SequenceLocation.Start = match.Start
		feature.SequenceLocation.End = match.End
		sequence.AddFeature(&feature)
	}
	return sequence
}

// RemoveSequence is a generator for a problematicSequenceFuncs for specific sequences.
func ForbiddenSequence(sequencesToRemove []string) func(string) []Match {
	return func(sequence string) []Match {
		var enzymes []string
		var matches []Match
		for _, enzyme := range sequencesToRemove {
			enzymes = []string{enzyme, transform.ReverseComplement(enzyme)}
			for _, site := range enzymes {
				re := regexp.MustCompile(site)
				locs := re.FindAllStringIndex(sequence, -1)
				for _, loc := range locs {
					matches = append(matches, Match{loc[0], loc[1], "Forbidden sequence: " + site})
				}
			}
		}
		return matches
	}
}
