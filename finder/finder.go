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

func MatchSequences(sequences map[string]string) func(string) []Match {
	return func(sequence string) []Match {
		var enzymes []string
		var matches []Match
		for pattern, comment := range sequences {
			enzymes = []string{pattern, transform.ReverseComplement(pattern)}
			for _, site := range enzymes {
				re := regexp.MustCompile(site)
				locs := re.FindAllStringIndex(sequence, -1)
				for _, loc := range locs {
					matches = append(matches, Match{loc[0], loc[1], comment + " | " + site})
				}
			}
		}
		return matches
	}
}

// RemoveRepeat is a generator to make a problematicSequenceFunc for repeats.
func RemoveRepeat(repeatLen int) func(string) []Match {
	return func(sequence string) []Match {
		// Get a kmer list
		var matches []Match
		kmers := make(map[string]bool)
		for i := 0; i < len(sequence)-repeatLen; i++ {
			_, alreadyFound := kmers[sequence[i:i+repeatLen]]
			if alreadyFound {
				matches = append(matches, Match{i, i + repeatLen, "Repeated sequence | " + sequence[i:i+repeatLen]})
			}
			kmers[sequence[i:i+repeatLen]] = true
		}
		return matches
	}
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
					matches = append(matches, Match{loc[0], loc[1], "Forbidden sequence | " + site})
				}
			}
		}
		return matches
	}
}

func GlobalRemoveRepeat(repeatLen int, globalSequence string) func(string) []Match {
	globalKmers := transform.GetKmerTable(repeatLen, globalSequence)
	return func(sequence string) []Match {
		var matches []Match
		for i := 0; i < len(sequence)-repeatLen; i++ {
			subsequence := sequence[i : i+repeatLen]
			reverseComplement := transform.ReverseComplement(subsequence)
			_, forwardGlobalRepeat := globalKmers[subsequence]
			_, reverseGlobalRepeat := globalKmers[reverseComplement]
			if forwardGlobalRepeat {
				matches = append(matches, Match{i, i + repeatLen, "Global repeated sequence | " + subsequence})
			}
			if reverseGlobalRepeat {
				matches = append(matches, Match{i, i + repeatLen, "Global repeated sequence | " + reverseComplement})
			}
		}
		return matches
	}
}
