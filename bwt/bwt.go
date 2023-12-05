package bwt

import (
	"strings"

	"golang.org/x/exp/slices"
)

const nullChar = "$"

// BWT Burrow Wheeler Transform
// Data structure that compactly represents any sequence of characters and
// allows for sub sequence querying.
type BWT struct {
	skipList []skipEntry
	l        waveletTree
}

func (bwt BWT) Count(pattern string) int {
	skip, ok := bwt.lookupSkip(pattern[len(pattern)-1])
	if !ok {
		return 0
	}
	nextRange := skip.openEndedInterval
	for i := 1; i < len(pattern); i++ {
		if nextRange.end-nextRange.start <= 0 {
			return 0
		}

		currChar := pattern[len(pattern)-1-i]

		currCharRangeStart := bwt.l.Rank(currChar, nextRange.start)
		currCharRangeEnd := bwt.l.Rank(currChar, nextRange.end)

		nextCharSkip, ok := bwt.lookupSkip(currChar)
		if !ok {
			return 0
		}

		nextRange.start = nextCharSkip.openEndedInterval.start + currCharRangeStart
		nextRange.end = nextCharSkip.openEndedInterval.start + currCharRangeEnd
	}
	return nextRange.end - nextRange.start
}

func (bwt BWT) lookupSkip(c byte) (entry skipEntry, ok bool) {
	for i := range bwt.skipList {
		if bwt.skipList[i].char == c {
			return bwt.skipList[i], true
		}
	}
	return skipEntry{}, false
}

type interval struct {
	start int
	end   int
}

type skipEntry struct {
	char              byte
	openEndedInterval interval
}

func New(sequence string) BWT {
	sequence += nullChar

	prefixArray := make([]string, len(sequence))
	for i := 0; i < len(sequence); i++ {
		prefixArray[i] = sequence[len(sequence)-i-1:]
	}

	slices.Sort(prefixArray)

	lastColBuilder := strings.Builder{}
	for i := 0; i < len(prefixArray); i++ {
		currChar := sequence[getBWTIndex(len(sequence), len(prefixArray[i]))]
		lastColBuilder.WriteByte(currChar)
	}

	return BWT{
		skipList: buildSkipList(prefixArray),
		l:        NewWaveletTreeFromString(lastColBuilder.String()),
	}
}

func buildSkipList(prefixArray []string) []skipEntry {
	prevChar := prefixArray[0][0]
	skipList := []skipEntry{{char: prevChar, openEndedInterval: interval{start: 0}}}
	for i := 1; i < len(prefixArray); i++ {
		currChar := prefixArray[i][0]
		if currChar != prevChar {
			skipList[len(skipList)-1].openEndedInterval.end = i
			skipList = append(skipList, skipEntry{
				char:              currChar,
				openEndedInterval: interval{start: i},
			})
			prevChar = currChar
		}
	}
	skipList[len(skipList)-1].openEndedInterval.end = len(prefixArray)
	return skipList
}

func getBWTIndex(lenOfSequenceBeingBuilt, lenOfSuffixArrayVisited int) int {
	bwtCharIndex := lenOfSequenceBeingBuilt - lenOfSuffixArrayVisited - 1
	if bwtCharIndex == -1 {
		bwtCharIndex = lenOfSequenceBeingBuilt - 1
	}
	return bwtCharIndex
}
