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
	searchRange := interval{start: 0, end: bwt.getLenOfOriginalString()}
	for i := 0; i < len(pattern); i++ {
		if searchRange.end-searchRange.start <= 0 {
			return 0
		}

		c := pattern[len(pattern)-1-i]
		skip, ok := bwt.lookupSkip(c)
		if !ok {
			return 0
		}
		searchRange.start = skip.openEndedInterval.start + bwt.l.Rank(c, searchRange.start)
		searchRange.end = skip.openEndedInterval.start + bwt.l.Rank(c, searchRange.end)
	}
	return searchRange.end - searchRange.start
}

func (bwt BWT) lookupSkip(c byte) (entry skipEntry, ok bool) {
	for i := range bwt.skipList {
		if bwt.skipList[i].char == c {
			return bwt.skipList[i], true
		}
	}
	return skipEntry{}, false
}

func (bwt BWT) getLenOfOriginalString() int {
	return bwt.skipList[len(bwt.skipList)-1].openEndedInterval.end
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
	fb := strings.Builder{}
	for i := 0; i < len(prefixArray); i++ {
		fb.WriteByte(prefixArray[i][0])
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
