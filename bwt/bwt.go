package bwt

import (
	"fmt"
	"strings"

	"golang.org/x/exp/slices"
)

const nullChar = "0"

// BWT Burrow Wheeler Transform
// Compresses and Indexes a given sequence so that it can be
// be used for search, alignment, and text extraction. This is
// useful for sequences so large that it would be beneficial
// to reduce its memory footprint while also maintaining a way
// to analyze and work with the sequence.
type BWT struct {
	// firstColumnSkipList is the first column of the BWT. It is
	// represented as a list of skipEntries because the first column of
	// the BWT is always lexographically ordered. This saves time and memory.
	firstColumnSkipList []skipEntry
	// lastCoulmn last column of the BWT- the actual textual representation
	// of the BWT.
	lastCoulmn waveletTree
	// suffixArray an array that allows us to map a posistion in the first
	// column to a position in the original sequence. This is needed to be
	// able to extract text from the BWT.
	suffixArray []int
}

// Count represents the number of times the provided pattern
// shows up in the original sequence.
func (bwt BWT) Count(pattern string) int {
	searchRange := bwt.lfSearch(pattern)
	return searchRange.end - searchRange.start
}

// Locate returns a list of offsets at which the begging
// of the provided pattern occurrs in the original
// sequence.
func (bwt BWT) Locate(pattern string) []int {
	searchRange := bwt.lfSearch(pattern)
	if searchRange.start >= searchRange.end {
		return nil
	}

	numOfOffsets := searchRange.end - searchRange.start
	offsets := make([]int, numOfOffsets)
	for i := 0; i < numOfOffsets; i++ {
		offsets[i] = bwt.suffixArray[searchRange.start+i]
	}

	return offsets
}

// Extract this allows us to extract parts of the original
// sequence from the BWT.
// start is the begging of the range of text to extract inclusive.
// end is the end of the range of text to extract exclusive.
// If either start or end are out of bounds, Extract will panic.
func (bwt BWT) Extract(start, end int) string {
	if end > bwt.getLenOfOriginalStringWithNullChar()-1 {
		msg := fmt.Sprintf("end [%d] exceeds the max range of the BWT [%d]", end, bwt.getLenOfOriginalStringWithNullChar()-1)
		panic(msg)
	}
	if start < 0 {
		msg := fmt.Sprintf("start [%d] exceeds the min range of the BWT [0]", start)
		panic(msg)
	}

	strB := strings.Builder{}
	for i := start; i < end; i++ {
		fPos := bwt.getFCharPosFromOriginalSequenceCharPos(i)
		skip := bwt.lookupSkipByOffset(fPos)
		strB.WriteByte(skip.char)
	}
	return strB.String()
}

// Len return the length of the sequence used to build the BWT
func (bwt BWT) Len() int {
	return bwt.getLenOfOriginalStringWithNullChar() - 1
}

// getFCharPosFromOriginalSequenceCharPos looks up mapping from the original position
// of the sequence to its corresponding posisiton in the First Column of the BWT
func (bwt BWT) getFCharPosFromOriginalSequenceCharPos(originalPos int) int {
	for i := range bwt.suffixArray {
		if bwt.suffixArray[i] == originalPos {
			return i
		}
	}
	panic("Unable to find the corresponding original position for a character in the original sequence in the suffix array. This should not be possible and indicates a malformed BWT.")
}

// lfSearch LF Search- Last First Search.
// Finds the valid range within the BWT index where the provided pattern is possible.
// If the final range is <= 0, then the pattern does not exist in the original sequence.
func (bwt BWT) lfSearch(pattern string) interval {
	searchRange := interval{start: 0, end: bwt.getLenOfOriginalStringWithNullChar()}
	for i := 0; i < len(pattern); i++ {
		if searchRange.end-searchRange.start <= 0 {
			return interval{}
		}

		c := pattern[len(pattern)-1-i]
		skip, ok := bwt.lookupSkipByChar(c)
		if !ok {
			return interval{}
		}
		searchRange.start = skip.openEndedInterval.start + bwt.lastCoulmn.Rank(c, searchRange.start)
		searchRange.end = skip.openEndedInterval.start + bwt.lastCoulmn.Rank(c, searchRange.end)
	}
	return searchRange
}

// lookupSkipByChar looks up a skipEntry by its character in the First Coulmn
func (bwt BWT) lookupSkipByChar(c byte) (entry skipEntry, ok bool) {
	for i := range bwt.firstColumnSkipList {
		if bwt.firstColumnSkipList[i].char == c {
			return bwt.firstColumnSkipList[i], true
		}
	}
	return skipEntry{}, false
}

// lookupSkipByOffset looks up a skipEntry based off of an
// offset of the Fist Coulmn of the BWT.
func (bwt BWT) lookupSkipByOffset(offset int) skipEntry {
	if offset > bwt.getLenOfOriginalStringWithNullChar()-1 {
		msg := fmt.Sprintf("offset [%d] exceeds the max bound of the BWT [%d]", offset, bwt.getLenOfOriginalStringWithNullChar()-1)
		panic(msg)
	}
	if offset < 0 {
		msg := fmt.Sprintf("offset [%d] exceeds the min bound of the BWT [0]", offset)
		panic(msg)
	}

	for i := range bwt.firstColumnSkipList {
		if bwt.firstColumnSkipList[i].openEndedInterval.start <= offset && offset < bwt.firstColumnSkipList[i].openEndedInterval.end {
			return bwt.firstColumnSkipList[i]
		}
	}
	panic("figure out what to do here")
}

func (bwt BWT) getLenOfOriginalStringWithNullChar() int {
	return bwt.firstColumnSkipList[len(bwt.firstColumnSkipList)-1].openEndedInterval.end
}

type interval struct {
	start int
	end   int
}

type skipEntry struct {
	char byte
	// openEndedInterval start is inclusive and end is exclusive
	openEndedInterval interval
}

// New returns a BWT of the provided sequence
// The provided sequence must not contain the nullChar
// defined in this package. If it does, New will return
// an error.
func New(sequence string) (BWT, error) {
	if strings.Contains(sequence, nullChar) {
		return BWT{}, fmt.Errorf("Provided sequence contains the nullChar %s. BWT cannot be constructed", nullChar)
	}

	sequence += nullChar

	prefixArray := make([]string, len(sequence))
	for i := 0; i < len(sequence); i++ {
		prefixArray[i] = sequence[len(sequence)-i-1:]
	}

	// TODO: at the time of writing, the nullChar is 0, this is to ensure correctness in most cases.
	// Do we want to roll our own sorting so we can make sure whatever is defined as the nullChar
	// will absolutely be defined as the least?
	slices.Sort(prefixArray)

	suffixArray := make([]int, len(sequence))
	lastColBuilder := strings.Builder{}
	for i := 0; i < len(prefixArray); i++ {
		currChar := sequence[getBWTIndex(len(sequence), len(prefixArray[i]))]
		lastColBuilder.WriteByte(currChar)

		suffixArray[i] = len(sequence) - len(prefixArray[i])
	}
	fb := strings.Builder{}
	for i := 0; i < len(prefixArray); i++ {
		fb.WriteByte(prefixArray[i][0])
	}

	return BWT{
		firstColumnSkipList: buildSkipList(prefixArray),
		lastCoulmn:          NewWaveletTreeFromString(lastColBuilder.String()),
		suffixArray:         suffixArray,
	}, nil
}

// buildSkipList compressed the First Column of the BWT into a skip list
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

// getBWTIndex returns the position of the character from the sequence used to build the BWT
// that corresponds the last character that would exist in the entry of the prefixArray that
// would be the last character if we were actually doing full rotations
func getBWTIndex(lenOfSequenceBeingBuilt, lenOfSuffixArrayVisited int) int {
	bwtCharIndex := lenOfSequenceBeingBuilt - lenOfSuffixArrayVisited - 1
	if bwtCharIndex == -1 {
		bwtCharIndex = lenOfSequenceBeingBuilt - 1
	}
	return bwtCharIndex
}
