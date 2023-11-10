package bwt

import (
	"golang.org/x/exp/slices"
)

const nullChar = "$"

type window struct {
	char  byte
	start int
	end   int
}

func (w window) includes(innerWindow window) bool {
	return w.start <= innerWindow.start && innerWindow.end <= w.end
}

// BWT Burrow Wheeler Transform
// Data structure that compactly represents any sequence of characters and
// allows for sub sequence querying.
type BWT struct {
	// First col
	f []window
	// Last col
	l map[byte][]window
	// index of the original sequence in the suffix array from BTW construction
	indexOfOriginalSequenceFromSuffixArray int
}

func New(sequence string) BWT {
	f, l, idx := build(sequence)
	return BWT{
		f:                                      f,
		l:                                      l,
		indexOfOriginalSequenceFromSuffixArray: idx,
	}
}

func (b BWT) QueryExistence(substr string) bool {
	for i := len(substr) - 1; i >= 0; i-- {
		win, ok := b.getFWindow(substr[i])
		if !ok {
			return false
		}

		if i == 0 && ok {
			return true
		}

		valid := b.charInLExistsInFWindow(win, substr[i-1])

		if !valid {
			return false
		}
	}

	// shouldn't be getting here
	return false
}

func (b BWT) charInLExistsInFWindow(w window, char byte) bool {
	if windows, ok := b.l[char]; ok {
		for i := range windows {
			if w.includes(windows[i]) {
				return true
			}
		}
	}
	return false
}

// Alphabets should be small
func (b BWT) getFWindow(char byte) (w window, ok bool) {
	for i := range b.f {
		if b.f[i].char == char {
			return b.f[i], true
		}
	}
	return window{}, false
}

func build(s string) (f []window, l map[byte][]window, indexOfOriginalSequenceInSuffixArray int) {
	s += nullChar
	prefixArray := make([]string, len(s))
	for i := 0; i < len(s); i++ {
		prefixArray[i] = s[len(s)-i-1:]
	}

	slices.Sort(prefixArray)

	l = make(map[byte][]window)
	prevFChar := prefixArray[0][0]
	prevFWin := window{char: prevFChar, start: 0}
	prevLChar := s[getBWTIndex(len(s), len(prefixArray[0]))]
	prevLWin := window{char: prevLChar, start: 0}
	for i := 1; i < len(prefixArray); i++ {
		currFChar := prefixArray[i][0]
		if prevFChar != currFChar {
			prevFWin.end = i - 1
			f = append(f, prevFWin)
			prevFChar = currFChar
			prevFWin = window{char: currFChar, start: i}
		}

		currLChar := s[getBWTIndex(len(s), len(prefixArray[i]))]
		if prevLChar != currLChar {
			prevLWin.end = i - 1
			if _, ok := l[prevLChar]; ok {
				l[prevLChar] = append(l[prevLChar], prevLWin)
			} else {
				l[prevLChar] = []window{prevLWin}
			}
			prevLChar = currLChar
			prevLWin = window{char: currLChar, start: i}
		}
		if len(s) == len(prefixArray[i]) {
			indexOfOriginalSequenceInSuffixArray = i
		}
	}
	prevFWin.end = len(s) - 1
	f = append(f, prevFWin)
	prevLWin.end = len(s) - 1
	if _, ok := l[prevLChar]; ok {
		l[prevLChar] = append(l[prevLChar], prevLWin)
	} else {
		l[prevLChar] = []window{prevLWin}
	}
	if indexOfOriginalSequenceInSuffixArray == 0 {
		indexOfOriginalSequenceInSuffixArray = len(s) - 1
	}

	return f, l, indexOfOriginalSequenceInSuffixArray
}

func getBWTIndex(lenOfSequenceBeingBuilt, lenOfSuffixArrayVisited int) int {
	bwtCharIndex := lenOfSequenceBeingBuilt - lenOfSuffixArrayVisited - 1
	if bwtCharIndex == -1 {
		bwtCharIndex = lenOfSequenceBeingBuilt - 1
	}
	return bwtCharIndex
}
