package bwt

import (
	"fmt"
	"math"

	"golang.org/x/exp/slices"
)

// waveletTree datastructure that allows us to
// conduct RSA queries on strings.
type waveletTree struct {
	root  *node
	alpha []charInfo
}

// Access will return the ith character of the original
// string used to build the waveletTree
func (wt waveletTree) Access(i int) byte {
	curr := wt.root
	for !curr.isLeaf() {
		bit := curr.data.Access(i)
		i = curr.data.Rank(bit, i)
		if bit {
			curr = curr.right
		} else {
			curr = curr.left
		}
	}
	return curr.char
}

// Rank allows us to get the rank of a specified character in
// the original string
func (wt waveletTree) Rank(char byte, i int) int {
	curr := wt.root
	ci := wt.lookupCharInfo(char)
	level := 0
	var rank int
	for !curr.isLeaf() {
		pathBit := ci.path.getBit(ci.path.len() - 1 - level)
		rank = curr.data.Rank(pathBit, i)
		if pathBit {
			curr = curr.right
		} else {
			curr = curr.left
		}
		level++
		i = rank
	}
	return rank
}

// Select allows us to get the corresponding posisiton of a character
// in the original string given its rank.
func (wt waveletTree) Select(char byte, rank int) int {
	curr := wt.root
	ci := wt.lookupCharInfo(char)
	level := 0

	for !curr.isLeaf() {
		pathBit := ci.path.getBit(ci.path.len() - 1 - level)
		if pathBit {
			curr = curr.right
		} else {
			curr = curr.left
		}
		level++
	}

	for curr.parent != nil {
		curr = curr.parent
		level--
		pathBit := ci.path.getBit(ci.path.len() - 1 - level)
		rank, ok := curr.data.Select(pathBit, rank)
		if !ok {
			msg := fmt.Sprintf("could not find a correspodning bit for node.Select(%t, %d) for characterInfo %+v", pathBit, rank, ci)
			panic(msg)
		}
	}

	return rank
}

func (wt waveletTree) lookupCharInfo(char byte) charInfo {
	for i := range wt.alpha {
		if wt.alpha[i].char == char {
			return wt.alpha[i]
		}
	}
	msg := fmt.Sprintf("could not find character %s in alphabet %+v. this should not be possible and indicates that the WaveletTree is malformed", string(char), wt.alpha)
	panic(msg)
}

type node struct {
	data   rsaBitVector
	char   byte
	parent *node
	left   *node
	right  *node
}

func (n node) isLeaf() bool {
	return n.char != 0
}

type charInfo struct {
	char    byte
	maxRank int
	path    bitvector
}

func NewWaveletTreeFromString(str string) waveletTree {
	bytes := []byte(str)

	alpha := getCharInfoDescByRank(bytes)
	root := buildWaveletTree(0, alpha, bytes)

	return waveletTree{
		root:  root,
		alpha: alpha,
	}
}

func buildWaveletTree(currentLevel int, alpha []charInfo, bytes []byte) *node {
	if len(alpha) == 0 {
		return nil
	}

	if len(alpha) == 1 {
		return &node{char: alpha[0].char}
	}

	leftAlpha, rightAlpha := partitionAlpha(currentLevel, alpha)

	var leftBytes []byte
	var rightBytes []byte

	bv := newBitVector(len(bytes))
	for i := range bytes {
		if isInAlpha(rightAlpha, bytes[i]) {
			bv.setBit(i, true)
			rightBytes = append(rightBytes, bytes[i])
		} else {
			leftBytes = append(leftBytes, bytes[i])
		}
	}

	root := &node{
		data: newRSABitVectorFromBitVector(bv),
	}

	leftTree := buildWaveletTree(currentLevel+1, leftAlpha, leftBytes)
	rightTree := buildWaveletTree(currentLevel+1, rightAlpha, rightBytes)

	root.left = leftTree
	root.right = rightTree

	if leftTree != nil {
		leftTree.parent = root
	}
	if rightTree != nil {
		rightTree.parent = root
	}

	return root
}

func isInAlpha(alpha []charInfo, b byte) bool {
	for _, a := range alpha {
		if a.char == b {
			return true
		}
	}
	return false
}

// partitionAlpha partitions the alaphabet in half based on whether its corresponding path bit
// is a 0 or 1. 0 with comprise the left tree while 1 will comprise the right. The alphabet
// should be sorted in such a way that we remove the most amount of characters nearest to the
// root of the tree to reduce the memory footprint as much as possible.
func partitionAlpha(currentLevel int, alpha []charInfo) (left []charInfo, right []charInfo) {
	for _, a := range alpha {
		if a.path.getBit(a.path.len() - 1 - currentLevel) {
			right = append(right, a)
		} else {
			left = append(left, a)
		}
	}

	return left, right
}

// getCharInfoDescByRank takes in the bytes of the original
// string and return a sorted list of character metadata descending
// by rank. The character metadata is important for building the rest
// of the tree along with quering it later on. The sorting is important
// because this allows us to build the tree in the most memory efficient
// way since the characters with the greatest counts will be removed first
// before build the subsequent nodes in the lower levels.
// NOTE: alphabets are expected to be small for real usecases
func getCharInfoDescByRank(b []byte) []charInfo {
	ranks := make(map[byte]int)
	for i := 0; i < len(b); i++ {
		if _, ok := ranks[b[i]]; ok {
			ranks[b[i]] += 1
		} else {
			ranks[b[i]] = 0
		}
	}

	var sortedInfo []charInfo
	for k := range ranks {
		sortedInfo = append(sortedInfo, charInfo{char: k, maxRank: ranks[k]})
	}

	slices.SortFunc(sortedInfo, func(a, b charInfo) bool {
		if a.maxRank == b.maxRank {
			return a.char < b.char
		}
		return a.maxRank > b.maxRank
	})

	numOfBits := getTreeHeight(sortedInfo)
	for i := range sortedInfo {
		bv := newBitVector(numOfBits)
		encodeCharPathIntoBitVector(bv, uint64(i))
		sortedInfo[i].path = bv
	}

	return sortedInfo
}

// encodeCharPathIntoBitVector important metadata to understand
// which character we are woring with in a given path in the tree.
// For example, given the alphabet A B C D E F G, a possible encoding is:
// A: 000
// B: 001
// C: 010
// D: 011
// E: 100
// F: 101
// G: 110
// H: 111
//
// If we wanted to get to the leaf that represent the character D, we'd
// take the path:
//
//	   root
//	  /
//	left
//	  \
//	 right
//	    \
//	   right
func encodeCharPathIntoBitVector(bv bitvector, n uint64) {
	shift := 0
	for n>>shift > 0 {
		if n>>shift%2 == 1 {
			bv.setBit(bv.len()-1-shift, true)
		} else {
			bv.setBit(bv.len()-1-shift, false)
		}
		shift++
	}
}

func getTreeHeight(alpha []charInfo) int {
	return int(math.Log2(float64(len(alpha)))) + 1
}
