package bwt

import (
	"math"

	"golang.org/x/exp/slices"
)

type waveletTree struct {
	root  *node
	alpha []charInfo
}

// TODO: figure out empty nodes case
// TODO: figure out out of bounds case
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

// TODO: deal with bad lookup char
// TODO: deal with somehow bad path
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

func (wt waveletTree) lookupCharInfo(char byte) charInfo {
	for i := range wt.alpha {
		if wt.alpha[i].char == char {
			return wt.alpha[i]
		}
	}
	panic("better messaging or handling")
}

// TODO: talk about how we could probably greaty improve performance with one big bit vector that
// represents the whole tree by concatenation the level order traversal of each node's bits
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

// TODO: talk about arranging OG alpha such that we minimize memory
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

func getLeft(nodePos int) int {
	return nodePos*2 + 1
}

func getRight(nodePos int) int {
	return nodePos*2 + 2
}

func getParent(nodePos int) int {
	return (nodePos + 1) / 2
}

// alphabets are expected to be small for real usecases
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

	slices.SortStableFunc(sortedInfo, func(a, b charInfo) bool {
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
