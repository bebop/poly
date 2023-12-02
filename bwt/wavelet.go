package bwt

import (
	"math"

	"golang.org/x/exp/slices"
)

type WaveletTree struct {
	nodes []node
	alpha []charInfo
}

type node struct {
	data rsaBitVector
	char byte
}

func (n node) isLeaf() bool {
	return n.char != 0
}

type charInfo struct {
	char    byte
	maxRank int
	path    bitvector
}

func NewWaveletTreeFromString(str string) WaveletTree {
	bytes := []byte(str)

	alpha := getCharInfoDescByRank(bytes)
	nodes := buildWaveletTree(0, alpha, bytes)

	return WaveletTree{
		nodes: nodes,
	}
}

func buildWaveletTree(currentLevel int, alpha []charInfo, bytes []byte) []node {
	if len(alpha) == 0 {
		return nil
	}

	if len(alpha) == 1 {
		return []node{
			{char: alpha[0].char},
		}
	}

	if len(alpha) == 2 {
		return []node{
			{char: alpha[0].char},
			{char: alpha[1].char},
		}
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

	n := node{
		data: newRSABitVectorFromBitVector(bv),
	}

	leftTree := buildWaveletTree(currentLevel+1, leftAlpha, leftBytes)
	rightTree := buildWaveletTree(currentLevel+1, rightAlpha, rightBytes)

	tree := append([]node{n}, leftTree...)
	tree = append(tree, rightTree...)

	return tree
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
	return nodePos*2 + 1
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

	slices.SortFunc(sortedInfo, func(a, b charInfo) bool {
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
