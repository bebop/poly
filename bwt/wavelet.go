package bwt

import (
	"fmt"
	"math"

	"golang.org/x/exp/slices"
)

/*

For the waveletTree's usage, please read the its
method documentation. To understand what it is and how
it works for either curiosity or maintenance, then read below.

# WaveletTree Components

## The Character's Path Encoding

One important component is a character's path encoding.
Which character we are working with in a given path in the tree.
For example, given the alphabet A B C D E F G, a possible encoding is:

A: 000
B: 001
C: 010
D: 011
E: 100
F: 101
G: 110
H: 111

If we wanted to get to the leaf that represent the character D, we'd
take the path:

   root
  /
left
  \
 right
    \
   right

## The Data Represented at each node

Let us consider the sequence "bananas"
It has the alphabet b, a, n, s
Let's say it has the encoding:
a: 00
n: 01
b: 10
s: 11
and that 0 if left and 1 is right
We can represent this tree with bitvectors:

     0010101
    /       \
  1000      001
 /    \    /   \
a      n  b     s

If we translate each bit vector to its corresponding string, then it becomes:

     bananas
    /       \
  baaa      nns
 /    \    /   \
a      b  n     s

## RSA Intuition

RSA stands for (R)ank, (S)elect, (A)ccess.

From here you may be able to build some intuition as to how we can take RSA queries given
a characters path encoding and which character we'd like to Rank, Select, and Access.

### Rank Example

To get WaveletTree.Rank(a, 4) of bananas where a's encoding is 00
1. root.Rank(0, 4) of 0010101 is 2
2. Visit Left Child
3. child.Rank(0, 2) of 1000 is 1
4. Visit Left Child
5. return 1

### Select Example

To get WaveletTree.Select(n, 1) of bananas where n's encoding is 01
1. Go down to n's leaf using the path encoding is 01
2. Go back to n's leaf's parent
3. parent.Select(0, 1) of 001 is 1
4. Go to the next parent
5. parent.Select(1, 1) of 0010101 is 4
6. return 4 since we are at the root.

### Access Example

If you've reached this point, then you must really be trying to understand how the
waveletTree works. I recommend thinking through how access could work with the example
above. HINT: rank might help.

NOTE: The waveletTree does not literally have to be a tree. There are other forms that it may
exist in like the concatenation of order level representation of all its node's bitvectors...
as one example. Please reference the implementation if you'd like to understand how this
specific waveletTree works.

*/

// waveletTree is a data structure that allows us to index a sequence
// in a memory efficient way that allows us to conduct RSA, (R)ank (S)elect (A)ccess
// queries on strings. This is very useful in situations where you'd like to understand
// certain aspects of a sequence like:
// * the number of times a character appears
// * counting how the frequency of a character up to certain offset
// * locating characters of certain rank within the sequence
// * accessing the character at a given position
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

// Select allows us to get the corresponding position of a character
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
		nextRank, ok := curr.data.Select(pathBit, rank)
		if !ok {
			msg := fmt.Sprintf("could not find a corresponding bit for node.Select(%t, %d) for characterInfo %+v", pathBit, rank, ci)
			panic(msg)
		}
		rank = nextRank
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

// partitionAlpha partitions the alphabet in half based on whether its corresponding path bit
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
// of the tree along with querying it later on. The sorting is important
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
