package bwt

import "math/bits"

// TODO: doc what rsa is, why these DSAs, and why we take in a bit vector
// TODO: clarks select
type rsaBitVector struct {
	bv                  bitvector
	totalOnesRank       int
	jrc                 []chunk
	jrSubChunksPerChunk int
	jrBitsPerChunk      int
	jrBitsPerSubChunk   int
	oneSelectMap        map[int]int
	zeroSelectMap       map[int]int
}

// TODO: talk about why bv should never be modidifed after building the RSA bit vector
func newRSABitVectorFromBitVector(bv bitvector) rsaBitVector {
	jacobsonRankChunks, jrSubChunksPerChunk, jrBitsPerSubChunk, totalOnesRank := buildJacobsonRank(bv)
	ones, zeros := buildSelectMaps(bv)

	return rsaBitVector{
		bv:                  bv,
		totalOnesRank:       totalOnesRank,
		jrc:                 jacobsonRankChunks,
		jrSubChunksPerChunk: jrSubChunksPerChunk,
		jrBitsPerChunk:      jrSubChunksPerChunk * jrBitsPerSubChunk,
		jrBitsPerSubChunk:   jrBitsPerSubChunk,
		oneSelectMap:        ones,
		zeroSelectMap:       zeros,
	}
}

// Rank returns the rank of the given value up to, but not including
// the ith bit. We count Rank starting a 0.
// For Example:
// Given the bitvector 001000100001
// Rank(true, 8) = 1
// Rank(false, 8) = 5
func (rsa rsaBitVector) Rank(val bool, i int) int {
	if i > rsa.bv.len()-1 {
		if val {
			return rsa.totalOnesRank
		}
		return rsa.bv.len() - rsa.totalOnesRank
	}

	chunkPos := (i / rsa.jrBitsPerChunk)
	chunk := rsa.jrc[chunkPos]

	subChunkPos := (i % rsa.jrBitsPerChunk) / rsa.jrBitsPerSubChunk
	subChunk := chunk.subChunks[subChunkPos]

	bitOffset := i % rsa.jrBitsPerSubChunk

	bitSet := rsa.bv.getBitSet(chunkPos*rsa.jrSubChunksPerChunk + subChunkPos)

	shiftRightAmount := uint64(rsa.jrBitsPerSubChunk - bitOffset)
	if val {
		remaining := bitSet >> shiftRightAmount
		return chunk.onesCumulativeRank + subChunk.onesCumulativeRank + bits.OnesCount64(remaining)
	}
	remaining := ^bitSet >> shiftRightAmount

	// cumulative ranks for 0 should just be the sum of the compliment of cumulative ranks for 1
	return (chunkPos*rsa.jrBitsPerChunk - chunk.onesCumulativeRank) + (subChunkPos*rsa.jrBitsPerSubChunk - subChunk.onesCumulativeRank) + bits.OnesCount64(remaining)
}

// Select returns the the position of the given value of a specified Rank
// For Example:
// Given the bitvector 001000100001
// Select(true, 1) = 6
// Rank(false, 5) = 7
func (rsa rsaBitVector) Select(val bool, rank int) (i int, ok bool) {
	if val {
		i, ok := rsa.oneSelectMap[rank]
		return i, ok
	} else {
		i, ok := rsa.zeroSelectMap[rank]
		return i, ok
	}
}

// Access returns the value of a bit at a given offset
func (rsa rsaBitVector) Access(i int) bool {
	return rsa.bv.getBit(i)
}

type chunk struct {
	subChunks          []subChunk
	onesCumulativeRank int
}

type subChunk struct {
	onesCumulativeRank int
}

// TODO: talk about easy to read instead vs perf
func buildJacobsonRank(inBv bitvector) (jacobsonRankChunks []chunk, numOfSubChunksPerChunk, numOfBitsPerSubChunk, totalRank int) {
	// TODO: talk about why this is probably good enough, improves as n grows, gets worse as n gets smaller, and how this fits into a machine instruction, and how this is "simple"
	numOfSubChunksPerChunk = 4

	totalRank = 0
	chunkCumulativeRank := 0
	subChunkCumulativeRank := 0

	var currSubChunks []subChunk
	for i := range inBv.bits {
		if len(currSubChunks) == numOfSubChunksPerChunk {
			jacobsonRankChunks = append(jacobsonRankChunks, chunk{
				subChunks:          currSubChunks,
				onesCumulativeRank: chunkCumulativeRank,
			})

			chunkCumulativeRank += subChunkCumulativeRank

			currSubChunks = nil
			subChunkCumulativeRank = 0
		}
		currSubChunks = append(currSubChunks, subChunk{
			onesCumulativeRank: subChunkCumulativeRank,
		})

		onesCount := bits.OnesCount64(inBv.getBitSet(i))
		subChunkCumulativeRank += onesCount
		totalRank += onesCount
	}

	if currSubChunks != nil {
		jacobsonRankChunks = append(jacobsonRankChunks, chunk{
			subChunks:          currSubChunks,
			onesCumulativeRank: chunkCumulativeRank,
		})
	}

	return jacobsonRankChunks, numOfSubChunksPerChunk, wordSize, totalRank
}

// TODO: talk about how this could be improved memory wise. Talk about how clarks select exists, but keeping it "simple for now" but maybe worth
// making succinct later
func buildSelectMaps(inBv bitvector) (oneSelectMap, zeroSelectMap map[int]int) {
	oneSelectMap = make(map[int]int)
	zeroSelectMap = make(map[int]int)
	oneCount := 0
	zeroCount := 0
	for i := 0; i < inBv.len(); i++ {
		bit := inBv.getBit(i)
		if bit {
			oneSelectMap[oneCount] = i
			oneCount++
		} else {
			zeroSelectMap[zeroCount] = i
			zeroCount++
		}
	}

	return oneSelectMap, zeroSelectMap
}
