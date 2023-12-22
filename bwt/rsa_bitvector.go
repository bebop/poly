package bwt

import "math/bits"

// rsaBitVector allows us to perform RSA: (R)ank, (S)elect, and (A)ccess
// queries in a memory performant and memory compact way.
// To learn about how Rank, Select, and Access work, take a look at the
// examples in each respective method.
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

// newRSABitVectorFromBitVector allows us to build the auxiliary components
// needed to perform RSA queries on top of the provided bitvector.
// WARNING: Do not modify the underlying bitvector. The rsaBitvector will
// get out of sync with the original bitvector.
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
// the ith bit.
// For Example:
// Given the bitvector 001000100001
// Rank(true, 1) = 0
// Rank(true, 2) = 0
// Rank(true, 3) = 1
// Rank(true, 8) = 2
// Rank(false, 8) = 6
func (rsa rsaBitVector) Rank(val bool, i int) int {
	if i == rsa.bv.len() {
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

// Select returns the position of the given value with the provided Rank
// For Example:
// Given the bitvector 001000100001
// Select(true, 1) = 2
// Rank(false, 5) = 5
// Rank(false, 1) = 1
// Rank(false, 0) = 0
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

/*
buildJacobsonRank Jacobson rank is a succinct data structure. This allows us to represent something
normally would require O(N) worth of memory with less that N memory. Jacobson Rank allows for
sub linear growth. Jacobson rank also allows us to lookup rank for some value of a bitvector in O(1)
time. Theoretically, Jacobson Rank Requires:
1. Creating log(N) "Chunks"
2. Creating 2log(N) "Sub Chunks"
3. Having "Sub Chunks" be 0.5log(N) in length
4. For each "Chunk", store the cumulative rank of set bits relative to the overall bitvector
5. For each "Sub Chunk", store the cumulative rank of set bits relative to the parent "Chunk"
6. We can One's count the N bit word if possible. We will only consider this possibility :)

For simplicity and all around decent results, we just have "Sub Chunks" of size 64 bits.

It is O(1) because given some offset i, all we have to do is calculate rank is:
rank = CumulativeRank(ChunkOfi(i))) + CumulativeRank(SubChunkOfi(i))) + OnesCount(SubChunkOfi(i))

To understand why it is sub linear in space, you can refer to Ben Langmead and other literature that
describes the space complexity.
https://www.youtube.com/watch?v=M1sUZxXVjG8&list=PL2mpR0RYFQsADmYpW2YWBrXJZ_6EL_3nu&index=7
*/
func buildJacobsonRank(inBv bitvector) (jacobsonRankChunks []chunk, numOfSubChunksPerChunk, numOfBitsPerSubChunk, totalRank int) {
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

// This is not good. We should find a better means of select- like Clark's Select
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

	// Account for the case where we need to find the
	// position for the max rank for both 0's and 1's
	oneSelectMap[oneCount] = inBv.len()
	zeroSelectMap[zeroCount] = inBv.len()

	return oneSelectMap, zeroSelectMap
}
