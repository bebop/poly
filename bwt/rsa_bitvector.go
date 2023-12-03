package bwt

import "math/bits"

// TODO: doc what rsa is, why these DSAs, and why we take in a bit vector
// TODO: clarks select
type rsaBitVector struct {
	bv                bitvector
	jrc               []chunk
	jrBitsPerChunk    int
	jrBitsPerSubChunk int
	oneSelectMap      map[int]int
	zeroSelectMap     map[int]int
}

// TODO: talk about why bv should never be modidifed after building the RSA bit vector
func newRSABitVectorFromBitVector(bv bitvector) rsaBitVector {
	jacobsonRankChunks, jrBitsPerChunk, jrBitsPerSubChunk := buildJacobsonRank(bv)
	ones, zeros := buildSelectMaps(bv)

	return rsaBitVector{
		bv:                bv,
		jrc:               jacobsonRankChunks,
		jrBitsPerChunk:    jrBitsPerChunk,
		jrBitsPerSubChunk: jrBitsPerSubChunk,
		oneSelectMap:      ones,
		zeroSelectMap:     zeros,
	}
}

func (rsa rsaBitVector) Rank(val bool, i int) int {
	rsa.bv.checkBounds(i)

	chunkPos := (i / rsa.jrBitsPerChunk)
	chunk := rsa.jrc[chunkPos]

	subChunkPos := (i % rsa.jrBitsPerChunk) / rsa.jrBitsPerSubChunk
	subChunk := chunk.subChunks[subChunkPos]

	bitOffset := i % rsa.jrBitsPerSubChunk

	bitSet := rsa.bv.getBitSet(chunkPos*len(rsa.jrc) + subChunkPos)

	shiftRightAmount := uint64(rsa.jrBitsPerSubChunk - bitOffset)
	if val {
		remaining := bitSet >> shiftRightAmount
		return chunk.onesCumulativeRank + subChunk.onesCumulativeRank + bits.OnesCount64(remaining)
	}
	remaining := ^bitSet >> shiftRightAmount

	// cumulative ranks for 0 should just be the sum of the compliment of cumulative ranks for 1
	return (chunkPos*rsa.jrBitsPerChunk - chunk.onesCumulativeRank) + (subChunkPos*rsa.jrBitsPerSubChunk - subChunk.onesCumulativeRank) + bits.OnesCount64(remaining)
}

func (rsa rsaBitVector) Select(val bool, rank int) (i int, ok bool) {
	if val {
		i, ok := rsa.oneSelectMap[rank]
		return i, ok
	} else {
		i, ok := rsa.zeroSelectMap[rank]
		return i, ok
	}
}

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
func buildJacobsonRank(inBv bitvector) (jacobsonRankChunks []chunk, numOfSubChunksPerChunk, numOfBitsPerSubChunk int) {
	// TODO: talk about why this is probably good enough, improves as n grows, gets worse as n gets smaller, and how this fits into a machine instruction, and how this is "simple"
	numOfSubChunksPerChunk = 4

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
	}

	if currSubChunks != nil {
		jacobsonRankChunks = append(jacobsonRankChunks, chunk{
			subChunks:          currSubChunks,
			onesCumulativeRank: chunkCumulativeRank,
		})
	}

	return jacobsonRankChunks, numOfSubChunksPerChunk * wordSize, wordSize
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
			oneCount++
			oneSelectMap[oneCount] = i
		} else {
			zeroCount++
			zeroSelectMap[zeroCount] = i
		}
	}

	return oneSelectMap, zeroSelectMap
}
