package bwt

import (
	"math"
	"math/bits"
)

// TODO: talk about why this is and why we approximate things to make them "simple enough"
const wordSize = 64

// TODO: document static size and why we differentiate between capacity and number of bits
type bitvector struct {
	bits             []uint64
	capacityInChunks int
	numberOfBits     int
}

func newBitVector(initialNumberOfBits int) bitvector {
	capacity := getNumOfBitSetsNeededForNumOfBits(initialNumberOfBits)
	bits := make([]uint64, capacity)
	return bitvector{
		bits:             bits,
		capacityInChunks: capacity,
		numberOfBits:     initialNumberOfBits,
	}
}

func (b bitvector) getNumOfBitSets() int {
	return getNumOfBitSetsNeededForNumOfBits(b.len())
}

func (b bitvector) getBitSet(i int) uint64 {
	return b.bits[i]
}

func (b bitvector) getBit(i int) bool {
	b.checkBounds(i)

	chunkStart := i / wordSize
	offset := i % wordSize

	return (b.bits[chunkStart] & (uint64(1) << offset)) != 0
}

func (b bitvector) setBit(i int, val bool) {
	b.checkBounds(i)

	chunkStart := i / wordSize
	offset := i % wordSize

	if val {
		b.bits[chunkStart] |= uint64(1) << offset
	} else {
		b.bits[chunkStart] &= ^(uint64(1) << offset)
	}
}
func (b bitvector) checkBounds(i int) {
	if i >= b.len() || i < 0 {
		panic("better out of bounds message")
	}
}

const factor1point2Threshold = 1e9
const factor1point5Threshold = 1e6

func (b *bitvector) push(val bool) {
	previousNumberOfBits := b.numberOfBits
	nextNumberOfBits := previousNumberOfBits + 1
	if getNumOfBitSetsNeededForNumOfBits(nextNumberOfBits) <= b.capacityInChunks {
		b.numberOfBits = nextNumberOfBits
		b.setBit(previousNumberOfBits, val)
		return
	}

	var numOfBitsForNextCapacity int
	switch true {
	case nextNumberOfBits >= factor1point2Threshold:
		numOfBitsForNextCapacity = int(math.Ceil(float64(previousNumberOfBits) * 1.2))
		break
	case nextNumberOfBits >= factor1point5Threshold:
		numOfBitsForNextCapacity = int(math.Ceil(float64(previousNumberOfBits) * 1.5))
		break
	default:
		numOfBitsForNextCapacity = previousNumberOfBits * 2
	}

	nextCapacity := getNumOfBitSetsNeededForNumOfBits(numOfBitsForNextCapacity)

	nextBits := make([]uint64, nextCapacity)
	copy(b.bits, nextBits)
	b.bits = nextBits

	b.numberOfBits = nextNumberOfBits
	b.capacityInChunks = nextCapacity

	b.setBit(previousNumberOfBits, val)
}

func (b bitvector) len() int {
	return b.numberOfBits
}

func (b bitvector) capacity() int {
	return b.capacityInChunks
}

func getNumOfBitSetsNeededForNumOfBits(n int) int {
	return int(math.Ceil(float64(n) / wordSize))
}

// TODO: doc what rsa is, why these DSAs, and why we take in a bit vector
// TODO: clarks select
type rsaBitVector struct {
	bv                bitvector
	jrc               []chunk
	jrBitsPerChunk    int
	jrBitsPerSubChunk int
}

// TODO: talk about why bv should never be modidifed after building the RSA bit vector
func newRSABitVector(bv bitvector) rsaBitVector {
	jacobsonRankChunks, jrBitsPerChunk, jrBitsPerSubChunk := buildJacobsonRank(bv)
	return rsaBitVector{
		bv:                bv,
		jrc:               jacobsonRankChunks,
		jrBitsPerChunk:    jrBitsPerChunk,
		jrBitsPerSubChunk: jrBitsPerSubChunk,
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

func (rsa rsaBitVector) Select(i int) bool {
	return true
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
