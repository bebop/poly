package bwt

import (
	"math"
	"math/bits"
)

// TODO: talk about why this is
const wordSize = bits.UintSize

// TODO: document static size and why we differentiate between capacity and number of bits
type bitvector struct {
	bits             []uint
	capacityInChunks int
	numberOfBits     int
}

func newBitVector(initialNumberOfBits int) bitvector {
	capacity := getNumOfBitSetsNeededForNumOfBits(initialNumberOfBits)
	bits := make([]uint, capacity)
	return bitvector{
		bits:             bits,
		capacityInChunks: capacity,
		numberOfBits:     initialNumberOfBits,
	}
}

func (b bitvector) getNumOfBitSets() int {
	return getNumOfBitSetsNeededForNumOfBits(b.len())
}

func (b bitvector) getBitSet(i int) uint {
	return b.bits[i]
}

func (b bitvector) getBit(i int) bool {
	if i >= b.len() || i < 0 {
		panic("better out of bounds message")
	}

	chunkStart := i / wordSize
	offset := i % wordSize

	return (b.bits[chunkStart] & (uint(1) << offset)) != 0
}

func (b bitvector) setBit(i int, val bool) {
	if i >= b.len() || i < 0 {
		panic("better out of bounds message")
	}

	chunkStart := i / wordSize
	offset := i % wordSize

	if val {
		b.bits[chunkStart] |= uint(1) << offset
	} else {
		b.bits[chunkStart] &= ^(uint(1) << offset)
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

	nextBits := make([]uint, nextCapacity)
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
type RSABitVector struct {
	numOfBits         int
	jrc               []chunk
	jrBitsPerChunk    int
	jrBitsPerSubChunk int
	clarkSelect       []bitvector
}

func newRSABitVector(b bitvector) RSABitVector {
	jacobsonRankChunks, jrBitsPerChunk, jrBitsPerSubChunk := buildJacobsonRank(b)
	return RSABitVector{
		numOfBits:         b.len(),
		jrc:               jacobsonRankChunks,
		jrBitsPerChunk:    jrBitsPerChunk,
		jrBitsPerSubChunk: jrBitsPerSubChunk,
		clarkSelect:       []bitvector{},
	}
}

// TODO: doc and mention some bit math
func (rsa RSABitVector) rank(val bool, i int) int {
	chunkPos := (i / rsa.jrBitsPerChunk)
	chunk := rsa.jrc[chunkPos]

	subChunkPos := (i % rsa.jrBitsPerChunk) / rsa.jrBitsPerSubChunk
	subChunk := chunk.subChunks[subChunkPos]

	bitOffset := i % rsa.jrBitsPerSubChunk

	shiftRightAmount := uint(rsa.jrBitsPerSubChunk - bitOffset)
	if val {
		remaining := subChunk.bitSet >> shiftRightAmount
		return chunk.onesCumulativeRank + subChunk.onesCumulativeRank + bits.OnesCount(remaining)
	}
	remaining := ^subChunk.bitSet >> shiftRightAmount
	return (chunkPos*rsa.jrBitsPerChunk - chunk.onesCumulativeRank) + (subChunkPos*rsa.jrBitsPerSubChunk - subChunk.onesCumulativeRank) + bits.OnesCount(remaining)
}

type chunk struct {
	subChunks          []subChunk
	onesCumulativeRank int
}

type subChunk struct {
	bitSet             uint
	onesCumulativeRank int
}

// TODO: talk about easy to read instead vs perf
func buildJacobsonRank(inBv bitvector) (jacobsonRankChunks []chunk, numOfSubChunksPerChunk, numOfBitsPerSubChunk int) {
	// TODO: talk about why this is probably good enough, improves as n grows, gets worse as n gets smaller, and how this fits into a machine instruction, and how this is "simple"
	numOfSubChunksPerChunk = 4

	chunkCumulativeRank := 0
	subChunkCumulativeRank := 0

	var currSubChunks []subChunk
	for _, bitSet := range inBv.bits {
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
			bitSet:             bitSet,
			onesCumulativeRank: subChunkCumulativeRank,
		})

		onesCount := bits.OnesCount(bitSet)
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
