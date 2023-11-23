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
	capacity := getCapacityNeededForNumberOfBits(initialNumberOfBits)
	bits := make([]uint, capacity)
	return bitvector{
		bits:             bits,
		capacityInChunks: capacity,
		numberOfBits:     initialNumberOfBits,
	}
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
	if getCapacityNeededForNumberOfBits(nextNumberOfBits) <= b.capacityInChunks {
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

	nextCapacity := getCapacityNeededForNumberOfBits(numOfBitsForNextCapacity)

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

func getCapacityNeededForNumberOfBits(n int) int {
	return int(math.Ceil(float64(n) / wordSize))
}

// TODO: doc what rsa is, why these DSAs, and why we take in a bit vector
type RSABitVector struct {
	jacobsonRank []chunk
	clarkSelect  []bitvector
}

type chunk struct {
	bits               bitvector
	onesCumulativeRank int
}

func newRSABitVector(b bitvector) RSABitVector {
	return RSABitVector{}
}

// TODO: doc how this is ugly and building is probably bad. talk about chunk size
func buildJacobsonRank(inBv bitvector) (int, []chunk) {
	// TODO: doc magic numbers and doc that this will always be a natural number
	uLen := uint(inBv.len())
	leading1Offset := bits.UintSize - bits.LeadingZeros(uLen)
	perfectSquare := int(uint(1) << uint(leading1Offset))
	chunkSize := int(math.Pow(math.Log2(float64(perfectSquare)), 2))
	numChunks := inBv.len() / chunkSize

	// TODO: doc why we have the plus 1
	jacobsonRank := make([]chunk, numChunks+1)
	onesCount := 0
	for i := 0; i < numChunks; i++ {
		chunkBv := newBitVector(chunkSize)
		for j := 0; j < chunkSize; j++ {
			val := inBv.getBit(i*wordSize + j)
			if val {
				onesCount++
			}
			chunkBv.setBit(j, val)
		}
		jacobsonRank[i] = chunk{
			bits:               chunkBv,
			onesCumulativeRank: onesCount,
		}
	}

	// TODO: doc the last chunk
	lastChunkSize := inBv.len() % int(perfectSquare)
	lastChunkBv := newBitVector(lastChunkSize)
	for i := 0; i < lastChunkSize; i++ {
		val := inBv.getBit(numChunks*wordSize + i)
		if val {
			onesCount++
		}
		lastChunkBv.setBit(i, val)
	}
	jacobsonRank[len(jacobsonRank)-1] = chunk{
		bits:               lastChunkBv,
		onesCumulativeRank: onesCount,
	}

	return chunkSize, jacobsonRank
}

// TODO: 15:16 sub chunk impl
// TODO: 17:25 sub chunk calc of rank lower than machine word impl.
