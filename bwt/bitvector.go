package bwt

import (
	"math"
)

// TODO: talk about why this is and why we approximate things to make them "simple enough"
const wordSize = 64

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

	return (b.bits[chunkStart] & (uint64(1) << (63 - offset))) != 0
}

func (b bitvector) setBit(i int, val bool) {
	b.checkBounds(i)

	chunkStart := i / wordSize
	offset := i % wordSize

	if val {
		b.bits[chunkStart] |= uint64(1) << (63 - offset)
	} else {
		b.bits[chunkStart] &= ^(uint64(1) << (63 - offset))
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
