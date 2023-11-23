package bwt

import "math"

const chunkSize = 8

// TODO: document static size and why we differentiate between capacity and number of bits
type bitvector struct {
	bits             []uint8
	capacityInChunks int
	numberOfBits     int
}

func newBitVector(initialNumberOfBits int) bitvector {
	capacity := getCapacityNeededForNumberOfBits(initialNumberOfBits)
	bits := make([]uint8, capacity)
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

	chunkStart := i / chunkSize
	offset := i % chunkSize

	return (b.bits[chunkStart] & (uint8(1) << offset)) != 0
}

func (b bitvector) setBit(i int, val bool) {
	if i >= b.len() || i < 0 {
		panic("better out of bounds message")
	}

	chunkStart := i / chunkSize
	offset := i % chunkSize

	if val {
		b.bits[chunkStart] |= uint8(1) << offset
	} else {
		b.bits[chunkStart] &= ^(uint8(1) << offset)
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
		numOfBitsForNextCapacity = int(math.Ceil(float64(b.numberOfBits) * 1.2))
		break
	case nextNumberOfBits >= factor1point5Threshold:
		numOfBitsForNextCapacity = int(math.Ceil(float64(b.numberOfBits) * 1.5))
		break
	default:
		numOfBitsForNextCapacity = b.numberOfBits * 2
	}

	nextCapacity := getCapacityNeededForNumberOfBits(numOfBitsForNextCapacity)

	nextBits := make([]uint8, nextCapacity)
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
	return int(math.Ceil(float64(n) / 8.0))
}
