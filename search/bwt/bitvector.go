package bwt

import (
	"fmt"
	"math"
)

const wordSize = 64

// bitvector a sequence of 1's and 0's. You can also think
// of this as an array of bits. This allows us to encode
// data in a memory efficient manner.
type bitvector struct {
	bits         []uint64
	numberOfBits int
}

// newBitVector will return an initialized bitvector with
// the specified number of zeroed bits.
func newBitVector(initialNumberOfBits int) bitvector {
	capacity := getNumOfBitSetsNeededForNumOfBits(initialNumberOfBits)
	bits := make([]uint64, capacity)
	return bitvector{
		bits:         bits,
		numberOfBits: initialNumberOfBits,
	}
}

// getBitSet gets the while word as some offset from the
// bitvector. Useful if you'd prefer to work with the
// word rather than with individual bits.
func (b bitvector) getBitSet(bitSetPos int) uint64 {
	return b.bits[bitSetPos]
}

// getBit returns the value of the bit at a given offset
// True represents 1
// False represents 0
func (b bitvector) getBit(i int) bool {
	b.checkBounds(i)

	chunkStart := i / wordSize
	offset := i % wordSize

	return (b.bits[chunkStart] & (uint64(1) << (63 - offset))) != 0
}

// setBit sets the value of the bit at a given offset
// True represents 1
// False represents 0
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
		msg := fmt.Sprintf("access of %d is out of bounds for bitvector with length %d", i, b.len())
		panic(msg)
	}
}

func (b bitvector) len() int {
	return b.numberOfBits
}

func getNumOfBitSetsNeededForNumOfBits(n int) int {
	return int(math.Ceil(float64(n) / wordSize))
}
