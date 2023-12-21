package bwt

import (
	"testing"
)

type GetBitTestCase struct {
	position int
	expected bool
}

func TestBitVector(t *testing.T) {
	initialNumberOfBits := wordSize*10 + 1

	bv := newBitVector(initialNumberOfBits)

	if bv.len() != initialNumberOfBits {
		t.Fatalf("expected len to be %d but got %d", initialNumberOfBits, bv.len())
	}

	for i := 0; i < initialNumberOfBits; i++ {
		bv.setBit(i, true)
	}

	bv.setBit(3, false)
	bv.setBit(11, false)
	bv.setBit(13, false)
	bv.setBit(23, false)
	bv.setBit(24, false)
	bv.setBit(25, false)
	bv.setBit(42, false)
	bv.setBit(63, false)
	bv.setBit(64, false)
	bv.setBit(255, false)
	bv.setBit(256, false)

	getBitTestCases := []GetBitTestCase{
		{0, true},
		{1, true},
		{2, true},
		{3, false},
		{4, true},
		{7, true},
		{8, true},
		{9, true},
		{10, true},
		{11, false},
		{12, true},
		{13, false},
		{23, false},
		{24, false},
		{25, false},
		{42, false},
		{15, true},
		{16, true},
		{62, true},
		{63, false},
		{64, false},
		// Test past the first word
		{65, true},
		{72, true},
		{79, true},
		{80, true},
		{255, false},
		{256, false},
		{511, true},
		{512, true},
	}

	for _, v := range getBitTestCases {
		actual := bv.getBit(v.position)
		if actual != v.expected {
			t.Fatalf("expected %dth bit to be %t but got %t", v.position, v.expected, actual)
		}
	}
}

func TestBitVectorBoundPanic_GetBit_Lower(t *testing.T) {
	defer func() { _ = recover() }()

	initialNumberOfBits := wordSize*10 + 1
	bv := newBitVector(initialNumberOfBits)
	bv.getBit(-1)

	t.Fatalf("expected get bit lower bound panic")
}

func TestBitVectorBoundPanic_GetBit_Upper(t *testing.T) {
	defer func() { _ = recover() }()
	initialNumberOfBits := wordSize*10 + 1
	bv := newBitVector(initialNumberOfBits)
	bv.getBit(initialNumberOfBits)

	t.Fatalf("expected get bit upper bound panic")
}

func TestBitVectorBoundPanic_SetBit_Lower(t *testing.T) {
	defer func() {
		if r := recover(); r != nil {
			return
		}
		t.Fatalf("expected set bit lower bound panic")
	}()
	initialNumberOfBits := wordSize*10 + 1
	bv := newBitVector(initialNumberOfBits)
	bv.setBit(-1, true)
}

func TestBitVectorBoundPanic_SetBit_Upper(t *testing.T) {
	defer func() {
		if r := recover(); r != nil {
			return
		}
		t.Fatalf("expected set bit upper bound panic")
	}()
	initialNumberOfBits := wordSize*10 + 1
	bv := newBitVector(initialNumberOfBits)
	bv.setBit(initialNumberOfBits, true)
}
