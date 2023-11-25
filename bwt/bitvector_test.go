package bwt

import (
	"fmt"
	"math/bits"
	"strconv"
	"testing"
)

type GetBitTestCase struct {
	position int
	expected bool
}

func TestBitVector(t *testing.T) {
	initialNumberOfBits := wordSize*10 + 1
	expectedCapacity := 11

	bv := newBitVector(initialNumberOfBits)

	if bv.capacity() != expectedCapacity {
		t.Fatalf("expected capacity to be %d but got %d", expectedCapacity, bv.capacity())
	}

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
		{72, true},
		{79, true},
		{80, true},
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

func TestBitVectorPush_NextPushLessThanCapacity_Single(t *testing.T) {
	initialNumberOfBits := wordSize*10 + 1
	bv := newBitVector(initialNumberOfBits)
	bv.push(true)

	expectedCapacity := 11
	if bv.capacity() != expectedCapacity {
		t.Fatalf("expected capacity to be %d but got %d", expectedCapacity, bv.capacity())
	}

	expectedLength := initialNumberOfBits + 1
	if bv.len() != expectedLength {
		t.Fatalf("expected len to be %d but got %d", expectedLength, bv.len())
	}

	if bv.getBit(initialNumberOfBits) != true {
		t.Fatalf("expected %dth bit to be %t but got %t", initialNumberOfBits, true, bv.getBit(initialNumberOfBits))
	}
}

func TestBitVectorPush_NextPushGreaterThanCapacity_Single(t *testing.T) {
	initialNumberOfBits := wordSize * 10
	bv := newBitVector(initialNumberOfBits)
	initialCapacity := bv.capacity()
	bv.push(true)

	if bv.capacity() <= initialCapacity {
		t.Fatalf("expected capacity to have grown. currently the capacity is %d and was previously %d", bv.capacity(), initialCapacity)
	}

	expectedLength := initialNumberOfBits + 1
	if bv.len() != expectedLength {
		t.Fatalf("expected len to be %d but got %d", expectedLength, bv.len())
	}

	if bv.getBit(initialNumberOfBits) != true {
		t.Fatalf("expected %dth bit to be %t but got %t", initialNumberOfBits, true, bv.getBit(initialNumberOfBits))
	}
}

type rsaRankTestCase struct {
	val          bool
	bitPosition  int
	expectedRank int
}

func TestRSARank_wordSize64_singleCompleteChunk(t *testing.T) {
	if wordSize != 64 {
		t.Skip()
	}

	initialNumberOfBits := wordSize * 4
	bv := newBitVector(initialNumberOfBits)

	w0 := uint(0x8000000000000001)
	fmt.Println(bits.OnesCount(^w0))
	fmt.Println(strconv.FormatUint(uint64(w0), 2))
	w1 := uint(0xff0f30fffacea80d)
	fmt.Println(bits.OnesCount(^w1))
	fmt.Println(strconv.FormatUint(uint64(w1), 2))
	w2 := uint(0x90e0a0e0b0e0cf0c)
	fmt.Println(bits.OnesCount(^w2))
	fmt.Println(strconv.FormatUint(uint64(w2), 2))
	w3 := uint(0x3d0f064f7206f717)
	fmt.Println(bits.OnesCount(^w3))
	fmt.Println(strconv.FormatUint(uint64(w3), 2))

	bv.bits = []uint{w0, w1, w2, w3}

	rsa := newRSABitVector(bv)

	testCases := []rsaRankTestCase{
		{true, 0, 0}, {false, 0, 0},
		{true, 1, 1}, {false, 1, 0},
		{true, 2, 1}, {false, 2, 1},
		{true, 3, 1}, {false, 3, 2},
		{true, 62, 1}, {false, 62, 61},
		{true, 63, 1}, {false, 63, 62},

		{true, 64, 2}, {false, 64, 62},
		{true, 65, 3}, {false, 65, 62},
		{true, 72, 10}, {false, 72, 62},
		{true, 127, 40}, {false, 127, 87},

		{true, 128, 41}, {false, 128, 87},
		{true, 129, 42}, {false, 129, 87},
		{true, 130, 42}, {false, 130, 88},
		{true, 131, 42}, {false, 131, 89},
		{true, 132, 43}, {false, 132, 89},
		{true, 133, 43}, {false, 133, 90},
		{true, 159, 51}, {false, 159, 108},
		{true, 160, 51}, {false, 160, 109},
		{true, 161, 52}, {false, 161, 109},
		{true, 162, 52}, {false, 162, 110},
		{true, 163, 53}, {false, 163, 110},
		{true, 164, 54}, {false, 164, 110},
		{true, 165, 54}, {false, 165, 111},
		{true, 176, 57}, {false, 176, 119},
		{true, 177, 58}, {false, 177, 119},
		{true, 178, 59}, {false, 178, 119},
		{true, 179, 59}, {false, 179, 120},
		{true, 180, 59}, {false, 180, 121},
		{true, 183, 62}, {false, 183, 121},
		{true, 184, 63}, {false, 184, 121},
		{true, 185, 63}, {false, 185, 122},
		{true, 186, 63}, {false, 186, 123},
		{true, 187, 63}, {false, 187, 124},
		{true, 188, 63}, {false, 188, 125},
		{true, 189, 64}, {false, 189, 125},
		{true, 190, 65}, {false, 190, 125},
		{true, 191, 65}, {false, 191, 126},

		{true, 192, 65}, {false, 192, 127},
		{true, 193, 65}, {false, 193, 128},
		{true, 194, 65}, {false, 194, 129},
		{true, 195, 66}, {false, 195, 129},
		{true, 196, 67}, {false, 196, 129},
		{true, 248, 94}, {false, 248, 154},
		{true, 249, 94}, {false, 249, 155},
		{true, 250, 94}, {false, 250, 156},
		{true, 251, 94}, {false, 251, 157},
		{true, 252, 95}, {false, 252, 157},
		{true, 253, 95}, {false, 253, 158},
		{true, 254, 96}, {false, 254, 158},
		{true, 255, 97}, {false, 255, 158},
	}

	for _, tc := range testCases {
		rank := rsa.rank(tc.val, tc.bitPosition)
		if rank != tc.expectedRank {
			t.Fatalf("expected rank(%t, %d) to be %d but got %d", tc.val, tc.bitPosition, tc.expectedRank, rank)
		}
	}

}
