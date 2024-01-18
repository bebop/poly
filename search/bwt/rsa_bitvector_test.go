package bwt

import (
	"testing"
)

type rsaRankTestCase struct {
	val          bool
	bitPosition  int
	expectedRank int
}

func TestRSARank_singlePartialChunk(t *testing.T) {
	if wordSize != 64 {
		t.Skip()
	}

	bitsToTruncate := 22
	initialNumberOfBits := wordSize*2 - bitsToTruncate

	rsa := newTestRSAFromWords(initialNumberOfBits,
		0xffffffff00000000,
		0x00000000ffc00000,
	)

	testCases := []rsaRankTestCase{
		{true, 0, 0}, {false, 0, 0},

		{true, 64, 32}, {false, 64, 32},

		{true, 96, 32}, {false, 96, 64},

		{true, 105, 41}, {false, 105, 64},
	}

	for _, tc := range testCases {
		rank := rsa.Rank(tc.val, tc.bitPosition)
		if rank != tc.expectedRank {
			t.Fatalf("expected rank(%t, %d) to be %d but got %d", tc.val, tc.bitPosition, tc.expectedRank, rank)
		}
	}
}

func TestRSARank_singleCompleteChunk_PastBounds_Ones(t *testing.T) {
	rsa := newTestRSAFromWords(64*4,
		0x0000000000000000,
		0xffffffffffffffff,
		0x0000000000000000,
		0xffffffffffffffff,
	)

	testCases := []rsaRankTestCase{
		{true, 0, 0}, {false, 0, 0},
		{true, 255, 127}, {false, 255, 128},
		{true, 256, 128}, {false, 256, 128},
	}

	for _, tc := range testCases {
		rank := rsa.Rank(tc.val, tc.bitPosition)
		if rank != tc.expectedRank {
			t.Fatalf("expected rank(%t, %d) to be %d but got %d", tc.val, tc.bitPosition, tc.expectedRank, rank)
		}
	}
}

func TestRSARank_singleCompleteChunk_PastBounds_Zeros(t *testing.T) {
	rsa := newTestRSAFromWords(64*4,
		0xffffffffffffffff,
		0x0000000000000000,
		0xffffffffffffffff,
		0x0000000000000000,
	)

	testCases := []rsaRankTestCase{
		{true, 0, 0}, {false, 0, 0},
		{true, 255, 128}, {false, 255, 127},
		{true, 256, 128}, {false, 256, 128},
	}

	for _, tc := range testCases {
		rank := rsa.Rank(tc.val, tc.bitPosition)
		if rank != tc.expectedRank {
			t.Fatalf("expected rank(%t, %d) to be %d but got %d", tc.val, tc.bitPosition, tc.expectedRank, rank)
		}
	}
}

func TestRSARank_singleCompleteChunk(t *testing.T) {
	initialNumberOfBits := wordSize * 4

	rsa := newTestRSAFromWords(initialNumberOfBits,
		0x8000000000000001,
		0xff0f30fffacea80d,
		0x90e0a0e0b0e0cf0c,
		0x3d0f064f7206f717,
	)

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
		rank := rsa.Rank(tc.val, tc.bitPosition)
		if rank != tc.expectedRank {
			t.Fatalf("expected rank(%t, %d) to be %d but got %d", tc.val, tc.bitPosition, tc.expectedRank, rank)
		}
	}
}

func TestRSARank_multipleChunks(t *testing.T) {
	rsa := newTestRSAFromWords((8*4+3)*64,
		0x0000000000000000,
		0xffffffffffffffff,
		0x0000000000000000,
		0xffffffffffffffff,

		0xffffffffffffffff,
		0x0000000000000000,
		0xffffffffffffffff,
		0x0000000000000000,

		0x0000000000000000,
		0xffffffffffffffff,
		0x0000000000000000,
		0xffffffffffffffff,

		0xffffffffffffffff,
		0x0000000000000000,
		0xffffffffffffffff,
		0x0000000000000000,

		// If Jacobson rank is still there, this should go past the first
		// chunk
		0xffffffffffffffff,
		0x0000000000000000,
		0xffffffffffffffff,
		0x0000000000000000,

		0xffffffffffffffff,
		0x0000000000000000,
		0xffffffffffffffff,
		0x0000000000000000,

		0xffffffffffffffff,
		0x0000000000000000,
		0xffffffffffffffff,
		0x0000000000000000,

		0xffffffffffffffff,
		0x0000000000000000,
		0xffffffffffffffff,
		0x0000000000000000,

		// If Jacobson rank is still there, this should go past the second
		// chunk
		0xffffffffffffffff,
		0x0000000000000000,
		0xffffffffffffffff,
	)

	testCases := []rsaRankTestCase{
		{true, 0, 0}, {false, 0, 0},

		{true, 64, 0}, {false, 64, 64},
		{true, 128, 64}, {false, 128, 64},
		{true, 192, 64}, {false, 192, 128},
		{true, 256, 128}, {false, 256, 128},

		{true, 320, 192}, {false, 256, 128},
		{true, 384, 192}, {false, 384, 192},
		{true, 448, 256}, {false, 448, 192},
		{true, 512, 256}, {false, 512, 256},

		{true, 576, 256}, {false, 576, 320},
		{true, 640, 320}, {false, 640, 320},
		{true, 704, 320}, {false, 704, 384},
		{true, 768, 384}, {false, 768, 384},

		{true, 832, 448}, {false, 832, 384},
		{true, 896, 448}, {false, 896, 448},

		{true, 1024, 512}, {false, 1024, 512},

		{true, 2048, 1024}, {false, 2048, 1024},
	}

	for _, tc := range testCases {
		rank := rsa.Rank(tc.val, tc.bitPosition)
		if rank != tc.expectedRank {
			t.Fatalf("expected rank(%t, %d) to be %d but got %d", tc.val, tc.bitPosition, tc.expectedRank, rank)
		}
	}
}

type rsaSelectTestCase struct {
	val              bool
	rank             int
	expectedPosition int
}

func TestRSASelect(t *testing.T) {
	bitsToTruncate := 17
	initialNumberOfBits := wordSize*4 - bitsToTruncate
	rsa := newTestRSAFromWords(initialNumberOfBits,
		0x8010000000010000, // 1Count = 3
		0xfff1ffffffffffff, // 1Count = 63
		0x0000010000000000, // 1Count = 1
		0xffffffffffffffff, // Possible 1Count = 47
	)

	testCases := []rsaSelectTestCase{
		{true, 0, 0},
		{true, 1, 11},
		{true, 2, 47},
		{false, 0, 1},
		{false, 1, 2},
		{false, 3, 4},
		{false, 8, 9},
		{false, 9, 10},
		{false, 10, 12},
		{false, 11, 13},
		{false, 60, 63},

		{true, 3, 64},
		{true, 9, 70},
		{true, 13, 74},
		{true, 14, 75},
		{true, 15, 79},
		{true, 16, 80},
		{true, 63, 127},
		{false, 61, 76},
		{false, 62, 77},
		{false, 63, 78},

		{true, 64, 151},
		{true, 65, 192},
		{true, 111, 238},
		{false, 64, 128},

		{false, 126, 191},

		// Select of penultimate ranks should be the positions at which they appear.
		{true, 111, rsa.bv.len() - 1},
		{false, 126, 191},

		// Max bitvector positions for the max rank should be at the ends of the bitvector
		{true, 112, rsa.bv.len()},
		{false, 127, rsa.bv.len()},
	}

	for _, tc := range testCases {
		position, ok := rsa.Select(tc.val, tc.rank)

		if !ok {
			t.Fatalf("expected select(%t, %d) to be %d but went out of range", tc.val, tc.rank, tc.expectedPosition)
		}

		if position != tc.expectedPosition {
			t.Fatalf("expected select(%t, %d) to be %d but got %d", tc.val, tc.rank, tc.expectedPosition, position)
		}
	}
}

func TestRSASelect_notOk(t *testing.T) {
	bitsToTruncate := 17
	initialNumberOfBits := wordSize*4 - bitsToTruncate
	rsa := newTestRSAFromWords(initialNumberOfBits,
		0x8010000000010000,
		0xfff1ffffffffffff,
		0x0000010000000000,
		0xffffffffffffffff,
	)

	if _, ok := rsa.Select(true, -1); ok {
		t.Fatalf("expected select(true, -1) to be not ok but somehow returned a value")
	}

	pos, ok := rsa.Select(true, 111)
	if !ok {
		t.Fatalf("expected select(true, 111) to be ok but somehow got not ok")
	}

	if pos != 238 {
		t.Fatalf("expected select(true, 111) to be 238 but got %d", pos)
	}

	if _, ok := rsa.Select(true, 239); ok {
		t.Fatalf("expected select(true, 239) to be not ok but somehow returned a value")
	}
}

func newTestRSAFromWords(sizeInBits int, wordsToCopy ...uint64) rsaBitVector {
	bv := newBitVector(sizeInBits)
	for i := 0; i < sizeInBits; i++ {
		w := wordsToCopy[i/64]
		mask := uint64(1) << uint64(63-i%64)
		bit := w&mask != 0
		bv.setBit(i, bit)
	}
	return newRSABitVectorFromBitVector(bv)
}
