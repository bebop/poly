package poly

import (
	_ "crypto/md5"
	_ "crypto/sha1"
	_ "crypto/sha256"
	_ "crypto/sha512"
	"encoding/hex"
	"hash"
	"io"
	"strings"

	_ "golang.org/x/crypto/blake2b"
	_ "golang.org/x/crypto/blake2s"
	_ "golang.org/x/crypto/ripemd160"
	_ "golang.org/x/crypto/sha3"
)

// Where each hash function comes from.
// MD5                         // import crypto/md5
// SHA1                        // import crypto/sha1
// SHA224                      // import crypto/sha256
// SHA256                      // import crypto/sha256
// SHA384                      // import crypto/sha512
// SHA512                      // import crypto/sha512
// MD5SHA1                     // no implementation; MD5+SHA1 used for TLS RSA
// RIPEMD160                   // import golang.org/x/crypto/ripemd160
// SHA3_224                    // import golang.org/x/crypto/sha3
// SHA3_256                    // import golang.org/x/crypto/sha3
// SHA3_384                    // import golang.org/x/crypto/sha3
// SHA3_512                    // import golang.org/x/crypto/sha3
// SHA512_224                  // import crypto/sha512
// SHA512_256                  // import crypto/sha512
// BLAKE2s_256                 // import golang.org/x/crypto/blake2s
// BLAKE2b_256                 // import golang.org/x/crypto/blake2b
// BLAKE2b_384                 // import golang.org/x/crypto/blake2b
// BLAKE2b_512                 // import golang.org/x/crypto/blake2b

// boothLeastRotation gets the least rotation of a circular string.
func boothLeastRotation(sequence string) int {

	// https://en.wikipedia.org/wiki/Lexicographically_minimal_string_rotation
	// this is generally over commented but I'm keeping it this way for now. - Tim

	// first concatenate the sequence to itself to avoid modular arithmateic
	sequence += sequence // maybe do this as a buffer just for speed? May get annoying with larger sequences.
	leastRotationIndex := 0

	//initializing failure slice.
	failureSlice := make([]int, len(sequence))
	for i := range failureSlice {
		failureSlice[i] = -1
	}
	// iterate through each character in the doubled over sequence
	for characterIndex := 1; characterIndex < len(sequence); characterIndex++ {
		// get character
		character := sequence[characterIndex]
		// get failure
		failure := failureSlice[characterIndex-leastRotationIndex-1]
		// while failure does not equal -1 and character does not equal the character found at the least rotation + failure + 1 <- why this?
		for failure != -1 && character != sequence[leastRotationIndex+failure+1] {

			// if character is lexically less than whatever is at the leastRotationIndex index update leastRotation index
			if character < sequence[leastRotationIndex+failure+1] {
				leastRotationIndex = characterIndex - failure - 1
			}
			// update failure using previous failure as index?
			failure = failureSlice[failure]
		}

		// if character does not equal whatever character is at leastRotationIndex plus failure.
		if character != sequence[leastRotationIndex+failure+1] {

			// if character is lexically less then what is rotated least leastRotatationIndex gets value of character index.
			if character < sequence[leastRotationIndex] {
				leastRotationIndex = characterIndex
			}
			// assign -1 to whatever is at the index of difference between character and rotation indeces.
			failureSlice[characterIndex-leastRotationIndex] = -1

			// if character does equal whatever character is at leastRotationIndex plus failure.
		} else {
			// assign failure + 1 at the index of difference between character and rotation indeces.
			failureSlice[characterIndex-leastRotationIndex] = failure + 1
		}
	} // end loop

	return leastRotationIndex
}

// RotateSequence rotates circular sequences to deterministic point.
func RotateSequence(sequence string) string {
	rotationIndex := boothLeastRotation(sequence)
	var sequenceBuilder strings.Builder

	// writing the same sequence twice. using build incase of very long circular genome.
	sequenceBuilder.WriteString(sequence)
	sequenceBuilder.WriteString(sequence)

	concatenatedSequence := sequenceBuilder.String()
	sequence = concatenatedSequence[rotationIndex : rotationIndex+len(sequence)]
	return sequence
}

// Hash is a method wrapper for hashing Sequence structs.
func (sequence Sequence) Hash(hash hash.Hash) string {
	if sequence.Meta.Locus.Circular {
		sequence.Sequence = RotateSequence(sequence.Sequence)
	}
	seqHash, _ := hashSequence(sequence.Sequence, hash)
	return seqHash
}

// Hash is a method wrapper for hashing sequences contained in Feature structs.
func (feature Feature) Hash(hash hash.Hash) string {
	seqHash, _ := hashSequence(feature.GetSequence(), hash)
	return seqHash
}

// hashSequence takes a string and a hashing function and returns a hashed string.
func hashSequence(sequence string, hash hash.Hash) (string, error) {
	io.WriteString(hash, strings.ToUpper(sequence))
	return hex.EncodeToString(hash.Sum(nil)), nil
}
