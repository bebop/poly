package main

import (
	"crypto"
	_ "crypto/md5"
	_ "crypto/sha1"
	_ "crypto/sha256"
	_ "crypto/sha512"
	"encoding/hex"
	"errors"
	"io"
	"strings"

	_ "golang.org/x/crypto/blake2b"
	_ "golang.org/x/crypto/blake2s"
	_ "golang.org/x/crypto/ripemd160"
	_ "golang.org/x/crypto/sha3"
	"lukechampine.com/blake3"
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

// BoothLeastRotation gets the least rotation of a circular string.
// https://en.wikipedia.org/wiki/Lexicographically_minimal_string_rotation
// this is generally over commented but I'm keeping it this way for now. - Tim
func BoothLeastRotation(sequence string) int {

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

//RotateSequence rotates circular sequences to deterministic point.
func RotateSequence(sequence string) string {
	rotationIndex := BoothLeastRotation(sequence)
	concatenatedSequence := sequence + sequence
	sequence = concatenatedSequence[rotationIndex : rotationIndex+len(sequence)]
	return sequence
}

// GenericSequenceHash takes an AnnotatedSequence and a hash function and hashes it.
// from https://stackoverflow.com/questions/32620290/how-to-dynamically-switch-between-hash-algorithms-in-golang <- this had a bug I had to fix! - Tim
func GenericSequenceHash(annotatedSequence AnnotatedSequence, hash crypto.Hash) (string, error) {
	if !hash.Available() {
		return "", errors.New("hash unavailable")
	}
	if annotatedSequence.Meta.Locus.Circular {
		annotatedSequence.Sequence.Sequence = RotateSequence(annotatedSequence.Sequence.Sequence)
	}
	h := hash.New()
	io.WriteString(h, strings.ToUpper(annotatedSequence.Sequence.Sequence))
	return hex.EncodeToString(h.Sum(nil)), nil
}

// Hash is a method wrapper for hashing annotatedSequence structs.
func (annotatedSequence AnnotatedSequence) Hash(hash crypto.Hash) string {
	seqHash, _ := GenericSequenceHash(annotatedSequence, hash)
	return seqHash
}

// Blake3SequenceHash Blake3 function doesn't use standard golang hash interface
// so we couldn't use it with the generic sequence hash.
func Blake3SequenceHash(annotatedSequence AnnotatedSequence) string {

	if annotatedSequence.Meta.Locus.Circular {
		annotatedSequence.Sequence.Sequence = RotateSequence(annotatedSequence.Sequence.Sequence)
	}

	b := blake3.Sum256([]byte(strings.ToUpper(annotatedSequence.Sequence.Sequence)))
	return hex.EncodeToString(b[:])
}

// Blake3Hash is a method wrapper for hashing annotatedSequence structs with Blake3.
func (annotatedSequence AnnotatedSequence) Blake3Hash() string {
	seqHash := Blake3SequenceHash(annotatedSequence)
	return seqHash
}
