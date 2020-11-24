package poly

import (
	_ "crypto/md5"
	_ "crypto/sha1"
	_ "crypto/sha256"
	_ "crypto/sha512"
	"encoding/hex"
	"errors"
	"hash"
	"io"
	"sort"
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

var sequenceCharMap = map[rune]rune{'A': 'T', 'T': 'A', 'U': 'A', 'G': 'C', 'C': 'G', 'Y': 'R', 'R': 'Y', 'S': 'S', 'W': 'W', 'K': 'M', 'M': 'K', 'B': 'V', 'D': 'H', 'H': 'D', 'V': 'B', 'N': 'N', 'Z': 'Z'}

func ReverseComplement(sequence string) string {
	n := len(sequence)
	complement := make([]rune, n)
	reverseComplement := make([]rune, n)
	for _, char := range strings.ToUpper(sequence) {
		complement = append(complement, sequenceCharMap[char])
	}
	for _, char := range complement {
		n--
		reverseComplement[n] = char
	}
	return string(reverseComplement[n:])
}

// SeqHash is a function to create SeqHashes, a specific kind of identifier
func SeqHash(sequence string, sequenceType string, circular bool, doubleStranded bool) (string, error) {
	// By definition, SeqHashes are of uppercase sequences
	sequence = strings.ToUpper(sequence)

	// Run checks on the input
	if sequenceType != "DNA" || sequenceType != "RNA" || sequenceType != "PROTEIN" {
		return "", errors.New("Only sequenceTypes of DNA, RNA, or PROTEIN allowed")
	}
	if sequenceType == "DNA" || sequenceType == "RNA" {
		for _, char := range sequence {
			if !strings.Contains("ATUGCYRSWKMBDHVNZ", string(char)) {
				return "", errors.New("Only letters ATUGCYRSWKMBDHVNZ are allowed for DNA/RNA. Got letter: " + string(char))
			}
		}
	}
	if sequenceType == "PROTEIN" {
		for _, char := range sequence {
			// Selenocysteine (Sec; U) and pyrrolysine (Pyl; O) are added
			// in accordance with https://www.uniprot.org/help/sequences
			if !strings.Contains("ACDEFGHIKLMNPQRSTVWYUO", string(char)) {
				return "", errors.New("Only letters ACDEFGHIKLMNPQRSTVWYUO are allowed for Proteins. Got letter: " + string(char))
			}
		}
	}
	if sequenceType == "PROTEIN" && doubleStranded == true {
		return "", errors.New("Proteins cannot be double stranded")
	}

	// Gets Deterministic sequence based off of metadata + sequence
	var deterministicSequence string
	switch {
	case circular == true && doubleStranded == true:
		potentialSequences := []string{RotateSequence(sequence), RotateSequence(ReverseComplement(sequence))}
		sort.Strings(potentialSequences)
		deterministicSequence = potentialSequences[0]
	case circular == true && doubleStranded == false:
		deterministicSequence = RotateSequence(sequence)
	case circular == false && doubleStranded == true:
		potentialSequences := []string{sequence, ReverseComplement(sequence)}
		sort.Strings(potentialSequences)
		deterministicSequence = potentialSequences[0]
	case circular == false && doubleStranded == false:
		deterministicSequence = sequence
	}

	// Build 3 letter metadata
	var sequenceTypeLetter string
	var circularLetter string
	var doubleStrandedLetter string
	// Get first letter. D for DNA, R for RNA, and P for Protein
	switch sequenceType {
	case "DNA":
		sequenceTypeLetter = "D"
	case "RNA":
		sequenceTypeLetter = "R"
	case "PROTEIN":
		sequenceTypeLetter = "P"
	}
	// Get 2nd letter. C for circular, L for Liner
	if circular == true {
		circularLetter = "C"
	} else {
		circularLetter = "L"
	}
	// Get 3rd letter. D for Double stranded, S for Single stranded
	if doubleStranded == true {
		doubleStrandedLetter = "D"
	} else {
		doubleStrandedLetter = "S"
	}

	newhash := blake3.Sum256([]byte(deterministicSequence))
	seqhash := "v1" + "_" + sequenceTypeLetter + circularLetter + doubleStrandedLetter + "_" + hex.EncodeToString(newhash[:])
	return seqhash, nil

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
