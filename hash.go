package main

import (
	"encoding/hex"
	"strings"

	"lukechampine.com/blake3"
)

// boothLeastRotation gets the least rotation of a circular string.
// https://en.wikipedia.org/wiki/Lexicographically_minimal_string_rotation
// this is generally over commented but I'm keeping it this way for now. - Tim
func boothLeastRotation(sequence string) int {

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

func seqhash(sequence string, circular bool) string {
	s := strings.ToUpper(sequence) // why to upper here?
	if circular {
		r := boothLeastRotation(s)
		concatSeq := s + s
		s = concatSeq[r : r+len(s)]
	}
	b := blake3.Sum256([]byte(s))
	return hex.EncodeToString(b[:])
}
