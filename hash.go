package poly

import (
	"encoding/hex"
	"errors"
	"sort"
	"strings"

	"lukechampine.com/blake3"
)

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

/******************************************************************************
Dec, 2, 2020

Seqhash stuff starts here.

There is a big problem with current sequence databases - they all use different
identifiers and accession numbers. This means cross-referencing databases is
a complicated exercise, especially as the quantity of databases increases, or if
you need to compare "wild" DNA sequences.

Seqhash is a simple algorithm to produce consistent identifiers for any genetic sequence. The
basic premise of the Seqhash algorithm is to hash sequences with the hash being a robust
cross-database identifier. Sequences themselves shouldn't be used as a database index
(often, they're too big), so a hash based off of a sequence is the next best thing.

Usability wise, you should be able to Seqhash any rotation of a sequence in any direction and
get a consistent hash.

The Seqhash algorithm makes several opinionated design choices, primarily to make working
with Seqhashes more consistent and nice. The Seqhash algorithm only uses a single hash function,
Blake3, and only operates on DNA, RNA, and Protein sequences. These identifiers will be seen
by human beings, so versioning and metadata is attached to the front of the hashes so that
a human operator can quickly identify problems with hashing.

If the sequence is DNA or RNA, the Seqhash algorithm needs to know whether or not the nucleic
acid is circular and/or double stranded. If circular, the sequence is rotated to a deterministic
point. If double stranded, the sequence is compared to its reverse complement, and the lexiographically
minimal sequence is taken (whether or not the min or max is used doesn't matter, just needs to
be consistent).

If the sequence is RNA, the sequence will be converted to DNA before hashing. While the full Seqhash
will still be different between RNA and DNA (due to the metadata string), the hash afterwards will be the same.
This makes it easy to cross reference DNA and RNA sequences. This fact is important for parts of Poly
store that relate to storing and searching large quantities of sequences - deduplication can easily
be used on those Seqhashes to save a lot of space.

For DNA or RNA sequences, only ATUGCYRSWKMBDHVNZ characters are allowed. For Proteins,
only ACDEFGHIKLMNPQRSTVWYUO*BXZ characters are allowed in sequences. Selenocysteine (Sec; U) and pyrrolysine
(Pyl; O) are included in the protein character set - usually U and O don't occur within protein sequences,
but for certain organisms they do, and it is certainly a relevant amino acid for those particular proteins.

A Seqhash is separated into 3 different elements divided by underscores. It looks like the following:

v1_DCD_4b0616d1b3fc632e42d78521deb38b44fba95cca9fde159e01cd567fa996ceb9

The first element is the version tag (v1 for version 1). If there is ever a Seqhash version 2, this tag
will differentiate seqhashes. The second element is the metadata tag, which has 3 letters. The first letter
codes for the sequenceType (D for DNA, R for RNA, and P for Protein). The second letter codes for whether or
not the sequence is circular (C for Circular, L for Linear). The final letter codes for whether or not the
sequence is double stranded (D for Double stranded, S for Single stranded). The final element is the blake3
hash of the sequence (once rotated and complemented, as stated above).

Seqhash is a simple algorithm that allows for much better indexing of genetic sequences than what is
currently available. I hope it will be widely adopted someday.

En Taro Adun,
Keoni

******************************************************************************/

// Seqhash is a function to create Seqhashes, a specific kind of identifier.
func Hash(sequence string, sequenceType string, circular bool, doubleStranded bool) (string, error) {
	// By definition, Seqhashes are of uppercase sequences
	sequence = strings.ToUpper(sequence)
	// If RNA, convert to a DNA sequence. The hash itself between a DNA and RNA sequence will not
	// be different, but their Seqhash will have a different metadata string (R vs D)
	if sequenceType == "RNA" {
		sequence = strings.ReplaceAll(sequence, "U", "T")
	}

	// Run checks on the input
	if sequenceType != "DNA" && sequenceType != "RNA" && sequenceType != "PROTEIN" {
		return "", errors.New("Only sequenceTypes of DNA, RNA, or PROTEIN allowed. Got sequenceType: " + sequenceType)
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
			// The release notes https://web.expasy.org/docs/relnotes/relstat.html
			// also state there are Asx (B), Glx (Z), and Xaa (X) amino acids, so
			// these are added in as well.
			if !strings.Contains("ACDEFGHIKLMNPQRSTVWYUO*BXZ", string(char)) {
				return "", errors.New("Only letters ACDEFGHIKLMNPQRSTVWYUO*BXZ are allowed for Proteins. Got letter: " + string(char))
			}
		}
	}
	// There is no check for circular proteins since proteins can be circular
	if sequenceType == "PROTEIN" && doubleStranded {
		return "", errors.New("Proteins cannot be double stranded")
	}

	// Gets Deterministic sequence based off of metadata + sequence
	var deterministicSequence string
	switch {
	case circular && doubleStranded:
		potentialSequences := []string{RotateSequence(sequence), RotateSequence(ReverseComplement(sequence))}
		sort.Strings(potentialSequences)
		deterministicSequence = potentialSequences[0]
	case circular && !doubleStranded:
		deterministicSequence = RotateSequence(sequence)
	case !circular && doubleStranded:
		potentialSequences := []string{sequence, ReverseComplement(sequence)}
		sort.Strings(potentialSequences)
		deterministicSequence = potentialSequences[0]
	case !circular && !doubleStranded:
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
	// Get 2nd letter. C for circular, L for Linear
	if circular {
		circularLetter = "C"
	} else {
		circularLetter = "L"
	}
	// Get 3rd letter. D for Double stranded, S for Single stranded
	if doubleStranded {
		doubleStrandedLetter = "D"
	} else {
		doubleStrandedLetter = "S"
	}

	newhash := blake3.Sum256([]byte(deterministicSequence))
	seqhash := "v1" + "_" + sequenceTypeLetter + circularLetter + doubleStrandedLetter + "_" + hex.EncodeToString(newhash[:])
	return seqhash, nil

}

// Hash is a method wrapper for hashing Sequence structs. Note that
// all sequence structs are, by default, double-stranded sequences,
// since Genbank does not track whether or not a given sequence in their
// database is single stranded or double stranded.
func (sequence Sequence) Hash() (string, error) {
	if sequence.Meta.Locus.MoleculeType == "" {
		return "", errors.New("No MoleculeType found for sequence")
	}
	var sequenceType string
	if strings.Contains(sequence.Meta.Locus.MoleculeType, "DNA") {
		sequenceType = "DNA"
	}
	if strings.Contains(sequence.Meta.Locus.MoleculeType, "RNA") {
		sequenceType = "RNA"
	}
	if (sequenceType != "DNA") && (sequenceType != "RNA") {
		return "", errors.New("SequenceType not found. Looking for MoleculeTypes with DNA or RNA, got: " + sequence.Meta.Locus.MoleculeType)
	}
	// If not explicitly circular, assume linear. All sequences are by default doubleStranded
	newSeqhash, err := Hash(sequence.Sequence, sequenceType, sequence.Meta.Locus.Circular, true)
	if err != nil {
		return "", err
	}
	return newSeqhash, nil

}
