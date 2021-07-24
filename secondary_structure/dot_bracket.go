package secondary_structure

import (
	"fmt"
	"strings"

	"github.com/TimothyStiles/poly/checks"
)

/******************************************************************************
July, 01, 2021

This file defines the structs and functions needed to get the 'annotated'
structure and `SecondaryStructure` of a RNA from its 'dot-bracket' notation.

'Dot-bracket' notation of a secondary structure is a string where each
character represents a base. Unpaired nucleotides are represented with a '.'
and base pairs are represented by paranthesis. '(' denotes the opening base
and ')' denotes the closing base of a base pair.
For example,
	dot-bracket structure: . . ( ( . . ) ) . .
									index: 0 1 2 3 4 5 6 7 8 9
denotes a hairpin where the bases at index 2 and 7,
and at index 3 and 6 are paired.

'Annotated' structure of a secondary structure is a string where each
character represents a base. In this notation, unpaired nucleotides are set
to a character depending on what part of the secondary structure is in.
The character mapping for the unpaired nucleotides is:
* On the exterior part of a RNA structure (not part of any loop): 'e'
* On the single stranded region of a hairpin: 'h'
* On a single stranded region of a multiloop: 'm'
* Part of an interior loop in a stem of a hairpin or multiloop: 'i'
For example,
	dot-bracket structure: . . ( ( . . ) ) . .
		annotated structure: e e ( ( h h ) ) e e
									index: 0 1 2 3 4 5 6 7 8 9
Note that the paranthesis surrounding one or more 'h' or 'm' characters form
the stem of that hairpin or multiloop respectively. Thus, in the example above,
the bases at index 2, 3, 6, and 7 form the stem of the hairpin loop enclosed
by the bases at index 3 and 6.

 ******************************************************************************/

const (
	dotBracketUnpairedNucleotide       byte = '.'
	dotBracketStemFivePrimeNucleotide  byte = '('
	dotBracketStemThreePrimeNucleotide byte = ')'

	exteriorLoopUnpairedNucleotide          byte = 'e'
	interiorLoopUnpairedNucleotide          byte = 'i'
	hairpinLoopNucleotide                   byte = 'h'
	multiLoopSingleStrandedRegionNucleotide byte = 'm'
)

// parseCompound holds all information needed to compute the annotated structure
// and `SecondaryStructure` of a RNA sequence.
type parseCompound struct {
	length             int   // length of `sequence`
	pairTable          []int // (see `pairTable()`)
	annotatedStructure []byte
}

// SecondaryStructureFromDotBracket returns the annotated structure and `SecondaryStructure` of
// a RNA sequence from its 'dot-bracket' structure.
func SecondaryStructureFromDotBracket(dotBracketStructure string) (annotatedStructure string, secondaryStructure *SecondaryStructure, err error) {
	// sanitize input and check if valid
	dotBracketStructure = strings.TrimSpace(dotBracketStructure)
	isValid, err := checks.IsValidDotBracketStructure(dotBracketStructure)
	if !isValid {
		return "", nil, err
	}

	pairTable, err := pairTableFromDotBracketStructure(dotBracketStructure)
	if err != nil {
		return "", nil, err
	}

	lenStructure := len(dotBracketStructure)
	pc := &parseCompound{
		length:             lenStructure,
		pairTable:          pairTable,
		annotatedStructure: make([]byte, lenStructure),
	}

	annotatedStructure, secondaryStructure = pc.evaluate()
	return
}

/**
* Returns a slice `pairTableFromDotBracketStructure` where `pairTableFromDotBracketStructure[i]` returns the index of the
* nucleotide that that the nucelotide at `i` is paired with, else -1.
* Examples -
* Index:   0  1  2  3  4  5  6  7  8  9 10 11 12 13 14 15 16 17 18 19 20 21 22 23 24 25 26 27 28 29
* Input: " .  .  (  (  (  (  .  .  .  )  )  )  )  .  .  .  (  (  .  .  .  .  .  .  .  .  )  )  .  ."
* Output:[-1 -1 12 11 10  9 -1 -1 -1  5  4  3  2 -1 -1 -1 27 26 -1 -1 -1 -1 -1 -1 -1 -1 17 16 -1 -1]
*
* Index:   0  1  2  3  4  5  6  7  8  9 10 11 12 13 14 15 16 17 18 19 20 21 22 23 24 25 26 27 28 29
* Input: " (  .  (  (  (  (  .  .  .  )  )  )  )  .  .  .  (  (  .  .  .  .  .  .  .  .  )  )  .  )"
* Output:[29 -1 12 11 10  9 -1 -1 -1  5  4  3  2 -1 -1 -1 27 26 -1 -1 -1 -1 -1 -1 -1 -1 17 16 -1  0]
 */
func pairTableFromDotBracketStructure(structure string) ([]int, error) {
	lenStructure := len(structure)
	pairTable := make([]int, lenStructure)

	// the characters of the opening and closing brackets in the structure string
	var openBracket, closeBracket byte = dotBracketStemFivePrimeNucleotide, dotBracketStemThreePrimeNucleotide

	// keeps track of the indexes of the opening base of a base pair.
	// Indexes of the opening base are pushed onto the stack and poped off when
	// the opening base's closing base is encountered
	var openingBaseIdxStack []int = make([]int, 0)

	// iterate through structure and create pair table
	for i := 0; i < lenStructure; i++ {

		if structure[i] == openBracket {
			// we've encountered an opening base of a pair so push onto stack
			openingBaseIdxStack = append(openingBaseIdxStack, i)
		} else if structure[i] == closeBracket {
			// we've encountered a closing pair

			lenStack := len(openingBaseIdxStack)
			if lenStack < 0 {
				return nil,
					fmt.Errorf("%v\nunbalanced brackets '%v%v' found while extracting base pairs",
						structure, openBracket, closeBracket)
			}

			// current index of one-indexed sequence
			lastElemIdx := lenStack - 1

			// the opening base index for the current closing base
			openingBaseIdx := openingBaseIdxStack[lastElemIdx]

			// update the pair table
			pairTable[i] = openingBaseIdx
			pairTable[openingBaseIdx] = i

			// update the stack
			openingBaseIdxStack = openingBaseIdxStack[:lastElemIdx]
		} else {
			// unpaired nucleotide encountered
			pairTable[i] = -1
		}
	}

	if len(openingBaseIdxStack) != 0 {
		return nil, fmt.Errorf("%v\nunbalanced brackets '%v%v' found while extracting base pairs",
			structure, openBracket, closeBracket)
	}

	return pairTable, nil
}

func (pc *parseCompound) evaluate() (annotatedStructure string, secondaryStructure *SecondaryStructure) {
	secondaryStructures := make([]interface{}, 0)
	pairTable := pc.pairTable

	// holds the length of the current exterior loop
	lenExteriorLoop := 0

	for i := 0; i < pc.length; i++ {

		if pairTable[i] == -1 {
			// unpaired nucleotide encountered
			pc.annotatedStructure[i] = exteriorLoopUnpairedNucleotide
			lenExteriorLoop++
			continue
		}

		if lenExteriorLoop != 0 {
			// add single stranded region of exterior loop to structures and reset
			// lenExteriorLoop for next iteration of for-loop
			ssr := &SingleStrandedRegion{
				FivePrimeIdx:  i - lenExteriorLoop,
				ThreePrimeIdx: i - 1,
			}
			secondaryStructures = append(secondaryStructures, ssr)
			lenExteriorLoop = 0
		}

		// we've found an opening base at `i` so evaluate that loop
		structures := evaluateLoop(pc, i)
		secondaryStructures = append(secondaryStructures, structures)

		// seek to end of current loop
		i = pairTable[i]
	}

	// add the single stranded region at the three prime end (if it exists)
	if lenExteriorLoop != 0 {
		ssr := &SingleStrandedRegion{
			FivePrimeIdx:  pc.length - lenExteriorLoop,
			ThreePrimeIdx: pc.length - 1,
		}
		secondaryStructures = append(secondaryStructures, ssr)
	}

	secondaryStructure = &SecondaryStructure{
		Structures: secondaryStructures,
		Length:     pc.length,
	}
	annotatedStructure = string(pc.annotatedStructure)

	return
}

// evaluateLoop evaluates and returns the loop enclosed by (closingFivePrimeIdx,
// closingThreePrimeIdx) which can be either a `Hairpin` or `MultiLoop`.
//
// evaluateLoop proceeds in a stack-wise manner from the closing base pair of
// a stem till the enclosing base pair of the stem. While iterating from the
// closing base pair to the enclosed base pair, each substructure of the
// stem is categorized and added to the list of substructures of the stem.
//
// Once the enclosing base pair of the stem is encountered, `evaluateLoop`
// evaluates whether the loop is a `Hairpin` or `Multiloop` and passes the
// computed stem to the relevant function (`hairpinLoop()` and `multiLoop()`
// for hairpins and multiloops respectively). The relevant function then
// returns the secondary structure (`Hairpin` or `Multiloop`) with the `Stem`
// field of the struct set to the stem computed in this function.
//
// The logic for how a loop is evaluated as a `Hairpin` or `Multiloop` is
// described below:
//
// evaluateLoop proceeds in a stack-wise manner from the 5' and 3' ends of the
// closing base pair until the iterators come across pairs that don't pair
// with each other
// For example,
// 			5' -- X · · U - A ...
// 		   			|    /    |
// 			3' -- Y · V · · B ...
// where (X,Y), (U,V), and (A,B) are base pairs and `·`s are arbitrary number
// of unpaired nucleotides. (X,Y) is the base pair which closes this loop so
// X would be the nucleotide at `closingFivePrimeIdx`.
//
// In this sequence, `evaluateLoop`'s for-loop would proceed normally and
// classify the `StemStructure` closed and enclosed by (X,Y) and (U,V) and
// add the stem structure to the stem.
//
// On the next iteration of the for loop, (U,V) and (A,B) will be classified
// and added to the stem structure of the stem.
//
// But if the continuation of the sequence was
// 			5' -- X · · U - A · ·
// 			 			|    /    |	   ·
// 			3' -- Y · V · · B · ·
// then the last encountered base pair would be (A,B) so `enclosedFivePrimeIdx`
// is A and `enclosedThreePrimeIdx` is B. On the next iteration of the for-loop,
// `enclosedThreePrimeIdx` will point to A and `enclosedThreePrimeIdx` will
// point to B. Thus, the for-loop breaks and since
// `enclosedFivePrimeIdx > enclosedThreePrimeIdx`, we know we have a hairpin.
//
// In the other scenario, suppose the continuation of the sequence was
// 									 			·   ·
// 									 			·   ·
// 								 	 			E - F ...
// 								  			·
// 			5' -- X · · U - A
// 			 			|    /    |
// 			3' -- Y · V · · B · · C - D ...
// 											 			·   ·
// 											 			·   ·
// then on the next iteration, `enclosedFivePrimeIdx` would point to E, and
// `enclosedThreePrimeIdx` would point to C, but since
// `pairTable[enclosedThreePrimeIdx] != enclosedFivePrimeIdx` (as C and E don't
// pair), we know we're in a multiloop.
func evaluateLoop(pc *parseCompound, closingFivePrimeIdx int) interface{} {
	pairTable := pc.pairTable
	closingThreePrimeIdx := pairTable[closingFivePrimeIdx]

	pc.annotatedStructure[closingFivePrimeIdx] = dotBracketStemFivePrimeNucleotide
	pc.annotatedStructure[closingThreePrimeIdx] = dotBracketStemThreePrimeNucleotide

	// init the `Stem` for this loop
	stem := Stem{
		ClosingFivePrimeIdx:  closingFivePrimeIdx,
		ClosingThreePrimeIdx: closingThreePrimeIdx,
	}
	stemStructures := make([]StemStructure, 0)

	// iterator from the 5' to 3' direction
	enclosedFivePrimeIdx := closingFivePrimeIdx

	// iterator from the 3' to 5' direction
	enclosedThreePrimeIdx := closingThreePrimeIdx

	for enclosedFivePrimeIdx < enclosedThreePrimeIdx {
		// process all enclosed `StemStructure`s

		// seek to a base pair from 5' end
		enclosedFivePrimeIdx++
		for pairTable[enclosedFivePrimeIdx] == -1 {
			enclosedFivePrimeIdx++
		}

		// seek to a base pair from 3' end
		enclosedThreePrimeIdx--
		for pairTable[enclosedThreePrimeIdx] == -1 {
			enclosedThreePrimeIdx--
		}

		if pairTable[enclosedThreePrimeIdx] != enclosedFivePrimeIdx || enclosedFivePrimeIdx > enclosedThreePrimeIdx {
			// enclosedFivePrimeIdx & enclosedThreePrimeIdx don't pair. Must have found hairpin or multi-loop.
			break
		} else {
			// We have found a `StemStructure` closed by (`closingFivePrimeIdx`,
			// `closingThreePrimeIdx`) and enclosed by (`enclosedFivePrimeIdx`,
			// `enclosedThreePrimeIdx`)

			// pass the base pairs on to get the StemStructure along with its type.
			stemStructure := stemStructure(pc,
				closingFivePrimeIdx, closingThreePrimeIdx,
				enclosedFivePrimeIdx, enclosedThreePrimeIdx,
			)
			stemStructures = append(stemStructures, stemStructure)

			pc.annotatedStructure[enclosedFivePrimeIdx] = dotBracketStemFivePrimeNucleotide
			pc.annotatedStructure[enclosedThreePrimeIdx] = dotBracketStemThreePrimeNucleotide

			closingFivePrimeIdx = enclosedFivePrimeIdx
			closingThreePrimeIdx = enclosedThreePrimeIdx
		}
	} // end for

	// Set remaining fields of the stem
	if closingFivePrimeIdx == stem.ClosingFivePrimeIdx {
		// stem doesn't have any StemStructures. Thus only consists of its closing
		// base pair.
		stem.EnclosedFivePrimeIdx = -1
		stem.EnclosedThreePrimeIdx = -1

		if len(stemStructures) > 0 {
			// sanity check to ensure there aren't any stem structures
			panic("stem contains StemStructures")
		}
	} else {
		// stem has stem structures
		stem.EnclosedFivePrimeIdx = closingFivePrimeIdx
		stem.EnclosedThreePrimeIdx = closingThreePrimeIdx
		stem.Structures = stemStructures
	}

	if enclosedFivePrimeIdx > enclosedThreePrimeIdx {
		// hairpin
		return hairpin(pc, closingFivePrimeIdx, closingThreePrimeIdx, stem)
	} else {
		// we have a multi-loop
		return multiLoop(pc, closingFivePrimeIdx, stem)
	}
}

// stemStructure sets the required interior loop nucleotides of a
// `parseCompound`'s annotatedStructure and returns the `StemStructure` closed
// by (`closingFivePrimeIdx`, `closingThreePrimeIdx`) and enclosed by
// (`enclosedFivePrimeIdx`, `enclosedThreePrimeIdx`)
func stemStructure(pc *parseCompound,
	closingFivePrimeIdx, closingThreePrimeIdx,
	enclosedFivePrimeIdx, enclosedThreePrimeIdx int) StemStructure {

	for i := closingFivePrimeIdx + 1; i < enclosedFivePrimeIdx; i++ {
		pc.annotatedStructure[i] = interiorLoopUnpairedNucleotide
	}

	for i := enclosedThreePrimeIdx + 1; i < closingThreePrimeIdx; i++ {
		pc.annotatedStructure[i] = interiorLoopUnpairedNucleotide
	}

	stemStructure := NewStemStructure(closingFivePrimeIdx, closingThreePrimeIdx,
		enclosedFivePrimeIdx, enclosedThreePrimeIdx)

	return stemStructure
}

// hairpin sets the required single stranded hairpin nucleotides of a
// `parseCompound`'s annotatedStructure and returns `Hairpin` closed
// by (`closingFivePrimeIdx`, `closingThreePrimeIdx`)
func hairpin(pc *parseCompound, closingFivePrimeIdx, closingThreePrimeIdx int,
	stem Stem) *Hairpin {

	// keep track of whether the hairpin has single stranded nucleotides
	hairpinHasSingleStrandedNucleotides := false
	for i := closingFivePrimeIdx + 1; i < closingThreePrimeIdx; i++ {
		pc.annotatedStructure[i] = hairpinLoopNucleotide
		hairpinHasSingleStrandedNucleotides = true
	}

	var singleStrandedFivePrimeIdx, singleStrandedThreePrimeIdx int
	if hairpinHasSingleStrandedNucleotides {
		singleStrandedFivePrimeIdx = closingFivePrimeIdx + 1
		singleStrandedThreePrimeIdx = closingThreePrimeIdx - 1
	} else {
		// There is no hairpin loop
		singleStrandedFivePrimeIdx = -1
		singleStrandedThreePrimeIdx = -1
	}

	return &Hairpin{
		Stem:                        stem,
		SingleStrandedFivePrimeIdx:  singleStrandedFivePrimeIdx,
		SingleStrandedThreePrimeIdx: singleStrandedThreePrimeIdx,
	}
}

// multiLoop sets the nucleotides present in a multi-loops's single stranded
// region in a `parseCompound`'s annotatedStructure and returns the
// `MultiLoop` closed by (`closingFivePrimeIdx`, `closingThreePrimeIdx`)
func multiLoop(pc *parseCompound, closingFivePrimeIdx int,
	stem Stem) *MultiLoop {
	pairTable := pc.pairTable

	// the substructures present in the multi-loop
	var substructures []interface{}

	if closingFivePrimeIdx >= pairTable[closingFivePrimeIdx] {
		panic("multiLoopEnergy: closingFivePrimeIdx is not 5' base of a pair that closes a loop")
	}

	// the multi-loop's closing three prime base
	var closingThreePrimeIdx int = pairTable[closingFivePrimeIdx]

	// iterator from the 5' to 3' direction
	enclosedFivePrimeIdx := closingFivePrimeIdx + 1

	// keep track of the number of bases in the current single stranded region
	lenMultiLoopSingleStrandedRegion := 0

	// seek to the first stem (i.e. the first enclosed base pair)
	for enclosedFivePrimeIdx <= closingThreePrimeIdx && pairTable[enclosedFivePrimeIdx] == -1 {
		pc.annotatedStructure[enclosedFivePrimeIdx] = multiLoopSingleStrandedRegionNucleotide
		lenMultiLoopSingleStrandedRegion++
		enclosedFivePrimeIdx++
	}

	// add single stranded region of multi-loop to structures and reset
	// lenMultiLoopSingleStrandedRegion for next iteration of for-loop
	if lenMultiLoopSingleStrandedRegion != 0 {
		ssr := &SingleStrandedRegion{
			FivePrimeIdx:  enclosedFivePrimeIdx - lenMultiLoopSingleStrandedRegion,
			ThreePrimeIdx: enclosedFivePrimeIdx - 1,
		}
		substructures = append(substructures, ssr)
		lenMultiLoopSingleStrandedRegion = 0
	}

	for enclosedFivePrimeIdx < closingThreePrimeIdx {
		// process the enclosed substructure
		substructure := evaluateLoop(pc, enclosedFivePrimeIdx)
		substructures = append(substructures, substructure)

		// seek to the next stem
		enclosedFivePrimeIdx = pairTable[enclosedFivePrimeIdx] + 1
		for enclosedFivePrimeIdx < closingThreePrimeIdx && pairTable[enclosedFivePrimeIdx] == -1 {
			pc.annotatedStructure[enclosedFivePrimeIdx] = multiLoopSingleStrandedRegionNucleotide
			lenMultiLoopSingleStrandedRegion++
			enclosedFivePrimeIdx++
		}

		if lenMultiLoopSingleStrandedRegion != 0 {
			ssr := &SingleStrandedRegion{
				FivePrimeIdx:  enclosedFivePrimeIdx - lenMultiLoopSingleStrandedRegion,
				ThreePrimeIdx: enclosedFivePrimeIdx - 1,
			}
			substructures = append(substructures, ssr)
			lenMultiLoopSingleStrandedRegion = 0
		}
	}

	return &MultiLoop{
		Stem:          stem,
		Substructures: substructures,
	}
}
