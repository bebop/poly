package secondary_structure

import (
	"fmt"
	"regexp"
)

/******************************************************************************
July, 01, 2021

This file defines the structs and functions needed to get the 'annotated'
structure and `SecondaryStructure` of a RNA from its 'dot-bracket' notation.

'Dot-bracket' notation of a secondary structure is a string where each
character represents a base. Unpaired nucleotides are represented with a '.'
and base pairs are represented by paranthesis. '(' denotes the opening base
and ')' denotes the closing base of a base pair.
For example, "..((..)).." denotes a hairpin where the bases at index 2 and 7,
and at index 3 and 6 are paired.

'Annotated' structure of a secondary structure is a string where each
character represents a base. In this notation, unpaired nucleotides are set
to a character depending on what part of the secondary structure is in.
The character mapping for the unpaired nucleotides is:
* On the exterior part of a RNA structure (not part of any loop): 'e'
* On the single stranded region of a hairpin: 'h'
* On a single stranded region of a multiloop: 'm'
* Part of an interior loop in a stem of a hairpin or multiloop: 'i'
For example, "..((..)).." would have an annotated structure of "ee((hh))ee".
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
// and `SecondaryStructure` of a RNA sequence
type parseCompound struct {
	length             int   // length of `sequence`
	pairTable          []int // (see `pairTable()`)
	annotatedStructure []byte
}

// FromDotBracket returns the annotated structure and `SecondaryStructure` of
// a RNA sequence from its 'dot-bracket' structure.
func FromDotBracket(dotBracketStructure string) (string, *SecondaryStructure, error) {
	err := ensureValidDotBracketStructure(dotBracketStructure)
	if err != nil {
		return "", nil, err
	}

	pairTable, err := pairTable(dotBracketStructure)
	if err != nil {
		return "", nil, err
	}

	lenStructure := len(dotBracketStructure)
	pc := &parseCompound{
		length:             lenStructure,
		pairTable:          pairTable,
		annotatedStructure: make([]byte, lenStructure),
	}

	secondaryStructure := evaluateParseCompound(pc)

	return string(pc.annotatedStructure), &secondaryStructure, nil
}

/**
* Returns a slice `pairTable` where `pairTable[i]` returns the index of the
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
func pairTable(structure string) ([]int, error) {
	var pairTable []int
	lenStructure := len(structure)

	pairTable = make([]int, lenStructure)

	// the characters of the opening and closing brackets in the structure string
	var openBracket, closeBracket byte = '(', ')'

	// keeps track of the indexes of open brackets. Indexes of open brackets are
	// pushed onto stack and poped off when a closing bracket is encountered
	var openBracketIdxStack []int = make([]int, lenStructure)
	var stackIdx int = 0

	// iterate through structure and create pair table
	for i := 0; i < lenStructure; i++ {
		if structure[i] == openBracket {
			// push index of open bracket onto stack
			openBracketIdxStack[stackIdx] = i
			stackIdx++
		} else if structure[i] == closeBracket {
			stackIdx--

			if stackIdx < 0 {
				return nil,
					fmt.Errorf("%v\nunbalanced brackets '%v%v' found while extracting base pairs",
						structure, openBracket, closeBracket)
			}

			openBracketIdx := openBracketIdxStack[stackIdx]
			// current index of one-indexed sequence
			pairTable[i] = openBracketIdx
			pairTable[openBracketIdx] = i
		} else {
			pairTable[i] = -1
		}
	}

	if stackIdx != 0 {
		return nil, fmt.Errorf("%v\nunbalanced brackets '%v%v' found while extracting base pairs",
			structure, openBracket, closeBracket)
	}

	return pairTable, nil
}

func ensureValidDotBracketStructure(structure string) error {
	dotBracketStructureRegex, _ := regexp.Compile("^[().]+")
	structureIdxs := dotBracketStructureRegex.FindStringIndex(structure)

	if len(structureIdxs) == 0 || structureIdxs[0] != 0 || structureIdxs[1] != len(structure) {
		return fmt.Errorf("found invalid characters in structure. Only dot-bracket notation allowed")
	}
	return nil
}

func evaluateParseCompound(pc *parseCompound) SecondaryStructure {
	secondaryStructures := make([]interface{}, 0)

	pairTable := pc.pairTable
	lenExteriorLoop := 0
	for i := 0; i < pc.length; i++ {
		if pairTable[i] == -1 {
			pc.annotatedStructure[i] = exteriorLoopUnpairedNucleotide
			lenExteriorLoop++
			continue
		}

		if lenExteriorLoop != 0 {
			// add single stranded region of exterior loop to structures and reset
			// lenExteriorLoop for next iteration of for-loop
			ssr := SingleStrandedRegion{
				FivePrimeIdx:  i - lenExteriorLoop,
				ThreePrimeIdx: i - 1,
			}
			secondaryStructures = append(secondaryStructures, ssr)
			lenExteriorLoop = 0
		}

		structures := evaluateLoop(pc, i)
		secondaryStructures = append(secondaryStructures, structures)

		// seek to end of current loop
		i = pairTable[i]
	}

	// add the single stranded region at the three prime end (if it exists)
	if lenExteriorLoop != 0 {
		ssr := SingleStrandedRegion{
			FivePrimeIdx:  pc.length - lenExteriorLoop,
			ThreePrimeIdx: pc.length - 1,
		}
		secondaryStructures = append(secondaryStructures, ssr)
	}

	return SecondaryStructure{
		Structures: secondaryStructures,
		Length:     pc.length,
	}
}

// evaluateLoop evaluates and returns the loop enclosed by (closingFivePrimeIdx,
// closingThreePrimeIdx) which can be either a `Hairpin` or `MultiLoop`.
//
// evaluateLoop proceeds in a stack-wise manner from the closing base pair of
// a stem till the enclosing base pair of the stem. While iterating from the
// closing base pair to the enclosed base pair, each substructure of the
// stem is categoried and added to the list of substructures of the stem.
//
// Once the enclosing base pair of the stem are encountered, `evaluateLoop`
// evaluates whether the loop is a `Hairpin` or `Multiloop` and passes the
// computed stem to the relevant function (`hairpinLoop()` and `multiLoop()`
// respectively). The relevant function then returns the secondary structure
// with the stem set as the stem computed in this function.
//
// The logic for how a loop is evaluated as a `Hairpin` or `Multiloop` is
// described below:

// evaluateLoop proceeds in a stack-wise manner from the 5' and 3' ends of the
// closing base pair until the iterators come across pairs that don't pair
// with each other
// For example,
// 5' -- X · · U - A ...
// 		   |    /    |
// 3' -- Y · V · · B ...
// where (X,Y), (U,V), and (A,B) are base pairs and `·`s are arbitrary number
// of unpaired nucleotides. (X,Y) is the base pair which closes this loop so
// X would be the nucleotide at `closingFivePrimeIdx`.
// In this sequence, `evaluateLoop`'s for-loop would proceed normally and pass
// on (X,Y) and (U,V) to `stackBulgeInteriorLoop` to classify the
// `StemStructure` closed and enclosed by (X,Y) and (U,V).
// On the next iteration of the for loop, (U,V) and (A,B) will be passed on to
// `stackBulgeInteriorLoop` to classify the  `StemStructure` closed and
// enclosed by (U,V) and (A,B).
//
// But if the continuation of the sequence was
// 5' -- X · · U - A · ·
// 			 |    /    |	   ·
// 3' -- Y · V · · B · ·
// then the last encounted base pair would be (A,B) so `enclosedFivePrimeIdx`
// is A and `enclosedThreePrimeIdx` is B. On the next iteration of the for-loop,
// `enclosedThreePrimeIdx` will point to A and `enclosedThreePrimeIdx` will
// point to B. Thus, the for-loop breaks and since
// `enclosedFivePrimeIdx > enclosedThreePrimeIdx`, we know we have a hairpin.
//
// In the other scenario, suppose the contiuation of the sequence was
// 									 ·   ·
// 									 ·   ·
// 								 	 E - F ...
// 								  ·
// 5' -- X · · U - A
// 			 |    /    |
// 3' -- Y · V · · B · · C - D ...
// 											 ·   ·
// 											 ·   ·
// then on the next iteration, `enclosedFivePrimeIdx` would point to E, and
// `enclosedThreePrimeIdx` would point to C, but since
// `pairTable[enclosedThreePrimeIdx] != enclosedFivePrimeIdx` (as C and E don't
// pair), we know we're in some sort of multi-loop.
func evaluateLoop(pc *parseCompound, closingFivePrimeIdx int) interface{} {

	pairTable := pc.pairTable
	closingThreePrimeIdx := pairTable[closingFivePrimeIdx]

	pc.annotatedStructure[closingFivePrimeIdx] = dotBracketStemFivePrimeNucleotide
	pc.annotatedStructure[closingThreePrimeIdx] = dotBracketStemThreePrimeNucleotide

	// init the stem structure for this loop
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
		// base pairs
		if len(stemStructures) > 0 {
			// sanity check
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
// `parseCompound`'s annotatedStructure and returns `StemStructure` closed
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

	stemStructure := newStemStructure(closingFivePrimeIdx, closingThreePrimeIdx,
		enclosedFivePrimeIdx, enclosedThreePrimeIdx)

	return stemStructure
}

// stemStructure sets the required single stranded hairpin nucleotides of a
// `parseCompound`'s annotatedStructure and returns `Hairpin` closed
// by (`closingFivePrimeIdx`, `closingThreePrimeIdx`)
func hairpin(pc *parseCompound, closingFivePrimeIdx, closingThreePrimeIdx int,
	stem Stem) Hairpin {

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
		// There is no hairpin
		singleStrandedFivePrimeIdx = -1
		singleStrandedThreePrimeIdx = -1
	}

	return Hairpin{
		Stem:                        stem,
		SingleStrandedFivePrimeIdx:  singleStrandedFivePrimeIdx,
		SingleStrandedThreePrimeIdx: singleStrandedThreePrimeIdx,
	}
}

// multiLoop sets the nucleotides present in a multi-loops's single stranded
// region in a `parseCompound`'s annotatedStructure and returns the
// `MultiLoop` closed by (`closingFivePrimeIdx`, `closingThreePrimeIdx`)
func multiLoop(pc *parseCompound, closingFivePrimeIdx int,
	stem Stem) MultiLoop {
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
		ssr := SingleStrandedRegion{
			FivePrimeIdx:  enclosedFivePrimeIdx - lenMultiLoopSingleStrandedRegion,
			ThreePrimeIdx: enclosedFivePrimeIdx - 1,
		}
		substructures = append(substructures, ssr)
		lenMultiLoopSingleStrandedRegion = 0
	}

	for enclosedFivePrimeIdx < closingThreePrimeIdx {
		// process all enclosed in the multi-loop
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
			ssr := SingleStrandedRegion{
				FivePrimeIdx:  enclosedFivePrimeIdx - lenMultiLoopSingleStrandedRegion,
				ThreePrimeIdx: enclosedFivePrimeIdx - 1,
			}
			substructures = append(substructures, ssr)
			lenMultiLoopSingleStrandedRegion = 0
		}
	}

	substructuresFivePrimeIdx, substructuresThreePrimeIdx := stem.EnclosedFivePrimeIdx+1, stem.EnclosedThreePrimeIdx-1
	return MultiLoop{
		Stem:                       stem,
		SubstructuresFivePrimeIdx:  stem.EnclosedFivePrimeIdx + 1,
		SubstructuresThreePrimeIdx: stem.EnclosedThreePrimeIdx - 1,
		Substructures: SecondaryStructure{
			Structures: substructures,
			Length:     substructuresThreePrimeIdx - substructuresFivePrimeIdx + 1,
		},
	}
}

/************************************************************************

The following are functions to help test the code the `FromDotBracket` function.

They are for internal use only, but also serve as an example for how to parse
a `SecondaryStructure` recursively.

***********************************************************************/

// doDotBracketFromSecondaryStructure returns the dot-bracket structure of
// a `SecondaryStructure`.
func dotBracketFromSecondaryStructure(secondaryStructure SecondaryStructure) string {
	return doDotBracketFromSecondaryStructure(secondaryStructure, 0)
}

// Since we'd like to call this code recursively for the substructures
// enclosed in a multi-loop, we include a `offset` param as the index fields of
// the substructures of a `MultiLoop` have an absolute reference to the indexes
// of the original `SecondaryStructure`, but when we process a `MultiLoop`
// recursively, we need to get the relative reference for the indexes for the
// output.
func doDotBracketFromSecondaryStructure(secondaryStructure SecondaryStructure, offset int) string {
	var dotBracket []byte = make([]byte, secondaryStructure.Length)

	for _, structure := range secondaryStructure.Structures {
		switch structure := structure.(type) {
		case SingleStrandedRegion:
			for i := structure.FivePrimeIdx; i <= structure.ThreePrimeIdx; i++ {
				dotBracket[i-offset] = dotBracketUnpairedNucleotide
			}
		case MultiLoop:
			dotBracketFromStem(&dotBracket, structure.Stem, offset)
			substructuresDotBracket := doDotBracketFromSecondaryStructure(structure.Substructures, structure.SubstructuresFivePrimeIdx)
			lenSubstructuresDotBracket := len(substructuresDotBracket)
			if lenSubstructuresDotBracket != structure.SubstructuresThreePrimeIdx-structure.SubstructuresFivePrimeIdx+1 {
				panic("len of dot bracket from substructures != len substructure")
			}
			for i, j := structure.SubstructuresFivePrimeIdx, 0; i <= structure.SubstructuresThreePrimeIdx; i++ {
				dotBracket[i] = substructuresDotBracket[j]
				j++
			}
		case Hairpin:
			dotBracketFromStem(&dotBracket, structure.Stem, offset)
			if structure.SingleStrandedFivePrimeIdx != -1 {
				for i := structure.SingleStrandedFivePrimeIdx; i <= structure.SingleStrandedThreePrimeIdx; i++ {
					dotBracket[i-offset] = dotBracketUnpairedNucleotide
				}
			}
		}
	}

	return string(dotBracket)
}

func dotBracketFromStem(dotBracket *[]byte, stem Stem, offset int) {
	(*dotBracket)[stem.ClosingFivePrimeIdx-offset] = dotBracketStemFivePrimeNucleotide
	(*dotBracket)[stem.ClosingThreePrimeIdx-offset] = dotBracketStemThreePrimeNucleotide
	for _, stemStructure := range stem.Structures {
		dotBracketFromStemStructure(dotBracket, stemStructure, offset)
	}
}

func dotBracketFromStemStructure(dotBracket *[]byte, stemStructure StemStructure, offset int) {
	for i := stemStructure.ClosingFivePrimeIdx + 1; i < stemStructure.EnclosedFivePrimeIdx; i++ {
		(*dotBracket)[i-offset] = dotBracketUnpairedNucleotide
	}
	(*dotBracket)[stemStructure.EnclosedFivePrimeIdx-offset] = dotBracketStemFivePrimeNucleotide

	for i := stemStructure.EnclosedThreePrimeIdx + 1; i < stemStructure.ClosingThreePrimeIdx; i++ {
		(*dotBracket)[i-offset] = dotBracketUnpairedNucleotide
	}
	(*dotBracket)[stemStructure.EnclosedThreePrimeIdx-offset] = dotBracketStemThreePrimeNucleotide
}
