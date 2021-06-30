package seconday_structure

import (
	"fmt"
	"regexp"
)

/******************************************************************************

 ******************************************************************************/

// parseCompound holds all information needed to compute the free energy of a
// RNA secondary structure (a RNA sequence along with its folded structure).
type parseCompound struct {
	length             int   // length of `sequence`
	pairTable          []int // (see `pairTable()`)
	annotatedStructure []byte
}

func FromDotBracket(dotBracketStructure string) (string, SecondaryStructure) {
	err := ensureValidAnnotatedStructure(dotBracketStructure)
	if err != nil {
		panic(err)
	}

	pairTable, err := pairTable(dotBracketStructure)
	if err != nil {
		panic(err)
	}

	lenStructure := len(dotBracketStructure)
	fc := &parseCompound{
		length:             lenStructure,
		pairTable:          pairTable,
		annotatedStructure: make([]byte, lenStructure),
	}

	secondaryStructure := evaluateParseCompound(fc)

	return string(fc.annotatedStructure), secondaryStructure
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

func ensureValidAnnotatedStructure(structure string) error {
	annotatedStructureRegex, _ := regexp.Compile("^[().]+")
	structureIdxs := annotatedStructureRegex.FindStringIndex(structure)

	if structureIdxs[0] != 0 || structureIdxs[1] != len(structure) {
		return fmt.Errorf("found invalid characters in structure. Only annotated structure notation allowed")
	}
	return nil
}

var (
	exteriorLoopUnpairedNucleotide byte = 'e'
	interiorLoopUnpairedNucleotide byte = 'i'
	hairpinLoopNucleotide          byte = 'h'
	multiLoopNucleotide            byte = 'm'
	fivePrimeStemNucleotide        byte = '('
	threePrimeStemNucleotide       byte = ')'
)

func evaluateParseCompound(fc *parseCompound) SecondaryStructure {
	secondaryStructures := make([]interface{}, 0)

	// get structure of exterior loop
	exteriorLoopStructure(fc)

	annotatedStructure := fc.annotatedStructure
	lenExteriorLoop := 0
	for i := 0; i < fc.length; i++ {
		if annotatedStructure[i] == exteriorLoopUnpairedNucleotide {
			lenExteriorLoop++
			continue
		}

		if lenExteriorLoop != 0 {
			ssr := SingleStrandedRegion{
				FivePrimeIdx:  i - lenExteriorLoop,
				ThreePrimeIdx: i - 1,
			}
			secondaryStructures = append(secondaryStructures, ssr)
			lenExteriorLoop = 0
		}

		structures := stackEnergy(fc, i)
		secondaryStructures = append(secondaryStructures, structures)
		// seek to end of current loop
		i = fc.pairTable[i]
	}

	if lenExteriorLoop != 0 {
		ssr := SingleStrandedRegion{
			FivePrimeIdx:  fc.length - lenExteriorLoop,
			ThreePrimeIdx: fc.length - 1,
		}
		secondaryStructures = append(secondaryStructures, ssr)
	}

	return SecondaryStructure{
		Structures: secondaryStructures,
		Length:     fc.length,
	}
}

/**
* Calculate the energy contribution of stabilizing dangling-ends/mismatches
* for all stems branching off the exterior loop.
* For example, if the structure is
* Structure: " .  .  (  (  (  (  .  .  .  )  )  )  )  .  .  .  (  (  .  .  .  .  .  .  .  .  )  )  .  ."
* Index:       0  1  2  3  4  5  6  7  8  9 10 11 12 13 14 15 16 17 18 19 20 21 22 23 24 25 26 27 28 29
* then the exterior loops in the structure are as follows: (2, 12) and (16, 27).
* See `exteriorStemEnergy()` for information of how exterior loops are evaluated.
* More information available at: https://rna.urmc.rochester.edu/NNDB/turner04/exterior.html
 */
func exteriorLoopStructure(fc *parseCompound) {
	pairTable := fc.pairTable
	length := fc.length
	pairFivePrimeIdx := 0

	// seek to opening base of first stem
	for pairFivePrimeIdx < length && pairTable[pairFivePrimeIdx] == -1 {
		fc.annotatedStructure[pairFivePrimeIdx] = exteriorLoopUnpairedNucleotide
		pairFivePrimeIdx++
	}

	for pairFivePrimeIdx < length {
		/* pairFivePrimeIdx must have a pairing partner */
		pairThreePrimeIdx := pairTable[pairFivePrimeIdx]

		// seek to the next stem
		pairFivePrimeIdx = pairThreePrimeIdx + 1
		for pairFivePrimeIdx < length && pairTable[pairFivePrimeIdx] == -1 {
			fc.annotatedStructure[pairFivePrimeIdx] = exteriorLoopUnpairedNucleotide
			pairFivePrimeIdx++
		}
	}
}

/**
    stackEnergy recursively calculate energy of substructure enclosed by
    (closingFivePrimeIdx, closingThreePrimeIdx).

    vivek: I think this function should be better named.
    Given a base pair (closingFivePrimeIdx, closingThreePrimeIdx), this function
    finds all loops that the base pair encloses. For each enclosed base pair,
    we find out if it belongs to one of the three groups:
    - stacking pair, bulge, or interior loop
    - hairpin loop
    - multi-loop
    and pass on the indexes of the closing pair and enclosed pair to the relevant
    functions.
    I think it's called stackEnergy because it proceeds in a stack-wise manner
    from the 5' and 3' ends until the iterators come across pairs that don't pair
    with each other
		For example,
    5' -- X · · U - A ...
					|    /    |
		3' -- Y · V · · B ...
		where (X,Y), (U,V), and (A,B) are base pairs and `·`s are arbitrary number
		of unpaired nucleotides. (X,Y) is the base pair which closes this loop so
		X would be the nucleotide at `closingFivePrimeIdx`.

		In this sequence, `stackEnergy`'s for-loop would proceed normally and pass on
		(U,V) and (A,B) to `stackBulgeInteriorLoopEnergy`.

		But if the continuation of the sequence was
		5' -- X · · U - A · ·
					|    /    |	   ·
		3' -- Y · V · · B · ·
		then the last encounted base pair would be (A,B) so `enclosedFivePrimeIdx`
		is A and `enclosedThreePrimeIdx` is B. On the next iteration of the for-loop,
		`enclosedThreePrimeIdx` will point to A and `enclosedThreePrimeIdx` will point
		to B. Thus, the for-loop breaks and since `enclosedFivePrimeIdx > enclosedThreePrimeIdx`,
		we know we have a hairpin loop.

		In the other scenario, suppose the contiuation of the sequence was
											·   ·
											·   ·
											E - F ...
										 ·
		5' -- X · · U - A
					|    /    |
		3' -- Y · V · · B · · C - D ...
													·   ·
													·   ·
		then on the next iteration, `enclosedFivePrimeIdx` would point to E, and
		`enclosedThreePrimeIdx` would point to C, but since
		`pairTable[enclosedThreePrimeIdx] != enclosedFivePrimeIdx` (as C and E don't
		pair), we know we're in some sort of multi-loop. Note that multi-loops
		have atleast two enclosed pairs.
*/
func stackEnergy(fc *parseCompound, closingFivePrimeIdx int) interface{} {
	// energyContributions := make([]StructuralMotifWithEnergy, 0)

	pairTable := fc.pairTable
	closingThreePrimeIdx := pairTable[closingFivePrimeIdx]

	fc.annotatedStructure[closingFivePrimeIdx] = fivePrimeStemNucleotide
	fc.annotatedStructure[closingThreePrimeIdx] = threePrimeStemNucleotide

	stem := Stem{
		ClosingFivePrimeIdx:  closingFivePrimeIdx,
		ClosingThreePrimeIdx: closingThreePrimeIdx,
	}

	stemStructures := make([]StemStructure, 0)

	// iterator from the 5' to 3' direction starting at `pairFivePrimeIdx`
	enclosedFivePrimeIdx := closingFivePrimeIdx

	// iterator from the 3' to 5' direction starting at `pairThreePrimeIdx`
	enclosedThreePrimeIdx := closingThreePrimeIdx

	for enclosedFivePrimeIdx < enclosedThreePrimeIdx {
		// process all enclosed stacks and interior loops

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
			// we have either a stacking pair, bulge, or interior loop
			stemStructure := stackBulgeInteriorLoopEnergy(fc,
				closingFivePrimeIdx, closingThreePrimeIdx,
				enclosedFivePrimeIdx, enclosedThreePrimeIdx,
			)
			stemStructures = append(stemStructures, stemStructure)

			fc.annotatedStructure[enclosedFivePrimeIdx] = fivePrimeStemNucleotide
			fc.annotatedStructure[enclosedThreePrimeIdx] = threePrimeStemNucleotide

			closingFivePrimeIdx = enclosedFivePrimeIdx
			closingThreePrimeIdx = enclosedThreePrimeIdx
		}
	} /* end for */

	// Finished with the stem
	stem.EnclosedFivePrimeIdx = closingFivePrimeIdx
	stem.EnclosedThreePrimeIdx = closingThreePrimeIdx
	stem.structures = stemStructures

	if enclosedFivePrimeIdx > enclosedThreePrimeIdx {
		// hairpin
		return hairpinLoopEnergy(fc, closingFivePrimeIdx, closingThreePrimeIdx, stem)
	} else {
		// we have a multi-loop
		return multiLoopEnergy(fc, closingFivePrimeIdx, stem)
	}
}

/**
 *  Evaluate the free energy contribution of an interior loop with delimiting
 *  base pairs (i,j) and (k,l). See `evaluateStackBulgeInteriorLoop()` for more details.
 */
func stackBulgeInteriorLoopEnergy(fc *parseCompound,
	closingFivePrimeIdx, closingThreePrimeIdx,
	enclosedFivePrimeIdx, enclosedThreePrimeIdx int) StemStructure {

	for i := closingFivePrimeIdx + 1; i < enclosedFivePrimeIdx; i++ {
		fc.annotatedStructure[i] = interiorLoopUnpairedNucleotide
	}

	for i := enclosedThreePrimeIdx + 1; i < closingThreePrimeIdx; i++ {
		fc.annotatedStructure[i] = interiorLoopUnpairedNucleotide
	}

	stemStructure := StemStructure{
		ClosingFivePrimeIdx: closingFivePrimeIdx, EnclosedFivePrimeIdx: enclosedFivePrimeIdx,
		EnclosedThreePrimeIdx: enclosedThreePrimeIdx, ClosingThreePrimeIdx: closingThreePrimeIdx,
	}

	stemStructure.setStructureType()

	return stemStructure
}

func hairpinLoopEnergy(fc *parseCompound, pairFivePrimeIdx, pairThreePrimeIdx int, stem Stem) HairpinLoop {

	hairpinHasSingleStrandedNucleotides := false
	for i := pairFivePrimeIdx + 1; i < pairThreePrimeIdx; i++ {
		hairpinHasSingleStrandedNucleotides = true
		fc.annotatedStructure[i] = hairpinLoopNucleotide
	}

	var singleStrandedFivePrimeIdx, singleStrandedThreePrimeIdx int
	if hairpinHasSingleStrandedNucleotides {
		singleStrandedFivePrimeIdx = pairFivePrimeIdx + 1
		singleStrandedThreePrimeIdx = pairThreePrimeIdx - 1
	} else {
		// There is no hairpin
		singleStrandedFivePrimeIdx = -1
		singleStrandedThreePrimeIdx = -1
	}

	return HairpinLoop{
		Stem:                        stem,
		StemFivePrimeIdx:            stem.ClosingFivePrimeIdx,
		StemThreePrimeIdx:           stem.ClosingThreePrimeIdx,
		SingleStrandedFivePrimeIdx:  singleStrandedFivePrimeIdx,
		SingleStrandedThreePrimeIdx: singleStrandedThreePrimeIdx,
	}
}

func multiLoopEnergy(fc *parseCompound, closingFivePrimeIdx int, stem Stem) MultiLoop {
	// energyContributions := make([]StructuralMotifWithEnergy, 0)
	pairTable := fc.pairTable

	var substructures []interface{}

	// (pairFivePrimeIdx, pairThreePrimeIdx) is exterior pair of multi loop
	if closingFivePrimeIdx >= pairTable[closingFivePrimeIdx] {
		panic("multiLoopEnergy: pairFivePrimeIdx is not 5' base of a pair that closes a loop")
	}

	var closingThreePrimeIdx int
	if closingFivePrimeIdx == 0 {
		closingThreePrimeIdx = fc.length + 1
	} else {
		closingThreePrimeIdx = pairTable[closingFivePrimeIdx]
	}

	// The index of the enclosed base pair at the 5' end
	enclosedFivePrimeIdx := closingFivePrimeIdx + 1

	lenMultiLoopSingleStrandedRegion := 0

	// seek to the first stem (i.e. the first enclosed base pair)
	for enclosedFivePrimeIdx <= closingThreePrimeIdx && pairTable[enclosedFivePrimeIdx] == -1 {
		fc.annotatedStructure[enclosedFivePrimeIdx] = multiLoopNucleotide
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

	for enclosedFivePrimeIdx < closingThreePrimeIdx {
		// add up the contributions of the substructures of the multi loop
		substructure := stackEnergy(fc, enclosedFivePrimeIdx)
		substructures = append(substructures, substructure)

		enclosedThreePrimeIdx := pairTable[enclosedFivePrimeIdx]

		// seek to the next stem
		enclosedFivePrimeIdx = enclosedThreePrimeIdx + 1

		for enclosedFivePrimeIdx < closingThreePrimeIdx && pairTable[enclosedFivePrimeIdx] == -1 {
			fc.annotatedStructure[enclosedFivePrimeIdx] = multiLoopNucleotide
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

		// add unpaired nucleotides
		// nbUnpairedNucleotides += enclosedFivePrimeIdx - enclosedThreePrimeIdx - 1
	}

	substructuresFivePrimeIdx, substructuresThreePrimeIdx := stem.EnclosedFivePrimeIdx+1, stem.EnclosedThreePrimeIdx-1
	return MultiLoop{
		Stem:                       stem,
		StemFivePrimeIdx:           stem.ClosingFivePrimeIdx,
		StemThreePrimeIdx:          stem.ClosingThreePrimeIdx,
		SubstructuresFivePrimeIdx:  stem.EnclosedFivePrimeIdx + 1,
		SubstructuresThreePrimeIdx: stem.EnclosedThreePrimeIdx - 1,
		Substructures: SecondaryStructure{
			Structures: substructures,
			Length:     substructuresThreePrimeIdx - substructuresFivePrimeIdx + 1,
		},
	}
}

const (
	unpairedNucleotide byte = '.'
)

func DotBracket(secondaryStructure SecondaryStructure, offset int) string {
	var dotBracket []byte = make([]byte, secondaryStructure.Length)
	for _, structure := range secondaryStructure.Structures {
		switch structure := structure.(type) {
		case SingleStrandedRegion:
			for i := structure.FivePrimeIdx; i <= structure.ThreePrimeIdx; i++ {
				dotBracket[i-offset] = unpairedNucleotide
			}
		case MultiLoop:
			dotBracketFromStem(&dotBracket, structure.Stem, offset)
			substructuresDotBracket := DotBracket(structure.Substructures, structure.SubstructuresFivePrimeIdx)
			lenSubstructuresDotBracket := len(substructuresDotBracket)
			if lenSubstructuresDotBracket != structure.SubstructuresThreePrimeIdx-structure.SubstructuresFivePrimeIdx+1 {
				panic("len of dot bracket from substructures != len substructure")
			}
			for i, j := structure.SubstructuresFivePrimeIdx, 0; i <= structure.SubstructuresThreePrimeIdx; i++ {
				dotBracket[i] = substructuresDotBracket[j]
				j++
			}
		case HairpinLoop:
			dotBracketFromStem(&dotBracket, structure.Stem, offset)
			if structure.SingleStrandedFivePrimeIdx != -1 {
				for i := structure.SingleStrandedFivePrimeIdx; i <= structure.SingleStrandedThreePrimeIdx; i++ {
					dotBracket[i-offset] = unpairedNucleotide
				}
			}
		}
	}

	return string(dotBracket)
}

func dotBracketFromStem(dotBracket *[]byte, stem Stem, offset int) {
	(*dotBracket)[stem.ClosingFivePrimeIdx-offset] = fivePrimeStemNucleotide
	(*dotBracket)[stem.ClosingThreePrimeIdx-offset] = threePrimeStemNucleotide
	for _, stemStructure := range stem.structures {
		dotBracketFromStemStructure(dotBracket, stemStructure, offset)
	}
}

func dotBracketFromStemStructure(dotBracket *[]byte, stemStructure StemStructure, offset int) {
	for i := stemStructure.ClosingFivePrimeIdx + 1; i < stemStructure.EnclosedFivePrimeIdx; i++ {
		(*dotBracket)[i-offset] = unpairedNucleotide
	}

	(*dotBracket)[stemStructure.EnclosedFivePrimeIdx-offset] = fivePrimeStemNucleotide

	for i := stemStructure.EnclosedThreePrimeIdx + 1; i < stemStructure.ClosingThreePrimeIdx; i++ {
		(*dotBracket)[i-offset] = unpairedNucleotide
	}

	(*dotBracket)[stemStructure.EnclosedThreePrimeIdx-offset] = threePrimeStemNucleotide
}
