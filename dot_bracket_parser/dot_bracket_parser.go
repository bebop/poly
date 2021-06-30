package dot_bracket_parser

import (
	"fmt"
	"regexp"
)

/******************************************************************************

 ******************************************************************************/

// A RNA's SecondaryStructure is composed of a list of MultiLoops, HairpinLoops,
// and SingleStrandedRegions. Note that since Go doesn't support inheritance,
// we use `interface{}` as the type for the list, but the only types that are
// allowed/used are `MultiLoop`, `HairpinLoop` and `SingleStrandedRegion`
type SecondaryStructure struct {
	structures []interface{}
	length     int
}

type MultiLoop struct {
	stemFivePrimeIdx, stemThreePrimeIdx                   int
	stem                                                  Stem
	subStructuresFivePrimeIdx, subStructuresThreePrimeIdx int
	subStructures                                         SecondaryStructure
}

type HairpinLoop struct {
	stemFivePrimeIdx, stemThreePrimeIdx                     int
	stem                                                    Stem
	singleStrandedFivePrimeIdx, singleStrandedThreePrimeIdx int
}

type SingleStrandedRegion struct {
	fivePrimeIdx, threePrimeIdx int
}

type Stem struct {
	closingFivePrimeIdx, enclosedFivePrimeIdx   int
	enclosedThreePrimeIdx, closingThreePrimeIdx int
	structures                                  []StemStructure
}

type StemStructureType int

const (
	StackingPair StemStructureType = iota
	Bulge
	Interior1x1Loop
	Interior2x1Loop
	Interior1xnLoop
	Interior2x2Loop
	Interior2x3Loop
	GenericInteriorLoop
)

type StemStructure struct {
	closingFivePrimeIdx, enclosedFivePrimeIdx   int
	enclosedThreePrimeIdx, closingThreePrimeIdx int
	structureType                               StemStructureType
}

func (structure *StemStructure) setStructureType() {
	nbUnpairedFivePrime := structure.enclosedFivePrimeIdx - structure.closingFivePrimeIdx - 1
	nbUnpairedThreePrime := structure.closingThreePrimeIdx - structure.enclosedThreePrimeIdx - 1

	var nbUnpairedLarger, nbUnpairedSmaller int

	if nbUnpairedFivePrime > nbUnpairedThreePrime {
		nbUnpairedLarger = nbUnpairedFivePrime
		nbUnpairedSmaller = nbUnpairedThreePrime
	} else {
		nbUnpairedLarger = nbUnpairedThreePrime
		nbUnpairedSmaller = nbUnpairedFivePrime
	}

	if nbUnpairedLarger == 0 {
		// stacking pair
		structure.structureType = StackingPair
		return
	}

	if nbUnpairedSmaller == 0 {
		// bulge
		structure.structureType = Bulge
		return
	} else {
		// interior loop
		if nbUnpairedSmaller == 1 {
			if nbUnpairedLarger == 1 {
				// 1x1 loop
				structure.structureType = Interior1x1Loop
				return
			}

			if nbUnpairedLarger == 2 {
				// 2x1 loop
				structure.structureType = Interior2x1Loop
				return
			} else {
				// 1xn loop
				structure.structureType = Interior1xnLoop
				return
			}
		} else if nbUnpairedSmaller == 2 {
			if nbUnpairedLarger == 2 {
				// 2x2 loop
				structure.structureType = Interior2x2Loop
				return
			} else if nbUnpairedLarger == 3 {
				// 2x3 loop
				structure.structureType = Interior2x3Loop
				return
			}
		}

		{
			/* generic interior loop (no else here!)*/
			structure.structureType = GenericInteriorLoop
			return
		}
	}
}

// parseCompound holds all information needed to compute the free energy of a
// RNA secondary structure (a RNA sequence along with its folded structure).
type parseCompound struct {
	length             int   // length of `sequence`
	pairTable          []int // (see `pairTable()`)
	annotatedStructure []byte
}

func SecondaryStructureFromDotBracket(dotBracketStructure string) (string, SecondaryStructure) {
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
				fivePrimeIdx:  i - lenExteriorLoop,
				threePrimeIdx: i - 1,
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
			fivePrimeIdx:  fc.length - lenExteriorLoop,
			threePrimeIdx: fc.length - 1,
		}
		secondaryStructures = append(secondaryStructures, ssr)
	}

	return SecondaryStructure{
		structures: secondaryStructures,
		length:     fc.length,
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
		closingFivePrimeIdx:  closingFivePrimeIdx,
		closingThreePrimeIdx: closingThreePrimeIdx,
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
	stem.enclosedFivePrimeIdx = closingFivePrimeIdx
	stem.enclosedThreePrimeIdx = closingThreePrimeIdx
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
		closingFivePrimeIdx: closingFivePrimeIdx, enclosedFivePrimeIdx: enclosedFivePrimeIdx,
		enclosedThreePrimeIdx: enclosedThreePrimeIdx, closingThreePrimeIdx: closingThreePrimeIdx,
	}

	stemStructure.setStructureType()

	return stemStructure

	// something happens here if the ends of the pair are on different strands
	// energy := evaluateStackBulgeInteriorLoop(nbUnpairedFivePrime, nbUnpairedThreePrime,
	// 	closingBasePairType, enclosedBasePairType,
	// 	fc.encodedSequence[closingFivePrimeIdx+1], fc.encodedSequence[closingThreePrimeIdx-1],
	// 	fc.encodedSequence[enclosedFivePrimeIdx-1], fc.encodedSequence[enclosedThreePrimeIdx+1], fc.energyParams)

	// energyContribution := StructuralMotifWithEnergy{
	// 	structure: StructuralMotif{
	// 		closingFivePrimeIdx:   closingFivePrimeIdx,
	// 		closingThreePrimeIdx:  closingThreePrimeIdx,
	// 		enclosedFivePrimeIdx:  enclosedFivePrimeIdx,
	// 		enclosedThreePrimeIdx: enclosedThreePrimeIdx,
	// 		loopType:              InteriorLoop,
	// 	},
	// }
	// return energy, energyContribution
}

/**
 *  Compute the energy of either a stacking pair, bulge, or interior loop.
 *  This function computes the free energy of a loop with the following
 *  structure:
 *        3'  5'
 *        |   |
 *        U - V
 *    a_n       b_1
 *     .        .
 *     .        .
 *     .        .
 *    a_1       b_m
 *        X - Y
 *        |   |
 *        5'  3'
 *  This general structure depicts an interior loop that is closed by the base pair (X,Y).
 *  The enclosed base pair is (V,U) which leaves the unpaired bases a_1-a_n and b_1-b_n
 *  that constitute the loop.
 *  The base pair (X,Y) will be refered to as `closingBasePair`, and the base pair
 *  (U,V) will be refered to as `enclosedBasePair`.
 *  In this example, the length of the interior loop is `n+m`
 *  where n or m may be 0 resulting in a bulge-loop or base pair stack.
 *  The mismatching nucleotides for the closing pair (X,Y) are:
 *  5'-mismatch: a_1
 *  3'-mismatch: b_m
 *  and for the enclosed base pair (V,U):
 *  5'-mismatch: b_1
 *  3'-mismatch: a_n
 *  @note Base pairs are always denoted in 5'->3' direction. Thus the `enclosedBasePair`
 *  must be 'turned around' when evaluating the free energy of the interior loop.
 *
 *  More information about
 * 		- Bulge Loops: https://rna.urmc.rochester.edu/NNDB/turner04/bulge.html
 *		- Interior Loops: https://rna.urmc.rochester.edu/NNDB/turner04/internal.html
 *
 *  @param  nbUnpairedLeftLoop      The size of the 'left'-loop (number of unpaired nucleotides)
 *  @param  nbUnpairedRightLoop     The size of the 'right'-loop (number of unpaired nucleotides)
 *  @param  closingBasePairType    	The encoded type of the closing base pair of the interior loop
 *  @param  enclosedBasePairType    The encoded type of the enclosed base pair
 *  @param  closingFivePrimeMismatch    The 5'-mismatching encoded nucleotide of the closing pair
 *  @param  closingThreePrimeMismatch   The 3'-mismatching encoded nucleotide of the closing pair
 *  @param  enclosedThreePrimeMismatch   The 3'-mismatching encoded nucleotide of the enclosed pair
 *  @param  enclosedFivePrimeMismatch    The 5'-mismatching encoded nucleotide of the enclosed pair
 *  @param  P       The datastructure containing scaled energy parameters
 *  @return The Free energy of the Interior-loop in dcal/mol
 */
// func evaluateStackBulgeInteriorLoop(nbUnpairedLeftLoop, nbUnpairedRightLoop int,
// 	closingBasePairType, enclosedBasePairType int,
// 	closingFivePrimeMismatch, closingThreePrimeMismatch,
// 	enclosedThreePrimeMismatch, enclosedFivePrimeMismatch int,
// 	energyParams *energyParams) int {
// 	/* compute energy of degree 2 loop (stack bulge or interior) */
// 	var nbUnpairedLarger, nbUnpairedSmaller int

// 	if nbUnpairedLeftLoop > nbUnpairedRightLoop {
// 		nbUnpairedLarger = nbUnpairedLeftLoop
// 		nbUnpairedSmaller = nbUnpairedRightLoop
// 	} else {
// 		nbUnpairedLarger = nbUnpairedRightLoop
// 		nbUnpairedSmaller = nbUnpairedLeftLoop
// 	}

// 	if nbUnpairedLarger == 0 {
// 		// stacking pair
// 		return energyParams.stackingPair[closingBasePairType][enclosedBasePairType]
// 	}

// 	if nbUnpairedSmaller == 0 {
// 		// bulge
// 		var energy int

// 		if nbUnpairedLarger <= maxLenLoop {
// 			energy = energyParams.bulge[nbUnpairedLarger]
// 		} else {
// 			energy = energyParams.bulge[maxLenLoop] + int(energyParams.lxc*math.Log(float64(nbUnpairedLarger)/float64(maxLenLoop)))
// 		}

// 		if nbUnpairedLarger == 1 {
// 			energy += energyParams.stackingPair[closingBasePairType][enclosedBasePairType]
// 		} else {
// 			if closingBasePairType > 2 {
// 				// The closing base pair is not a GC or CG pair (could be GU, UG, AU, UA, or non-standard)
// 				energy += energyParams.terminalAUPenalty
// 			}

// 			if enclosedBasePairType > 2 {
// 				// The encosed base pair is not a GC or CG pair (could be GU, UG, AU, UA, or non-standard)
// 				energy += energyParams.terminalAUPenalty
// 			}
// 		}

// 		return energy
// 	} else {
// 		// interior loop
// 		if nbUnpairedSmaller == 1 {
// 			if nbUnpairedLarger == 1 {
// 				// 1x1 loop
// 				return energyParams.interior1x1Loop[closingBasePairType][enclosedBasePairType][closingFivePrimeMismatch][closingThreePrimeMismatch]
// 			}

// 			if nbUnpairedLarger == 2 {
// 				// 2x1 loop
// 				// vivek: shouldn't it be a 1x2 loop (following the convention that the first integer corresponds to `nbUnpairedSmaller`)? No consistency in naming here.
// 				if nbUnpairedLeftLoop == 1 {
// 					return energyParams.interior2x1Loop[closingBasePairType][enclosedBasePairType][closingFivePrimeMismatch][enclosedFivePrimeMismatch][closingThreePrimeMismatch]
// 				} else {
// 					return energyParams.interior2x1Loop[enclosedBasePairType][closingBasePairType][enclosedFivePrimeMismatch][closingFivePrimeMismatch][enclosedThreePrimeMismatch]
// 				}
// 			} else {
// 				// 1xn loop

// 				var energy int
// 				// vivek: Why do we add 1 here?
// 				if nbUnpairedLarger+1 <= maxLenLoop {
// 					energy = energyParams.interiorLoop[nbUnpairedLarger+1]
// 				} else {
// 					energy = energyParams.interiorLoop[maxLenLoop] + int(energyParams.lxc*math.Log((float64(nbUnpairedLarger)+1.0)/float64(maxLenLoop)))
// 				}
// 				energy += min(energyParams.maxNinio, (nbUnpairedLarger-nbUnpairedSmaller)*energyParams.ninio)
// 				energy += energyParams.mismatch1xnInteriorLoop[closingBasePairType][closingFivePrimeMismatch][closingThreePrimeMismatch] + energyParams.mismatch1xnInteriorLoop[enclosedBasePairType][enclosedFivePrimeMismatch][enclosedThreePrimeMismatch]
// 				return energy
// 			}
// 		} else if nbUnpairedSmaller == 2 {
// 			if nbUnpairedLarger == 2 {
// 				// 2x2 loop
// 				return energyParams.interior2x2Loop[closingBasePairType][enclosedBasePairType][closingFivePrimeMismatch][enclosedThreePrimeMismatch][enclosedFivePrimeMismatch][closingThreePrimeMismatch]
// 			} else if nbUnpairedLarger == 3 {
// 				// 2x3 loop

// 				var energy int
// 				energy = energyParams.interiorLoop[nbNucleobase+1] + energyParams.ninio
// 				energy += energyParams.mismatch2x3InteriorLoop[closingBasePairType][closingFivePrimeMismatch][closingThreePrimeMismatch] + energyParams.mismatch2x3InteriorLoop[enclosedBasePairType][enclosedFivePrimeMismatch][enclosedThreePrimeMismatch]
// 				return energy
// 			}
// 		}

// 		{
// 			/* generic interior loop (no else here!)*/
// 			nbUnpairedNucleotides := nbUnpairedLarger + nbUnpairedSmaller

// 			var energy int
// 			if nbUnpairedNucleotides <= maxLenLoop {
// 				energy = energyParams.interiorLoop[nbUnpairedNucleotides]
// 			} else {
// 				energy = energyParams.interiorLoop[maxLenLoop] + int(energyParams.lxc*math.Log(float64(nbUnpairedNucleotides)/float64(maxLenLoop)))
// 			}

// 			energy += min(energyParams.maxNinio, (nbUnpairedLarger-nbUnpairedSmaller)*energyParams.ninio)

// 			energy += energyParams.mismatchInteriorLoop[closingBasePairType][closingFivePrimeMismatch][closingThreePrimeMismatch] + energyParams.mismatchInteriorLoop[enclosedBasePairType][enclosedFivePrimeMismatch][enclosedThreePrimeMismatch]
// 			return energy
// 		}
// 	}
// }

/**
 *  Evaluate free energy of a hairpin loop encosed by the base pair
 *  (pairFivePrimeIdx, pairThreePrimeIdx).
 *  See `evaluateHairpinLoop()` for more information.
 *
 *  @param  fc                 The foldCompund for the particular energy evaluation
 *  @param  pairFivePrimeIdx   5'-position of the base pair enclosing the hairpin loop
 *  @param  pairThreePrimeIdx  3'-position of the base pair enclosing the hairpin loop
 *  @returns                   Free energy of the hairpin loop closed by (pairFivePrimeIdx,pairThreePrimeIdx) in dcal/mol
 */
func hairpinLoopEnergy(fc *parseCompound, pairFivePrimeIdx, pairThreePrimeIdx int, stem Stem) HairpinLoop {
	// nbUnpairedNucleotides := pairThreePrimeIdx - pairFivePrimeIdx - 1 // also the size of the hairpin loop

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
		stem:                        stem,
		stemFivePrimeIdx:            stem.closingFivePrimeIdx,
		stemThreePrimeIdx:           stem.closingThreePrimeIdx,
		singleStrandedFivePrimeIdx:  singleStrandedFivePrimeIdx,
		singleStrandedThreePrimeIdx: singleStrandedThreePrimeIdx,
	}

	// energyContribution := StructuralMotifWithEnergy{
	// 	structure: StructuralMotif{
	// 		closingFivePrimeIdx:  pairFivePrimeIdx,
	// 		closingThreePrimeIdx: pairThreePrimeIdx,
	// 		loopType:             HairpinLoop,
	// 	},
	// 	energy: energy,
	// }
	// return energy, energyContribution
}

/**
 *  Compute the energy of a hairpin loop.
 *
 *  To evaluate the free energy of a hairpin loop, several parameters have to be known.
 *  A general hairpin-loop has this structure:
 *        a3 a4
 *      a2     a5
 *      a1     a6
 *        X - Y
 *        |   |
 *        5'  3'
 *  where X-Y marks the closing pair [e.g. a (G,C) pair]. The length of this loop is 6 as there are
 *  six unpaired nucleotides (a1-a6) enclosed by (X,Y). The 5' mismatching nucleotide is
 *  a1 while the 3' mismatch is a6. The nucleotide sequence of this loop is [a1 a2 a3 a4 a5 a6]
 *  @note The parameter `sequence` should contain the sequence of the loop in capital letters of the nucleic acid
 *  alphabet if the loop size is below 7. This is useful for unusually stable tri-, tetra- and hexa-loops
 *  which are treated differently (based on experimental data) if they are tabulated.
 *
 *  More information available at: https://rna.urmc.rochester.edu/NNDB/turner04/hairpin.html
 *
 *  @param  size                 The size of the hairpin loop (number of unpaired nucleotides)
 *  @param  basePairType         The pair type of the base pair closing the hairpin
 *  @param  fivePrimeMismatch    The 5'-mismatching encoded nucleotide
 *  @param  threePrimeMismatch   The 3'-mismatching encoded nucleotide
 *  @param  sequence						 The sequence of the loop
 *  @param  energyParams     		 The datastructure containing scaled energy parameters
 *  @return The free energy of the hairpin loop in dcal/mol
 */
// func evaluateHairpinLoop(size, basePairType, fivePrimeMismatch, threePrimeMismatch int, sequence string, energyParams *energyParams) int {
// 	var energy int

// 	if size <= maxLenLoop {
// 		energy = energyParams.hairpinLoop[size]
// 	} else {
// 		energy = energyParams.hairpinLoop[maxLenLoop] + int(energyParams.lxc*math.Log(float64(size)/float64(maxLenLoop)))
// 	}

// 	if size < 3 {
// 		// should only be the case when folding alignments
// 		return energy
// 	}

// 	if size == 4 {
// 		tetraLoop := sequence[:6]
// 		tetraLoopEnergy, present := energyParams.tetraLoop[tetraLoop]
// 		if present {
// 			return tetraLoopEnergy
// 		}
// 	} else if size == 6 {
// 		hexaLoop := sequence[:8]
// 		hexaLoopEnergy, present := energyParams.hexaLoop[hexaLoop]
// 		if present {
// 			return hexaLoopEnergy
// 		}
// 	} else if size == 3 {
// 		triLoop := sequence[:5]
// 		triLoopEnergy, present := energyParams.triLoop[triLoop]
// 		if present {
// 			return triLoopEnergy
// 		}

// 		if basePairType > 2 {
// 			// (X,Y) is not a GC or CG pair (could be GU, UG, AU, UA, or non-standard)
// 			return energy + energyParams.terminalAUPenalty
// 		}

// 		// (X,Y) is a GC or CG pair
// 		return energy
// 	}

// 	energy += energyParams.mismatchHairpinLoop[basePairType][fivePrimeMismatch][threePrimeMismatch]

// 	return energy
// }

/**
 *  Compute the energy of a multi-loop.
 *  A multi-loop has this structure:
 *				 	 		·   ·
 * 				 	 		·   ·
 * 				 	 		·   ·
 * 				 	 		A - B
 * 				 	 	 ·		 	 ·
 * 				 	  ·				C · · ·
 * 				 	 ·					|
 * 		5' -- X					D · · ·
 * 				  |				 ·
 * 		3' -- Y · V - U
 * 				 	 	  ·   ·
 * 				 	 	  ·   ·
 * 				 	 	  ·   ·
 * where X-Y marks the closing pair (X is the nucleotide at `closingFivePrimeIdx`)
 * of the multi-loop, and `·`s are an arbitrary number of unparied nucelotides
 * (can also be 0). (A,B), (C,D), and (V,U) are the base pairs enclosed in this
 * multi-loop. Note that there is at minimum two base pairs enclosed in a
 * multi-loop.
 *
 * multiLoopEnergy iterates through the multi-loop and finds each base pair
 * enclosed in the multi-loop. For each enclosed base pair, we add to the total
 * the energy due to that base pair and due to the substructure closed by the
 * enclosed base pair.
 * We also add a bonus energy based on the number of unpaired nucleotides in
 * the multi-loop (though it defaults to 0 right now).
 *
 * More information about multi-loops: https://rna.urmc.rochester.edu/NNDB/turner04/mb.html
 */
func multiLoopEnergy(fc *parseCompound, closingFivePrimeIdx int, stem Stem) MultiLoop {
	// energyContributions := make([]StructuralMotifWithEnergy, 0)
	pairTable := fc.pairTable

	var subStructures []interface{}

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
			fivePrimeIdx:  enclosedFivePrimeIdx - lenMultiLoopSingleStrandedRegion,
			threePrimeIdx: enclosedFivePrimeIdx - 1,
		}
		subStructures = append(subStructures, ssr)
		lenMultiLoopSingleStrandedRegion = 0
	}

	for enclosedFivePrimeIdx < closingThreePrimeIdx {
		// add up the contributions of the substructures of the multi loop
		subStructure := stackEnergy(fc, enclosedFivePrimeIdx)
		subStructures = append(subStructures, subStructure)

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
				fivePrimeIdx:  enclosedFivePrimeIdx - lenMultiLoopSingleStrandedRegion,
				threePrimeIdx: enclosedFivePrimeIdx - 1,
			}
			subStructures = append(subStructures, ssr)
			lenMultiLoopSingleStrandedRegion = 0
		}

		// add unpaired nucleotides
		// nbUnpairedNucleotides += enclosedFivePrimeIdx - enclosedThreePrimeIdx - 1
	}

	subStructuresFivePrimeIdx, subStructuresThreePrimeIdx := stem.enclosedFivePrimeIdx+1, stem.enclosedThreePrimeIdx-1
	return MultiLoop{
		stem:                       stem,
		stemFivePrimeIdx:           stem.closingFivePrimeIdx,
		stemThreePrimeIdx:          stem.closingThreePrimeIdx,
		subStructuresFivePrimeIdx:  stem.enclosedFivePrimeIdx + 1,
		subStructuresThreePrimeIdx: stem.enclosedThreePrimeIdx - 1,
		subStructures: SecondaryStructure{
			structures: subStructures,
			length:     subStructuresThreePrimeIdx - subStructuresFivePrimeIdx + 1,
		},
	}
}

const (
	unpairedNucleotide byte = '.'
)

func DotBracket(secondaryStructure SecondaryStructure, offset int) string {
	var dotBracket []byte = make([]byte, secondaryStructure.length)
	for _, structure := range secondaryStructure.structures {
		switch structure.(type) {
		case SingleStrandedRegion:
			ssr := structure.(SingleStrandedRegion)
			for i := ssr.fivePrimeIdx; i <= ssr.threePrimeIdx; i++ {
				dotBracket[i-offset] = unpairedNucleotide
			}
		case MultiLoop:
			multiLoop := structure.(MultiLoop)
			dotBracketFromStem(&dotBracket, multiLoop.stem, offset)
			subStructuresDotBracket := DotBracket(multiLoop.subStructures, multiLoop.subStructuresFivePrimeIdx)
			lenSubStructuresDotBracket := len(subStructuresDotBracket)
			if lenSubStructuresDotBracket != multiLoop.subStructuresThreePrimeIdx-multiLoop.subStructuresFivePrimeIdx+1 {
				panic("len of dot bracket from substructures != len substructure")
			}
			for i, j := multiLoop.subStructuresFivePrimeIdx, 0; i <= multiLoop.subStructuresThreePrimeIdx; i++ {
				dotBracket[i] = subStructuresDotBracket[j]
				j++
			}
		case HairpinLoop:
			hairpinLoop := structure.(HairpinLoop)
			dotBracketFromStem(&dotBracket, hairpinLoop.stem, offset)
			if hairpinLoop.singleStrandedFivePrimeIdx != -1 {
				for i := hairpinLoop.singleStrandedFivePrimeIdx; i <= hairpinLoop.singleStrandedThreePrimeIdx; i++ {
					dotBracket[i-offset] = unpairedNucleotide
				}
			}
		}
	}

	return string(dotBracket)
}

func dotBracketFromStem(dotBracket *[]byte, stem Stem, offset int) {
	(*dotBracket)[stem.closingFivePrimeIdx-offset] = fivePrimeStemNucleotide
	(*dotBracket)[stem.closingThreePrimeIdx-offset] = threePrimeStemNucleotide
	for _, stemStructure := range stem.structures {
		dotBracketFromStemStructure(dotBracket, stemStructure, offset)
	}
}

func dotBracketFromStemStructure(dotBracket *[]byte, stemStructure StemStructure, offset int) {
	for i := stemStructure.closingFivePrimeIdx + 1; i < stemStructure.enclosedFivePrimeIdx; i++ {
		(*dotBracket)[i-offset] = unpairedNucleotide
	}

	(*dotBracket)[stemStructure.enclosedFivePrimeIdx-offset] = fivePrimeStemNucleotide

	for i := stemStructure.enclosedThreePrimeIdx + 1; i < stemStructure.closingThreePrimeIdx; i++ {
		(*dotBracket)[i-offset] = unpairedNucleotide
	}

	(*dotBracket)[stemStructure.enclosedThreePrimeIdx-offset] = threePrimeStemNucleotide
}
