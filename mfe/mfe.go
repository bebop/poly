package mfe

import (
	"errors"
	"fmt"
	"math"
	"regexp"
	"strings"

	"github.com/TimothyStiles/poly/energy_params"
	. "github.com/TimothyStiles/poly/secondary_structure"
)

/******************************************************************************
May, 26, 2021

The MFE calculation functionality (including comments and description of
	functions) has been taken directly from ViennaRNA with the follwing changes:
- ViennaRNA includes the ability to specify a dangle model (more info available
	at [src/bin/RNAeval.ggo#L153](https://github.com/ViennaRNA/ViennaRNA/blob/d6fbaf2b33d65944f9249a91ed5ab4b3277f7d06/src/bin/RNAeval.ggo#L153)).
	This implementation keeps it simple and defaults to the value of -d2.
- ViennaRNA includes the ability to specify hard and soft constraints (more info
	available from [their ppt explaining hard and soft constrains](http://benasque.org/2015rna/talks_contr/2713_lorenz_benasque_2015.pdf),
	[official docs](https://www.tbi.univie.ac.at/RNA/ViennaRNA/doc/html/group__constraints.html),
	and [their paper](https://almob.biomedcentral.com/articles/10.1186/s13015-016-0070-z)).
	This implementation keeps it simple and defaults to no hard or soft
	constraints.
- ViennaRNA includes the ability to calculate the minimum free energy of
	co-folded sequences. This implementation keeps it simple and defaults to
	calculating the mfe of only single sequences.

File is structured as so:

	Structs:
		foldCompound - holds all information needed to compute the free energy of a
			RNA secondary structure (a RNA sequence along with its folded structure).
		EnergyContribution - holds information about the energy contribution of a
			single loop, and the indexes of the base pairs that delimit the loop in
			the RNA secondary structure.

	Big functions from this file:

		MinimumFreeEnergy - given a RNA sequence, it's structure, and a temperature,
			it returns the minimum free energy of the RNA secondary structure at the
			specified temperature. A slice of `EnergyContribution` is also returned
			which allows for in-depth examination of the energy contribution of each
			loop present in the secondary structure.

This file contains everything you need to caluclate minimum free energy of a
RNA sequence (except for the energy parameters which are located in
`energy_params.go`).

TODO: Biological context

 ******************************************************************************/

const (
	// DefaultTemperature is the temperature in Celcius for free energy evaluation
	DefaultTemperature float64 = 37.0
)

type DanglingEndsModel int

const (
	// NoDanglingEnds specifies no dangling end energies to be added to energy calcualtions
	NoDanglingEnds DanglingEndsModel = 0
	// DoubleDanglingEnds specifies energies due to dangling ends (on both five and three prime sides)
	// should be added to energy calcualtions
	DoubleDanglingEnds DanglingEndsModel = 2
	// DefaultDanglingEndsModel defaults to DoubleDangles
	DefaultDanglingEndsModel = DoubleDanglingEnds
)

// foldCompound holds all information needed to compute the free energy of a
// RNA secondary structure (a RNA sequence along with its folded structure).
type foldCompound struct {
	length              int                         // length of `sequence`
	energyParams        *energy_params.EnergyParams // The precomputed free energy contributions for each type of loop
	sequence            string                      // The input sequence string
	encodedSequence     []int                       // Numerical encoding of the sequence (see `encodeSequence()` for more information)
	pairTable           []int                       // (see `pairTable()`)
	basePairEncodedType map[byte]map[byte]int       // (see `basePairEncodedTypeMap()`)
	dangleModel         DanglingEndsModel
}

//  MinimumFreeEnergy returns the free energy of an already folded RNA and a list of
//  energy contribution of each loop.
//
//  @param sequence         		A RNA sequence
//  @param structure        		Secondary structure in dot-bracket notation
//  @param temperature      		Temperature at which to evaluate the free energy
//															 of the structure
//  @param danglingEndsModel  	Specify whether to include energy from dangling
// 															ends (NoDanglingEnds or DoubleDangles)
//  @return                 		The free energy of the input structure given the
// 															input sequence in kcal/mol, and a slice of
// 															`StructuralMotifWithEnergy` (which gives
// 															information about the energy contribution of
// 															each loop present in the secondary structure)
func MinimumFreeEnergy(sequence, structure string, temperature float64, energyParamsSet energy_params.EnergyParamsSet, danglingEndsModel DanglingEndsModel) (float64, *SecondaryStructure, error) {
	lenSequence := len(sequence)
	lenStructure := len(structure)

	if lenSequence != lenStructure {
		return 0, nil, fmt.Errorf("length of sequence (%v) != length of structure (%v)", lenSequence, lenStructure)
	} else if lenStructure == 0 {
		return 0, nil, errors.New("lengths of sequence and structure cannot be 0")
	}

	sequence = strings.ToUpper(sequence)

	err := ensureValidRNA(sequence)
	if err != nil {
		return 0, nil, err
	}

	err = ensureValidStructure(structure)
	if err != nil {
		return 0, nil, err
	}

	pairTable, err := pairTable(structure)
	if err != nil {
		return 0, nil, err
	}

	fc := &foldCompound{
		length:              lenSequence,
		energyParams:        energy_params.NewEnergyParams(energyParamsSet, temperature),
		sequence:            sequence,
		encodedSequence:     encodeSequence(sequence),
		pairTable:           pairTable,
		basePairEncodedType: basePairEncodedTypeMap(),
		dangleModel:         danglingEndsModel,
	}

	secondaryStructure := evaluateFoldCompound(fc)

	return float64(secondaryStructure.Energy) / 100.0, &secondaryStructure, nil
}

func ensureValidRNA(sequence string) error {
	rnaRegex, _ := regexp.Compile("^[ACGU]+")
	rnaIdxs := rnaRegex.FindStringIndex(sequence)

	if rnaIdxs[0] != 0 || rnaIdxs[1] != len(sequence) {
		return fmt.Errorf("found invalid characters in RNA sequence. Only A, C, G, and U allowed")
	}
	return nil
}

func ensureValidStructure(structure string) error {
	dotBracketStructureRegex, _ := regexp.Compile("^[().]+")
	structureIdxs := dotBracketStructureRegex.FindStringIndex(structure)

	if structureIdxs[0] != 0 || structureIdxs[1] != len(structure) {
		return fmt.Errorf("found invalid characters in structure. Only dot-bracket notation allowed")
	}
	return nil
}

var nucleotideRuneEncodedIntMap map[byte]int = map[byte]int{
	'A': 1,
	'C': 2,
	'G': 3,
	'U': 4,
}

// Encodes a sequence into its numerical representaiton
func encodeSequence(sequence string) []int {
	lenSequence := len(sequence)
	encodedSequence := make([]int, lenSequence)

	// encode the sequence based on nucleotideRuneMap
	for i := 0; i < lenSequence; i++ {
		encodedSequence[i] = nucleotideRuneEncodedIntMap[sequence[i]]
	}

	return encodedSequence
}

/**
Returns a map that encodes a base pair to its numerical representation.
Various loop energy parameters depend in general on the pairs closing the loops.
Internally, the library distinguishes 8 types of pairs (CG=1, GC=2, GU=3, UG=4,
AU=5, UA=6, nonstandard=7, 0= no pair).
The map is:
   _  A  C  G  U
_ {0, 0, 0, 0, 0}
A {0, 0, 0, 0, 5}
C {0, 0, 0, 1, 0}
G {0, 0, 2, 0, 3}
U {0, 6, 0, 4, 0}
For example, map['A']['U'] is 5.
The encoded numerical representation of a base pair carries no meaning in
itself except for where to find the relevent energy contributions of the base
pair in the matrices of the energy paramaters (found in `energy_params.go`).
Thus, any change to this map must be reflected in the energy
parameters matrices in `energy_params.go`, and vice versa.
*/
func basePairEncodedTypeMap() map[byte]map[byte]int {
	nucelotideAEncodedTypeMap := map[byte]int{'U': 5}
	nucelotideCEncodedTypeMap := map[byte]int{'G': 1}
	nucelotideGEncodedTypeMap := map[byte]int{'C': 2, 'U': 3}
	nucelotideUEncodedTypeMap := map[byte]int{'A': 6, 'G': 4}

	basePairMap := map[byte]map[byte]int{
		'A': nucelotideAEncodedTypeMap,
		'C': nucelotideCEncodedTypeMap,
		'G': nucelotideGEncodedTypeMap,
		'U': nucelotideUEncodedTypeMap,
	}

	return basePairMap
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

func evaluateFoldCompound(fc *foldCompound) SecondaryStructure {
	// get the energy contributions due to exterior loops
	var energy int = exteriorLoopsEnergy(fc)
	secondaryStructure := SecondaryStructure{
		ExteriorLoopsEnergy: energy,
		Length:              fc.length,
	}

	secondaryStructures := make([]interface{}, 0)
	pairTable := fc.pairTable
	lenExteriorLoop := 0
	for i := 0; i < fc.length; i++ {
		if pairTable[i] == -1 {
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

		en, structure := evaluateLoop(fc, i)
		energy += en
		secondaryStructures = append(secondaryStructures, structure)

		// seek to end of current loop
		i = pairTable[i]
	}

	// add the single stranded region at the three prime end (if it exists)
	if lenExteriorLoop != 0 {
		ssr := SingleStrandedRegion{
			FivePrimeIdx:  fc.length - lenExteriorLoop,
			ThreePrimeIdx: fc.length - 1,
		}
		secondaryStructures = append(secondaryStructures, ssr)
	}

	secondaryStructure.Energy = energy
	secondaryStructure.Structures = secondaryStructures

	return secondaryStructure
}

/**
* (see `basePairEncodedTypeMap()`)
 */
func encodedBasePairType(fc *foldCompound, basePairFivePrimeIdx, basePairThreePrimeIdx int) int {
	return fc.basePairEncodedType[fc.sequence[basePairFivePrimeIdx]][fc.sequence[basePairThreePrimeIdx]]
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
func exteriorLoopsEnergy(fc *foldCompound) int {
	pairTable := fc.pairTable
	var energy int = 0
	length := fc.length
	pairFivePrimeIdx := 0

	// seek to opening base of first stem
	for pairFivePrimeIdx < length && pairTable[pairFivePrimeIdx] == -1 {
		pairFivePrimeIdx++
	}

	for pairFivePrimeIdx < length {
		/* pairFivePrimeIdx must have a pairing partner */
		pairThreePrimeIdx := pairTable[pairFivePrimeIdx]

		/* get encoded type of base pair (pairFivePrimeIdx, pairThreePrimeIdx) */
		basePairType := encodedBasePairType(fc, pairFivePrimeIdx,
			pairThreePrimeIdx)

		switch fc.dangleModel {
		case NoDanglingEnds:
			energy += EvaluateExteriorStem(basePairType, -1, -1, fc.energyParams)

		case DoubleDanglingEnds:
			var fivePrimeMismatch, threePrimeMismatch int
			if pairFivePrimeIdx > 0 {
				fivePrimeMismatch = fc.encodedSequence[pairFivePrimeIdx-1]
			} else {
				fivePrimeMismatch = -1
			}

			if pairThreePrimeIdx < length-1 {
				threePrimeMismatch = fc.encodedSequence[pairThreePrimeIdx+1]
			} else {
				threePrimeMismatch = -1
			}

			energy += EvaluateExteriorStem(basePairType, fivePrimeMismatch, threePrimeMismatch, fc.energyParams)

		}

		// seek to the next stem
		pairFivePrimeIdx = pairThreePrimeIdx + 1
		for pairFivePrimeIdx < length && pairTable[pairFivePrimeIdx] == -1 {
			pairFivePrimeIdx++
		}
	}

	return energy
}

// Evaluate a stem branching off the exterior loop.
//
// Given a base pair (i,j) (encoded by `basePairType`), compute the
// energy contribution including dangling-end/terminal-mismatch contributions.
//
// You can prevent taking 5'-, 3'-dangles or mismatch contributions into
// account by passing -1 for `fivePrimeMismatch` and/or `threePrimeMismatch`
// respectively.
//
// @param  basePairType   The encoded base pair type of (i, j) (see `basePairType()`)
// @param  fivePrimeMismatch     The encoded nucleotide directly adjacent (in the 5' direction) to i (may be -1 if index of i is 0)
// @param  threePrimeMismatch    The encoded nucleotide directly adjacent (in the 3' direction) to j (may be -1 if index of j is len(sequence) - 1)
// @param  energyParams          The pre-computed energy parameters
// @return                       The energy contribution of the introduced exterior-loop stem
//
func EvaluateExteriorStem(basePairType int, fivePrimeMismatch int, threePrimeMismatch int, energyParams *energy_params.EnergyParams) int {
	var energy int = 0

	if fivePrimeMismatch >= 0 && threePrimeMismatch >= 0 {
		energy += energyParams.MismatchExteriorLoop[basePairType][fivePrimeMismatch][threePrimeMismatch]
	} else if fivePrimeMismatch >= 0 {
		// `j` is the last nucleotide of the sequence
		energy += energyParams.Dangle5[basePairType][fivePrimeMismatch]
	} else if threePrimeMismatch >= 0 {
		// `i` is the first nucleotide of the sequence
		energy += energyParams.Dangle3[basePairType][threePrimeMismatch]
	}

	if basePairType > 2 {
		// The closing base pair is not a GC or CG pair (could be GU, UG, AU, UA, or non-standard)
		energy += energyParams.TerminalAUPenalty
	}

	return energy
}

/**
    evaluateLoop recursively calculate energy of substructure enclosed by
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
    I think it's called evaluateLoop because it proceeds in a stack-wise manner
    from the 5' and 3' ends until the iterators come across pairs that don't pair
    with each other
		For example,
    5' -- X · · U - A ...
					|    /    |
		3' -- Y · V · · B ...
		where (X,Y), (U,V), and (A,B) are base pairs and `·`s are arbitrary number
		of unpaired nucleotides. (X,Y) is the base pair which closes this loop so
		X would be the nucleotide at `closingFivePrimeIdx`.

		In this sequence, `evaluateLoop`'s for-loop would proceed normally and pass on
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
func evaluateLoop(fc *foldCompound, closingFivePrimeIdx int) (int, interface{}) {

	pairTable := fc.pairTable
	var stemEnergy int = 0
	closingThreePrimeIdx := pairTable[closingFivePrimeIdx]

	if encodedBasePairType(fc, closingFivePrimeIdx, closingThreePrimeIdx) == 0 {
		panic(fmt.Sprintf("bases %v and %v (%v%v) can't pair!",
			closingFivePrimeIdx, closingThreePrimeIdx,
			string(fc.sequence[closingFivePrimeIdx]),
			string(fc.sequence[closingThreePrimeIdx])))
	}

	// init the stem for this loop
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
			// We have found a `StemStructure` closed by (`closingFivePrimeIdx`,
			// `closingThreePrimeIdx`) and enclosed by (`enclosedFivePrimeIdx`,
			// `enclosedThreePrimeIdx`)

			if encodedBasePairType(fc, enclosedThreePrimeIdx, enclosedFivePrimeIdx) == 0 {
				panic(fmt.Sprintf("bases %d and %d (%c%c) can't pair!",
					enclosedFivePrimeIdx, enclosedThreePrimeIdx,
					fc.sequence[enclosedThreePrimeIdx],
					fc.sequence[enclosedFivePrimeIdx]))
			}

			// pass the base pairs on to get the StemStructure along with its type.
			stemStructure := stemStructure(fc,
				closingFivePrimeIdx, closingThreePrimeIdx,
				enclosedFivePrimeIdx, enclosedThreePrimeIdx,
			)
			stemStructures = append(stemStructures, stemStructure)
			stemEnergy += stemStructure.Energy

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
		stem.Energy = stemEnergy
	}

	if enclosedFivePrimeIdx > enclosedThreePrimeIdx {
		// hairpin
		hairpin := hairpin(fc, closingFivePrimeIdx, closingThreePrimeIdx, stem)
		return stemEnergy + hairpin.Energy, hairpin
	} else {
		// we have a multi-loop
		multiLoop := multiLoop(fc, closingFivePrimeIdx, stem)
		// multiloop's total energy consists of the energy of the substructures, the loop and the stem
		return stemEnergy + multiLoop.Energy + multiLoop.SubstructuresEnergy, multiLoop
	}
}

/**
 *  Evaluate the free energy contribution of an interior loop with delimiting
 *  base pairs (i,j) and (k,l). See `evaluateStackBulgeInteriorLoop()` for more details.
 */
func stemStructure(fc *foldCompound,
	closingFivePrimeIdx, closingThreePrimeIdx,
	enclosedFivePrimeIdx, enclosedThreePrimeIdx int) StemStructure {

	stemStructure := NewStemStructure(closingFivePrimeIdx, closingThreePrimeIdx,
		enclosedFivePrimeIdx, enclosedThreePrimeIdx)

	stemStructure.Energy = stemStructureEnergy(stemStructure, fc)

	return stemStructure
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
func stemStructureEnergy(stemStructure StemStructure, fc *foldCompound) int {

	closingBasePairType := encodedBasePairType(fc, stemStructure.ClosingFivePrimeIdx, stemStructure.ClosingThreePrimeIdx)
	enclosedBasePairType := encodedBasePairType(fc, stemStructure.EnclosedThreePrimeIdx, stemStructure.EnclosedFivePrimeIdx)

	closingFivePrimeMismatch := fc.encodedSequence[stemStructure.ClosingFivePrimeIdx+1]
	closingThreePrimeMismatch := fc.encodedSequence[stemStructure.ClosingThreePrimeIdx-1]
	enclosedThreePrimeMismatch := fc.encodedSequence[stemStructure.EnclosedFivePrimeIdx-1]
	enclosedFivePrimeMismatch := fc.encodedSequence[stemStructure.EnclosedThreePrimeIdx+1]

	return EvaluateStemStructure(stemStructure,
		closingBasePairType, enclosedBasePairType,
		closingFivePrimeMismatch, closingThreePrimeMismatch,
		enclosedThreePrimeMismatch, enclosedFivePrimeMismatch,
		fc.energyParams)

}

func EvaluateStemStructure(stemStructure StemStructure,
	closingBasePairType, enclosedBasePairType,
	closingFivePrimeMismatch, closingThreePrimeMismatch,
	enclosedThreePrimeMismatch, enclosedFivePrimeMismatch int,
	energyParams *energy_params.EnergyParams) int {
	nbUnpairedFivePrime := stemStructure.EnclosedFivePrimeIdx - stemStructure.ClosingFivePrimeIdx - 1
	nbUnpairedThreePrime := stemStructure.ClosingThreePrimeIdx - stemStructure.EnclosedThreePrimeIdx - 1
	var nbUnpairedLarger, nbUnpairedSmaller int

	if nbUnpairedFivePrime > nbUnpairedThreePrime {
		nbUnpairedLarger = nbUnpairedFivePrime
		nbUnpairedSmaller = nbUnpairedThreePrime
	} else {
		nbUnpairedLarger = nbUnpairedThreePrime
		nbUnpairedSmaller = nbUnpairedFivePrime
	}

	var energy int
	switch stemStructure.Type {
	case StackingPair:
		energy = energyParams.StackingPair[closingBasePairType][enclosedBasePairType]
	case Bulge:
		if nbUnpairedLarger <= energy_params.MaxLenLoop {
			energy = energyParams.Bulge[nbUnpairedLarger]
		} else {
			energy = energyParams.Bulge[energy_params.MaxLenLoop] + int(energyParams.LXC*math.Log(float64(nbUnpairedLarger)/float64(energy_params.MaxLenLoop)))
		}

		if nbUnpairedLarger == 1 {
			energy += energyParams.StackingPair[closingBasePairType][enclosedBasePairType]
		} else {
			if closingBasePairType > 2 {
				// The closing base pair is not a GC or CG pair (could be GU, UG, AU, UA, or non-standard)
				energy += energyParams.TerminalAUPenalty
			}

			if enclosedBasePairType > 2 {
				// The encosed base pair is not a GC or CG pair (could be GU, UG, AU, UA, or non-standard)
				energy += energyParams.TerminalAUPenalty
			}
		}
	case Interior1x1Loop:
		energy = energyParams.Interior1x1Loop[closingBasePairType][enclosedBasePairType][closingFivePrimeMismatch][closingThreePrimeMismatch]
	case Interior2x1Loop:
		if nbUnpairedFivePrime == 1 {
			energy = energyParams.Interior2x1Loop[closingBasePairType][enclosedBasePairType][closingFivePrimeMismatch][enclosedFivePrimeMismatch][closingThreePrimeMismatch]
		} else {
			energy = energyParams.Interior2x1Loop[enclosedBasePairType][closingBasePairType][enclosedFivePrimeMismatch][closingFivePrimeMismatch][enclosedThreePrimeMismatch]
		}
	case Interior1xnLoop:
		if nbUnpairedLarger+1 <= energy_params.MaxLenLoop {
			energy = energyParams.InteriorLoop[nbUnpairedLarger+1]
		} else {
			energy = energyParams.InteriorLoop[energy_params.MaxLenLoop] + int(energyParams.LXC*math.Log((float64(nbUnpairedLarger)+1.0)/float64(energy_params.MaxLenLoop)))
		}
		energy += min(energyParams.MaxNinio, (nbUnpairedLarger-nbUnpairedSmaller)*energyParams.Ninio)
		energy += energyParams.Mismatch1xnInteriorLoop[closingBasePairType][closingFivePrimeMismatch][closingThreePrimeMismatch] + energyParams.Mismatch1xnInteriorLoop[enclosedBasePairType][enclosedFivePrimeMismatch][enclosedThreePrimeMismatch]
	case Interior2x2Loop:
		energy = energyParams.Interior2x2Loop[closingBasePairType][enclosedBasePairType][closingFivePrimeMismatch][enclosedThreePrimeMismatch][enclosedFivePrimeMismatch][closingThreePrimeMismatch]
	case Interior2x3Loop:
		energy = energyParams.InteriorLoop[energy_params.NBBases+1] + energyParams.Ninio
		energy += energyParams.Mismatch2x3InteriorLoop[closingBasePairType][closingFivePrimeMismatch][closingThreePrimeMismatch] + energyParams.Mismatch2x3InteriorLoop[enclosedBasePairType][enclosedFivePrimeMismatch][enclosedThreePrimeMismatch]
	case GenericInteriorLoop:
		nbUnpairedNucleotides := nbUnpairedLarger + nbUnpairedSmaller

		if nbUnpairedNucleotides <= energy_params.MaxLenLoop {
			energy = energyParams.InteriorLoop[nbUnpairedNucleotides]
		} else {
			energy = energyParams.InteriorLoop[energy_params.MaxLenLoop] + int(energyParams.LXC*math.Log(float64(nbUnpairedNucleotides)/float64(energy_params.MaxLenLoop)))
		}

		energy += min(energyParams.MaxNinio, (nbUnpairedLarger-nbUnpairedSmaller)*energyParams.Ninio)

		energy += energyParams.MismatchInteriorLoop[closingBasePairType][closingFivePrimeMismatch][closingThreePrimeMismatch] + energyParams.MismatchInteriorLoop[enclosedBasePairType][enclosedFivePrimeMismatch][enclosedThreePrimeMismatch]
	}

	return energy
}

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
func hairpin(fc *foldCompound, closingFivePrimeIdx, closingThreePrimeIdx int, stem Stem) Hairpin {
	nbUnpairedNucleotides := closingThreePrimeIdx - closingFivePrimeIdx - 1 // also the size of the hairpin loop
	basePairType := encodedBasePairType(fc, closingFivePrimeIdx, closingThreePrimeIdx)

	energy := EvaluateHairpinLoop(nbUnpairedNucleotides, basePairType,
		fc.encodedSequence[closingFivePrimeIdx+1],
		fc.encodedSequence[closingThreePrimeIdx-1],
		fc.sequence[closingFivePrimeIdx-1:], fc.energyParams)

	var singleStrandedFivePrimeIdx, singleStrandedThreePrimeIdx int
	if nbUnpairedNucleotides > 0 {
		singleStrandedFivePrimeIdx = closingFivePrimeIdx + 1
		singleStrandedThreePrimeIdx = closingThreePrimeIdx - 1
	} else {
		// There is no single stranded region on the hairpin
		singleStrandedFivePrimeIdx = -1
		singleStrandedThreePrimeIdx = -1
	}

	return Hairpin{
		Stem:                        stem,
		SingleStrandedFivePrimeIdx:  singleStrandedFivePrimeIdx,
		SingleStrandedThreePrimeIdx: singleStrandedThreePrimeIdx,
		Energy:                      energy,
	}
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
func EvaluateHairpinLoop(size, basePairType, fivePrimeMismatch, threePrimeMismatch int, sequence string, energyParams *energy_params.EnergyParams) int {
	var energy int

	if size <= energy_params.MaxLenLoop {
		energy = energyParams.HairpinLoop[size]
	} else {
		energy = energyParams.HairpinLoop[energy_params.MaxLenLoop] + int(energyParams.LXC*math.Log(float64(size)/float64(energy_params.MaxLenLoop)))
	}

	if size < 3 {
		// should only be the case when folding alignments
		return energy
	}

	if size == 4 {
		tetraLoop := sequence[:6]
		tetraLoopEnergy, present := energyParams.TetraLoop[tetraLoop]
		if present {
			return tetraLoopEnergy
		}
	} else if size == 6 {
		hexaLoop := sequence[:8]
		hexaLoopEnergy, present := energyParams.HexaLoop[hexaLoop]
		if present {
			return hexaLoopEnergy
		}
	} else if size == 3 {
		triLoop := sequence[:5]
		triLoopEnergy, present := energyParams.TriLoop[triLoop]
		if present {
			return triLoopEnergy
		}

		if basePairType > 2 {
			// (X,Y) is not a GC or CG pair (could be GU, UG, AU, UA, or non-standard)
			return energy + energyParams.TerminalAUPenalty
		}

		// (X,Y) is a GC or CG pair
		return energy
	}

	energy += energyParams.MismatchHairpinLoop[basePairType][fivePrimeMismatch][threePrimeMismatch]

	return energy
}

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
 * multiLoop iterates through the multi-loop and finds each base pair
 * enclosed in the multi-loop. For each enclosed base pair, we add to the total
 * the energy due to that base pair and due to the substructure closed by the
 * enclosed base pair.
 * We also add a bonus energy based on the number of unpaired nucleotides in
 * the multi-loop (though it defaults to 0 right now).
 *
 * More information about multi-loops: https://rna.urmc.rochester.edu/NNDB/turner04/mb.html
 */
func multiLoop(fc *foldCompound, closingFivePrimeIdx int, stem Stem) MultiLoop {
	pairTable := fc.pairTable

	// the substructures present in the multi-loop
	var substructures []interface{}

	// energetic penalty imposed when a base pair encloses a multi loop
	multiLoopEnergy := fc.energyParams.MultiLoopClosingPenalty

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
		enclosedFivePrimeIdx++
		lenMultiLoopSingleStrandedRegion++
	}

	// add inital unpaired nucleotides
	nbUnpairedNucleotides := enclosedFivePrimeIdx - closingFivePrimeIdx - 1

	// the energy due to the substructures enclosed in the multi-loop
	// it is kept seperate so that we can distinguish between the energy
	// contribution due to the multi-loop vs energy contribution due to
	// substructure that branch out from this multi-loop
	var substructuresEnergy int = 0

	for enclosedFivePrimeIdx < closingThreePrimeIdx {

		// add up the contributions of the substructures of the multi loop
		en, substructure := evaluateLoop(fc, enclosedFivePrimeIdx)
		substructuresEnergy += en
		substructures = append(substructures, substructure)

		enclosedThreePrimeIdx := pairTable[enclosedFivePrimeIdx]
		enclosedPairType := encodedBasePairType(fc, enclosedFivePrimeIdx, enclosedThreePrimeIdx)

		switch fc.dangleModel {
		case NoDanglingEnds:
			multiLoopEnergy += EvaluateMultiLoopStem(enclosedPairType, -1, -1, fc.energyParams)

		case DoubleDanglingEnds:
			enclosedFivePrimeMismatch := fc.encodedSequence[enclosedFivePrimeIdx-1]
			enclosedThreePrimeMismatch := fc.encodedSequence[enclosedThreePrimeIdx+1]

			multiLoopEnergy += EvaluateMultiLoopStem(enclosedPairType, enclosedFivePrimeMismatch,
				enclosedThreePrimeMismatch, fc.energyParams)
		}
		// seek to the next stem
		enclosedFivePrimeIdx = enclosedThreePrimeIdx + 1

		for enclosedFivePrimeIdx < closingThreePrimeIdx && pairTable[enclosedFivePrimeIdx] == -1 {
			enclosedFivePrimeIdx++
			lenMultiLoopSingleStrandedRegion++
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

		// add unpaired nucleotides
		nbUnpairedNucleotides += enclosedFivePrimeIdx - enclosedThreePrimeIdx - 1
	}

	if closingFivePrimeIdx > 0 {
		// actual closing pair
		closingPairType := encodedBasePairType(fc, closingThreePrimeIdx, closingFivePrimeIdx)

		switch fc.dangleModel {
		case NoDanglingEnds:
			multiLoopEnergy += EvaluateMultiLoopStem(closingPairType, -1, -1, fc.energyParams)

		case DoubleDanglingEnds:
			closingFivePrimeMismatch := fc.encodedSequence[closingThreePrimeIdx-1]
			closingThreePrimeMismatch := fc.encodedSequence[closingFivePrimeIdx+1]

			multiLoopEnergy += EvaluateMultiLoopStem(closingPairType, closingFivePrimeMismatch,
				closingThreePrimeMismatch, fc.energyParams)
		}

	} else {
		// virtual closing pair
		multiLoopEnergy += EvaluateMultiLoopStem(0, -1, -1, fc.energyParams)
	}

	// add bonus energies for unpaired nucleotides
	multiLoopEnergy += nbUnpairedNucleotides * fc.energyParams.MultiLoopUnpairedNucelotideBonus

	multiLoop := MultiLoop{
		Stem:                stem,
		SubstructuresEnergy: substructuresEnergy,
		Substructures:       substructures,
		Energy:              multiLoopEnergy,
	}

	return multiLoop
}

/**
 *  Compute the energy contribution of a multi-loop stem.
 *
 *  Given a base pair (i,j) (encoded by `basePairType`), compute the
 *  energy contribution including dangling-end/terminal-mismatch contributions.
 *
 *  @param  basePairType          The encoded base pair type of (i, j) (see `basePairType()`)
 *  @param  fivePrimeMismatch     The encoded nucleotide directly adjacent (in the 5' direction) to i (may be -1 if index of i is 0)
 *  @param  threePrimeMismatch    The encoded nucleotide directly adjacent (in the 3' direction) to j (may be -1 if index of j is len(sequence) - 1)
 *  @param  energyParams          The pre-computed energy parameters
 *  @return                       The energy contribution of the introduced multi-loop stem
 */
func EvaluateMultiLoopStem(basePairType, fivePrimeMismatch, threePrimeMismatch int, energyParams *energy_params.EnergyParams) int {
	var energy int
	if fivePrimeMismatch >= 0 && threePrimeMismatch >= 0 {
		energy += energyParams.MismatchMultiLoop[basePairType][fivePrimeMismatch][threePrimeMismatch]
	} else if fivePrimeMismatch >= 0 {
		energy += energyParams.Dangle5[basePairType][fivePrimeMismatch]
	} else if threePrimeMismatch >= 0 {
		energy += energyParams.Dangle3[basePairType][threePrimeMismatch]
	}

	energy += energyParams.MultiLoopIntern[basePairType]

	if basePairType > 2 {
		// The closing base pair is not a GC or CG pair (could be GU, UG, AU, UA, or non-standard)
		energy += energyParams.TerminalAUPenalty
	}

	return energy
}

// Returns the minimum of two ints
func min(a, b int) int {
	if a < b {
		return a
	}
	return b
}

/**
PrintEnergyContributions pretty prints a list of energy contributions returned
from `MinimumFreeEnergy`.
*/
// func PrintEnergyContributions(energyContribution []StructuralMotifWithEnergy, sequence string) {
// 	for _, c := range energyContribution {
// 		structure := c.structure
// 		switch structure.loopType {
// 		case ExteriorLoop:
// 			fmt.Printf("External loop                           : %v\n", c.energy)
// 		case InteriorLoop:
// 			var i, j int = structure.closingFivePrimeIdx, structure.closingThreePrimeIdx
// 			var k, l int = structure.enclosedFivePrimeIdx, structure.enclosedThreePrimeIdx
// 			fmt.Printf("Interior loop (%v,%v) %v%v; (%v,%v) %v%v: %v\n",
// 				i, j,
// 				string(sequence[i]), string(sequence[j]),
// 				k, l,
// 				string(sequence[k]), string(sequence[l]),
// 				c.energy)
// 		case HairpinLoop:
// 			var i, j int = structure.closingFivePrimeIdx, structure.closingThreePrimeIdx
// 			fmt.Printf("Hairpin loop  (%v,%v) %v%v              : %v\n",
// 				i, j,
// 				string(sequence[i]), string(sequence[j]),
// 				c.energy)
// 		case MultiLoop:
// 			var i, j int = structure.closingFivePrimeIdx, structure.closingThreePrimeIdx
// 			fmt.Printf("Multi   loop (%v,%v) %v%v              : %v\n",
// 				i, j,
// 				string(sequence[i]), string(sequence[j]),
// 				c.energy)
// 		}
// 	}
// }
