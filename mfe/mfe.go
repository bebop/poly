package mfe

import (
	"errors"
	"fmt"
	"log"
	"math"
	"regexp"
	"strings"
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
	and [thier paper](https://almob.biomedcentral.com/articles/10.1186/s13015-016-0070-z)).
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
	// Arbitrary values to denote different types of loops. Used with the
	// `EnergyContribution` struct
	ExteriorLoop = iota
	InteriorLoop
	HairpinLoop
	MultiLoop

	/**
	* The number of distinguishable base pairs:
	* CG, GC, GU, UG, AU, UA, & non-standard (see `energy_params.md`)
	 */
	nbPairs int = 7
	// The maximum loop length
	maxLoop int = 30
	// Default temperature for free energy evaluation
	DefaultTemperature float64 = 37.0
)

// foldCompound holds all information needed to compute the free energy of a
// RNA secondary structure (a RNA sequence along with its folded structure).
type foldCompound struct {
	length              int                   // length of `sequence`
	energyParams        *energyParams         // The precomputed free energy contributions for each type of loop
	sequence            string                // The input sequence string
	encodedSequence     []int                 // Numerical encoding of the sequence (see `encodeSequence()` for more information)
	pairTable           []int                 // (see `pairTable()`)
	basePairEncodedType map[byte]map[byte]int // (see `basePairEncodedTypeMap()`)
}

/**
* EnergyContribution contains the energy contribution of a loop
* The fields `loopType` and `energy` will always have a value for all loop types.
* Except for `loopType == ExteriorLoop`, all loops will have a
* `closingFivePrimeIdx` and `closingThreePrimeIdx`.
* Only loops with `loopType == InteriorLoop` will have `enclosedFivePrimeIdx`
* and `enclosedThreePrimeIdx`.
 */
type EnergyContribution struct {
	closingFivePrimeIdx, closingThreePrimeIdx   int // The closing base pair is the base pair that closes the loop
	loopType                                    int // The type of loop (ExteriorLoop, InteriorLoop, HairpinLoop, or MultiLoop)
	energy                                      int // The free energy contribution of the loop in dcal / mol
	enclosedFivePrimeIdx, enclosedThreePrimeIdx int // The base pair (that is enclosed by a closing base pair) which delimits an interior loop (only available for interior loops)
}

/**
 *  MinimumFreeEnergy returns the free energy of an already folded RNA and a list of
 *  energy contribution of each loop.
 *
 *  @param sequence         A RNA sequence
 *  @param structure        Secondary structure in dot-bracket notation
 *  @param temperature      Temperature at which to evaluate the free energy of the structure
 *  @return                 The free energy of the input structure given the input sequence in kcal/mol,
														and a slice of `EnergyContribution` (which gives information about the
														energy contribution of each loop present in the secondary structure)
*/
func MinimumFreeEnergy(sequence, structure string, temperature float64) (float64, []EnergyContribution, error) {
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
		energyParams:        scaleEnergyParams(temperature),
		sequence:            sequence,
		encodedSequence:     encodeSequence(sequence),
		pairTable:           pairTable,
		basePairEncodedType: basePairEncodedTypeMap(),
	}

	energyInt, energyContributions := evaluateFoldCompound(fc)
	energy := float64(energyInt) / 100.0

	return energy, energyContributions, nil
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

// Encodes a sequence into its numerical representaiton
func encodeSequence(sequence string) []int {
	// Make a nucleotide byte -> int map
	nucleotideRuneMap := map[byte]int{
		'A': 1,
		'C': 2,
		'G': 3,
		'U': 4,
	}

	lenSequence := len(sequence)
	encodedSequence := make([]int, lenSequence)

	// encode the sequence based on nucleotideRuneMap
	for i := 0; i < lenSequence; i++ {
		encodedSequence[i] = nucleotideRuneMap[sequence[i]]
	}

	return encodedSequence
}

/**
* Returns a map that encodes a base pair to its numerical representation.
* The map is:
*    _  A  C  G  U
*	_ {0, 0, 0, 0, 0}
*	A {0, 0, 0, 0, 5}
*	C {0, 0, 0, 1, 0}
*	G {0, 0, 2, 0, 3}
*	U {0, 6, 0, 4, 0}
* For example, map['A']['U'] is 5.
* The encoded numerical representation of a base pair carries no meaning in
* itself except for where to find the relevent energy contributions of the base
* pair in the matrices of the energy paramaters (found in `energy_params.go`).
* Thus, any change to this map must be reflected in the energy
* parameters matrices in `energy_params.go`, and vice versa (see
* `energy_params.md` for more information.)
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

// Contains all the energy parameters needed for the free energy calculations
type energyParams struct {
	stackingPair [nbPairs + 1][nbPairs + 1]int
	hairpinLoop  [31]int
	bulge        [maxLoop + 1]int
	/* This field was previous called internal_loop, but has been renamed to
	interior loop to remain consistent with the names of other interior loops */
	interiorLoop                                       [maxLoop + 1]int
	mismatchExteriorLoop                               [nbPairs + 1][5][5]int
	mismatchInteriorLoop                               [nbPairs + 1][5][5]int
	mismatch1xnInteriorLoop                            [nbPairs + 1][5][5]int
	mismatch2x3InteriorLoop                            [nbPairs + 1][5][5]int
	mismatchHairpinLoop                                [nbPairs + 1][5][5]int
	mismatchMultiLoop                                  [nbPairs + 1][5][5]int
	dangle5                                            [nbPairs + 1][5]int
	dangle3                                            [nbPairs + 1][5]int
	interior1x1Loop                                    [nbPairs + 1][nbPairs + 1][5][5]int
	interior2x1Loop                                    [nbPairs + 1][nbPairs + 1][5][5][5]int
	interior2x2Loop                                    [nbPairs + 1][nbPairs + 1][5][5][5][5]int
	ninio                                              [5]int
	lxc                                                float64
	multiLoopBase, multiLoopClosingPenalty, terminalAU int
	multiLoopIntern                                    [nbPairs + 1]int
	tetraloops                                         string
	tetraloop                                          [200]int
	triloops                                           string
	triloop                                            [40]int
	hexaloops                                          string
	hexaloop                                           [40]int
	tripleC, multipleCA, multipleCB                    int
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

	for i := 0; i < len(structure); i++ {
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

func evaluateFoldCompound(fc *foldCompound) (int, []EnergyContribution) {
	energyContributions := make([]EnergyContribution, 0)

	energy, contribution := exteriorLoopEnergy(fc)
	energyContributions = append(energyContributions, contribution)

	pairTable := fc.pairTable
	for i := 0; i < fc.length; i++ {
		if pairTable[i] == -1 {
			continue
		}

		en, contributions := stackEnergy(fc, i)
		energy += en
		energyContributions = append(energyContributions, contributions...)
		i = pairTable[i]
	}

	return energy, energyContributions
}

/**
* (see `basePairEncodedTypeMap()`)
 */
func encodedBasePairType(fc *foldCompound, basePairFivePrime, basePairThreePrime int) int {
	return fc.basePairEncodedType[fc.sequence[basePairFivePrime]][fc.sequence[basePairThreePrime]]
}

/**
* Calculate the energy contribution of stabilizing dangling-ends/mismatches
* for all stems branching off the exterior loop.
* For example, if the structure is
* Structure: " .  .  (  (  (  (  .  .  .  )  )  )  )  .  .  .  (  (  .  .  .  .  .  .  .  .  )  )  .  ."
* Index:       0  1  2  3  4  5  6  7  8  9 10 11 12 13 14 15 16 17 18 19 20 21 22 23 24 25 26 27 28 29
* then the exterior loops in the structure are as follows: (2, 12) and (16, 27).
 */
func exteriorLoopEnergy(fc *foldCompound) (int, EnergyContribution) {
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

		energy += exteriorStemEnergy(basePairType, fivePrimeMismatch, threePrimeMismatch, fc.energyParams)

		// seek to the next stem
		pairFivePrimeIdx = pairThreePrimeIdx + 1
		for pairFivePrimeIdx < length && pairTable[pairFivePrimeIdx] == -1 {
			pairFivePrimeIdx++
		}
	}

	return energy, EnergyContribution{energy: energy}
}

/**
 *  Evaluate a stem branching off the exterior loop.
 *
 *  Given a base pair (i,j) (encoded by `basePairType`), compute the
 *  energy contribution including dangling-end/terminal-mismatch contributions.
 *
 *  You can prevent taking 5'-, 3'-dangles or mismatch contributions into
 *  account by passing -1 for `fivePrimeMismatch` and/or `threePrimeMismatch`
 *  respectively.
 *
 *  @param  basePairType   The encoded base pair type of (i, j) (see `basePairType()`)
 *  @param  fivePrimeMismatch     The encoded nucleotide directly adjacent (in the 5' direction) to i (may be -1 if index of i is 0)
 *  @param  threePrimeMismatch    The encoded nucleotide directly adjacent (in the 3' direction) to j (may be -1 if index of j is len(sequence) - 1)
 *  @param  energyParams          The pre-computed energy parameters
 *  @return                       The energy contribution of the introduced exterior-loop stem
 */
func exteriorStemEnergy(basePairType int, fivePrimeMismatch int, threePrimeMismatch int, energyParams *energyParams) int {
	var energy int = 0

	if fivePrimeMismatch >= 0 && threePrimeMismatch >= 0 {
		energy += energyParams.mismatchExteriorLoop[basePairType][fivePrimeMismatch][threePrimeMismatch]
	} else if fivePrimeMismatch >= 0 {
		// `j` is the last nucleotide of the sequence
		energy += energyParams.dangle5[basePairType][fivePrimeMismatch]
	} else if threePrimeMismatch >= 0 {
		// `i` is the first nucleotide of the sequence
		energy += energyParams.dangle3[basePairType][threePrimeMismatch]
	}

	if basePairType > 2 {
		// The closing base pair is not a GC or CG pair (could be GU, UG, AU, UA, or non-standard)
		energy += energyParams.terminalAU
	}

	return energy
}

/**
 * stackEnergy recursively calculate energy of substructure enclosed by
 * (closingFivePrimeIdx, closingThreePrimeIdx).
 */
func stackEnergy(fc *foldCompound, closingFivePrimeIdx int) (int, []EnergyContribution) {
	// recursively calculate energy of substructure enclosed by (pairFivePrimeIdx, pairThreePrimeIdx)
	energyContributions := make([]EnergyContribution, 0)

	pairTable := fc.pairTable
	energy := 0
	closingThreePrimeIdx := pairTable[closingFivePrimeIdx]

	if encodedBasePairType(fc, closingFivePrimeIdx, closingThreePrimeIdx) == 0 {
		panic(fmt.Sprintf("bases %v and %v (%v%v) can't pair!",
			closingFivePrimeIdx, closingThreePrimeIdx,
			string(fc.sequence[closingFivePrimeIdx]),
			string(fc.sequence[closingThreePrimeIdx])))
	}

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

			if encodedBasePairType(fc, enclosedThreePrimeIdx, enclosedFivePrimeIdx) == 0 {
				panic(fmt.Sprintf("bases %d and %d (%c%c) can't pair!",
					enclosedFivePrimeIdx, enclosedThreePrimeIdx,
					fc.sequence[enclosedThreePrimeIdx],
					fc.sequence[enclosedFivePrimeIdx]))
			}

			en, contribution := stackBulgeInteriorLoopEnergy(fc,
				closingFivePrimeIdx, closingThreePrimeIdx,
				enclosedFivePrimeIdx, enclosedThreePrimeIdx,
			)
			energy += en
			energyContributions = append(energyContributions, contribution)

			closingFivePrimeIdx = enclosedFivePrimeIdx
			closingThreePrimeIdx = enclosedThreePrimeIdx
		}
	} /* end for */

	if enclosedFivePrimeIdx > enclosedThreePrimeIdx {
		// hairpin
		en, hairpinContribution := hairpinLoopEnergy(fc, closingFivePrimeIdx, closingThreePrimeIdx)
		energy += en
		energyContributions = append(energyContributions, hairpinContribution)
		return energy, energyContributions
	} else {
		// we have a multi-loop

		en, multiLoopContributions := multiLoopEnergy(fc, closingFivePrimeIdx)
		energy += en
		energyContributions = append(energyContributions, multiLoopContributions...)
	}

	return energy, energyContributions
}

/**
 *  Evaluate the free energy contribution of an interior loop with delimiting
 *  base pairs (i,j) and (k,l). See `evaluateInteriorLoop()` for more details.
 */
func stackBulgeInteriorLoopEnergy(fc *foldCompound,
	closingFivePrimeIdx, closingThreePrimeIdx,
	enclosedFivePrimeIdx, enclosedThreePrimeIdx int) (int, EnergyContribution) {

	nbUnpairedFivePrime := enclosedFivePrimeIdx - closingFivePrimeIdx - 1
	nbUnpairedThreePrime := closingThreePrimeIdx - enclosedThreePrimeIdx - 1

	closingBasePairType := encodedBasePairType(fc, closingFivePrimeIdx, closingThreePrimeIdx)
	enclosedBasePairType := encodedBasePairType(fc, enclosedThreePrimeIdx, enclosedFivePrimeIdx)

	energy := evaluateStackBulgeInteriorLoop(nbUnpairedFivePrime, nbUnpairedThreePrime,
		closingBasePairType, enclosedBasePairType,
		fc.encodedSequence[closingFivePrimeIdx+1], fc.encodedSequence[closingThreePrimeIdx-1],
		fc.encodedSequence[enclosedFivePrimeIdx-1], fc.encodedSequence[enclosedThreePrimeIdx+1], fc.energyParams)

	energyContribution := EnergyContribution{
		closingFivePrimeIdx:   closingFivePrimeIdx,
		closingThreePrimeIdx:  closingThreePrimeIdx,
		enclosedFivePrimeIdx:  enclosedFivePrimeIdx,
		enclosedThreePrimeIdx: enclosedThreePrimeIdx,
		energy:                energy,
		loopType:              InteriorLoop,
	}
	return energy, energyContribution
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
 *  must be 'turned around' when evaluating the free energy of the interior loop
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
func evaluateStackBulgeInteriorLoop(nbUnpairedLeftLoop, nbUnpairedRightLoop int,
	closingBasePairType, enclosedBasePairType int,
	closingFivePrimeMismatch, closingThreePrimeMismatch,
	enclosedThreePrimeMismatch, enclosedFivePrimeMismatch int,
	energyParams *energyParams) int {
	/* compute energy of degree 2 loop (stack bulge or interior) */
	var nbUnpairedLarger, nbUnpairedSmaller int

	if nbUnpairedLeftLoop > nbUnpairedRightLoop {
		nbUnpairedLarger = nbUnpairedLeftLoop
		nbUnpairedSmaller = nbUnpairedRightLoop
	} else {
		nbUnpairedLarger = nbUnpairedRightLoop
		nbUnpairedSmaller = nbUnpairedLeftLoop
	}

	if nbUnpairedLarger == 0 {
		// stacking pair
		return energyParams.stackingPair[closingBasePairType][enclosedBasePairType]
	}

	if nbUnpairedSmaller == 0 {
		// bulge
		var energy int

		if nbUnpairedLarger <= maxLoop {
			energy = energyParams.bulge[nbUnpairedLarger]
		} else {
			energy = energyParams.bulge[30] + int(energyParams.lxc*math.Log(float64(nbUnpairedLarger)/30.0))
		}

		if nbUnpairedLarger == 1 {
			energy += energyParams.stackingPair[closingBasePairType][enclosedBasePairType]
		} else {
			if closingBasePairType > 2 {
				// The closing base pair is not a GC or CG pair (could be GU, UG, AU, UA, or non-standard)
				energy += energyParams.terminalAU
			}

			if enclosedBasePairType > 2 {
				// The encosed base pair is not a GC or CG pair (could be GU, UG, AU, UA, or non-standard)
				energy += energyParams.terminalAU
			}
		}

		return energy
	} else {
		// interior loop
		if nbUnpairedSmaller == 1 {
			if nbUnpairedLarger == 1 {
				// 1x1 loop
				return energyParams.interior1x1Loop[closingBasePairType][enclosedBasePairType][closingFivePrimeMismatch][closingThreePrimeMismatch]
			}

			if nbUnpairedLarger == 2 {
				// 2x1 loop
				// vivek: shouldn't it be a 1x2 loop (following the convention that the first integer corresponds to `nbUnpairedSmaller`)? No consistency in naming here.
				if nbUnpairedLeftLoop == 1 {
					return energyParams.interior2x1Loop[closingBasePairType][enclosedBasePairType][closingFivePrimeMismatch][enclosedFivePrimeMismatch][closingThreePrimeMismatch]
				} else {
					return energyParams.interior2x1Loop[enclosedBasePairType][closingBasePairType][enclosedFivePrimeMismatch][closingFivePrimeMismatch][enclosedThreePrimeMismatch]
				}
			} else {
				// 1xn loop

				var energy int
				// vivek: Why do we add 1 here?
				if nbUnpairedLarger+1 <= maxLoop {
					energy = energyParams.interiorLoop[nbUnpairedLarger+1]
				} else {
					energy = energyParams.interiorLoop[30] + int(energyParams.lxc*math.Log((float64(nbUnpairedLarger)+1.0)/30.0))
				}
				energy += min(MAX_NINIO, (nbUnpairedLarger-nbUnpairedSmaller)*energyParams.ninio[2])
				energy += energyParams.mismatch1xnInteriorLoop[closingBasePairType][closingFivePrimeMismatch][closingThreePrimeMismatch] + energyParams.mismatch1xnInteriorLoop[enclosedBasePairType][enclosedFivePrimeMismatch][enclosedThreePrimeMismatch]
				return energy
			}
		} else if nbUnpairedSmaller == 2 {
			if nbUnpairedLarger == 2 {
				// 2x2 loop
				return energyParams.interior2x2Loop[closingBasePairType][enclosedBasePairType][closingFivePrimeMismatch][enclosedThreePrimeMismatch][enclosedFivePrimeMismatch][closingThreePrimeMismatch]
			} else if nbUnpairedLarger == 3 {
				// 2x3 loop

				var energy int
				energy = energyParams.interiorLoop[5] + energyParams.ninio[2]
				energy += energyParams.mismatch2x3InteriorLoop[closingBasePairType][closingFivePrimeMismatch][closingThreePrimeMismatch] + energyParams.mismatch2x3InteriorLoop[enclosedBasePairType][enclosedFivePrimeMismatch][enclosedThreePrimeMismatch]
				return energy
			}
		}

		{
			/* generic interior loop (no else here!)*/
			nbUnpairedNucleotides := nbUnpairedLarger + nbUnpairedSmaller

			var energy int
			if nbUnpairedNucleotides <= maxLoop {
				energy = energyParams.interiorLoop[nbUnpairedNucleotides]
			} else {
				energy = energyParams.interiorLoop[30] + int(energyParams.lxc*math.Log(float64(nbUnpairedNucleotides)/30.0))
			}

			energy += min(MAX_NINIO, (nbUnpairedLarger-nbUnpairedSmaller)*energyParams.ninio[2])

			energy += energyParams.mismatchInteriorLoop[closingBasePairType][closingFivePrimeMismatch][closingThreePrimeMismatch] + energyParams.mismatchInteriorLoop[enclosedBasePairType][enclosedFivePrimeMismatch][enclosedThreePrimeMismatch]
			return energy
		}
	}
}

/**
 *  Evaluate free energy of a hairpin loop encosed by the base pair
 *  (pairFivePrimeIdx, pairThreePrimeIdx).
 *
 *  @param  fc                 The foldCompund for the particular energy evaluation
 *  @param  pairFivePrimeIdx   5'-position of the base pair enclosing the hairpin loop
 *  @param  pairThreePrimeIdx  3'-position of the base pair enclosing the hairpin loop
 *  @returns                   Free energy of the hairpin loop closed by (pairFivePrimeIdx,pairThreePrimeIdx) in dcal/mol
 */
func hairpinLoopEnergy(fc *foldCompound, pairFivePrimeIdx, pairThreePrimeIdx int) (int, EnergyContribution) {
	nbUnpairedNucleotides := pairThreePrimeIdx - pairFivePrimeIdx - 1 // also the size of the hairpin loop
	basePairType := encodedBasePairType(fc, pairFivePrimeIdx, pairThreePrimeIdx)

	energy := evaluateHairpinLoop(nbUnpairedNucleotides, basePairType,
		fc.encodedSequence[pairFivePrimeIdx+1],
		fc.encodedSequence[pairThreePrimeIdx-1],
		fc.sequence[pairFivePrimeIdx-1:], fc.energyParams)

	energyContribution := EnergyContribution{
		closingFivePrimeIdx:  pairFivePrimeIdx,
		closingThreePrimeIdx: pairThreePrimeIdx,
		energy:               energy,
		loopType:             HairpinLoop,
	}
	return energy, energyContribution
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
 *  @param  size                 The size of the hairpin loop (number of unpaired nucleotides)
 *  @param  basePairType         The pair type of the base pair closing the hairpin
 *  @param  fivePrimeMismatch    The 5'-mismatching encoded nucleotide
 *  @param  threePrimeMismatch   The 3'-mismatching encoded nucleotide
 *  @param  sequence						 The sequence of the loop
 *  @param  energyParams     		 The datastructure containing scaled energy parameters
 *  @return The free energy of the hairpin loop in dcal/mol
 */
func evaluateHairpinLoop(size, basePairType, fivePrimeMismatch, threePrimeMismatch int, sequence string, energyParams *energyParams) int {
	var energy int

	if size <= 30 {
		energy = energyParams.hairpinLoop[size]
	} else {
		energy = energyParams.hairpinLoop[30] + int(energyParams.lxc*math.Log(float64(size)/30.0))
	}

	if size < 3 {
		// should only be the case when folding alignments
		return energy
	}

	if size == 4 {
		tetraloop := sequence[:6]
		idx := strings.Index(string(energyParams.tetraloops), tetraloop)
		if idx != -1 {
			return energyParams.tetraloop[idx/7]
		}
	} else if size == 6 {
		hexaloop := sequence[:8]
		idx := strings.Index(string(energyParams.hexaloops), hexaloop)
		if idx != -1 {
			return energyParams.hexaloop[idx/9]
		}
	} else if size == 3 {
		triloop := sequence[:5]
		idx := strings.Index(string(energyParams.triloops), triloop)
		if idx != -1 {
			return energyParams.triloop[idx/6]
		}

		if basePairType > 2 {
			// (X,Y) is not a GC or CG pair (could be GU, UG, AU, UA, or non-standard)
			return energy + energyParams.terminalAU
		}

		// (X,Y) is a GC or CG pair
		return energy
	}

	// size is 5 or greater than 6
	energy += energyParams.mismatchHairpinLoop[basePairType][fivePrimeMismatch][threePrimeMismatch]

	return energy
}

/**
 *  Compute the energy of a multi-loop.
 *  A multi-loop has this structure:
 *							·   ·
 * 							·   ·
 * 							·   ·
 * 							A - B
 * 						·		 	 ·
 * 					 ·				C · · ·
 * 					·					|
 * 		5' - X					D · · ·
 * 				 |				 ·
 * 		3' - Y - V - U
 * 						 ·   ·
 * 						 ·   ·
 * 						 ·   ·
 * where X-Y marks the closing pair (X is the nucleotide at `closingFivePrimeIdx`)
 * of the multi-loop, and `·`s are an arbitrary number of unparied nucelotides
 * (can also be 0). (A,B), (C,D), and (V,U) are the base pairs enclosed in this
 * multi-loop (vivek: I think there can be an arbitraty number of enclosed
 * base pairs. I'm not sure if there's a limit on the least or greatest numbers
 * of base pairs that can be enclosed).
 *
 * multiLoopEnergy iterates through the multi-loop and finds each base pair
 * enclosed in the multi-loop. For each enclosed base pair, we add to the total
 * the energy due to that base pair and due to the substructure closed by the
 * enclosed base pair.
 * We also add a bonus energy based on the number of unpaired nucleotides in
 * the multi-loop (though it defaults to 0 right now).
 */
func multiLoopEnergy(fc *foldCompound, closingFivePrimeIdx int) (int, []EnergyContribution) {
	energyContributions := make([]EnergyContribution, 0)
	pairTable := fc.pairTable

	// energetic penalty imposed when a base pair encloses a multi loop
	multiLoopEnergy := fc.energyParams.multiLoopClosingPenalty

	// (pairFivePrimeIdx, pairThreePrimeIdx) is exterior pair of multi loop
	if closingFivePrimeIdx >= pairTable[closingFivePrimeIdx] {
		panic("multiLoopEnergy: pairFivePrimeIdx is not 5' base of a pair that closes a loop")
	}

	closingThreePrimeIdx := pairTable[closingFivePrimeIdx]

	// add energy due to closing pair of multi-loop
	// vivek: Is this a bug? Why are the five and three prime mismatches opposite?
	// Created an issue in the original repo: https://github.com/ViennaRNA/ViennaRNA/issues/126
	basePairType := encodedBasePairType(fc, closingThreePrimeIdx, closingFivePrimeIdx)

	closingFivePrimeMismatch := fc.encodedSequence[closingThreePrimeIdx-1]
	closingThreePrimeMismatch := fc.encodedSequence[closingFivePrimeIdx+1]
	multiLoopEnergy += multiLoopStemEnergy(basePairType, closingFivePrimeMismatch,
		closingThreePrimeMismatch, fc.energyParams)

	// The index of the enclosed base pair at the 5' end
	enclosedFivePrimeIdx := closingFivePrimeIdx + 1

	// seek to the first stem (i.e. the first enclosed base pair)
	for enclosedFivePrimeIdx <= closingThreePrimeIdx && pairTable[enclosedFivePrimeIdx] == -1 {
		enclosedFivePrimeIdx++
	}

	// add inital unpaired nucleotides
	nbUnpairedNucleotides := enclosedFivePrimeIdx - closingFivePrimeIdx - 1

	// the energy due to the substructures enclosed in the multi-loop
	// it is kept seperate so that we can distinguish between the energy
	// contribution due to the multi-loop vs energy contribution due to
	// substructure that branch out from this multi-loop
	substructuresEnergy := 0

	// iterate through structure and find all enclosed base pairs
	for enclosedFivePrimeIdx < closingThreePrimeIdx {

		// add up the contributions of the substructures of the multi loop
		en, substructureContributions := stackEnergy(fc, enclosedFivePrimeIdx)
		substructuresEnergy += en
		energyContributions = append(energyContributions, substructureContributions...)

		/* enclosedFivePrimeIdx must have a pairing partner */
		enclosedThreePrimeIdx := pairTable[enclosedFivePrimeIdx]

		/* get type of the enclosed base pair (enclosedFivePrimeIdx,enclosedThreePrimeIdx) */
		basePairType := encodedBasePairType(fc, enclosedFivePrimeIdx, enclosedThreePrimeIdx)

		enclosedFivePrimeMismatch := fc.encodedSequence[enclosedFivePrimeIdx-1]
		enclosedThreePrimeMismatch := fc.encodedSequence[enclosedThreePrimeIdx+1]

		multiLoopEnergy += multiLoopStemEnergy(basePairType, enclosedFivePrimeMismatch,
			enclosedThreePrimeMismatch, fc.energyParams)

		// seek to the next stem
		enclosedFivePrimeIdx = enclosedThreePrimeIdx + 1

		for enclosedFivePrimeIdx < closingThreePrimeIdx && pairTable[enclosedFivePrimeIdx] == -1 {
			enclosedFivePrimeIdx++
		}

		// add unpaired nucleotides
		nbUnpairedNucleotides += enclosedFivePrimeIdx - enclosedThreePrimeIdx - 1
	}

	// add bonus energies for unpaired nucleotides
	multiLoopEnergy += nbUnpairedNucleotides * fc.energyParams.multiLoopBase

	// Add a energy contribution for this multi-loop
	contribution := EnergyContribution{
		closingFivePrimeIdx:  closingFivePrimeIdx,
		closingThreePrimeIdx: closingThreePrimeIdx,
		energy:               multiLoopEnergy,
		loopType:             MultiLoop,
	}
	energyContributions = append(energyContributions, contribution)
	return multiLoopEnergy + substructuresEnergy, energyContributions
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
func multiLoopStemEnergy(basePairType, fivePrimeMismatch, threePrimeMismatch int, energyParams *energyParams) int {
	var energy int = energyParams.mismatchMultiLoop[basePairType][fivePrimeMismatch][threePrimeMismatch] +
		energyParams.multiLoopIntern[basePairType]

	if basePairType > 2 {
		// The closing base pair is not a GC or CG pair (could be GU, UG, AU, UA, or non-standard)
		energy += energyParams.terminalAU
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

func scaleEnergyParams(temperature float64) *energyParams {
	var tempf int = int((temperature + K0) / Tmeasure)

	var params *energyParams = &energyParams{
		lxc:                     lxc37 * float64(tempf),
		tripleC:                 rescaleDg(TripleC37, TripleCdH, tempf),
		multipleCA:              rescaleDg(MultipleCA37, MultipleCAdH, tempf),
		multipleCB:              rescaleDg(TerminalAU37, TerminalAUdH, tempf),
		terminalAU:              rescaleDg(TerminalAU37, TerminalAUdH, tempf),
		multiLoopBase:           rescaleDg(ML_BASE37, ML_BASEdH, tempf),
		multiLoopClosingPenalty: rescaleDg(ML_closing37, ML_closingdH, tempf),
	}

	params.ninio[2] = rescaleDg(ninio37, niniodH, tempf)

	for i := 0; i < 31; i++ {
		params.hairpinLoop[i] = rescaleDg(hairpin37[i], hairpindH[i], tempf)
	}

	var i int
	for i = 0; i <= min(30, maxLoop); i++ {
		params.bulge[i] = rescaleDg(bulge37[i], bulgedH[i], tempf)
		params.interiorLoop[i] = rescaleDg(internal_loop37[i], internal_loopdH[i], tempf)
	}

	for ; i <= maxLoop; i++ {
		params.bulge[i] = params.bulge[30] +
			int(params.lxc*math.Log(float64(i)/30.0))
		params.interiorLoop[i] = params.interiorLoop[30] +
			int(params.lxc*math.Log(float64(i)/30.0))
	}

	// This could be 281 or only the amount of chars
	for i := 0; (i * 7) < len(Tetraloops); i++ {
		params.tetraloop[i] = rescaleDg(Tetraloop37[i], TetraloopdH[i], tempf)
	}

	for i := 0; (i * 5) < len(Triloops); i++ {
		params.triloop[i] = rescaleDg(Triloop37[i], TriloopdH[i], tempf)
	}

	for i := 0; (i * 9) < len(Hexaloops); i++ {
		params.hexaloop[i] = rescaleDg(Hexaloop37[i], HexaloopdH[i], tempf)
	}

	for i := 0; i <= nbPairs; i++ {
		params.multiLoopIntern[i] = rescaleDg(ML_intern37, ML_interndH, tempf)
	}

	/* stacks    G(T) = H - [H - G(T0)]*T/T0 */
	for i := 0; i <= nbPairs; i++ {
		for j := 0; j <= nbPairs; j++ {
			params.stackingPair[i][j] = rescaleDg(stack37[i][j],
				stackdH[i][j],
				tempf)
		}
	}

	/* mismatches */
	for i := 0; i <= nbPairs; i++ {
		for j := 0; j < 5; j++ {
			for k := 0; k < 5; k++ {
				var mm int
				params.mismatchInteriorLoop[i][j][k] = rescaleDg(mismatchI37[i][j][k],
					mismatchIdH[i][j][k],
					tempf)
				params.mismatchHairpinLoop[i][j][k] = rescaleDg(mismatchH37[i][j][k],
					mismatchHdH[i][j][k],
					tempf)
				params.mismatch1xnInteriorLoop[i][j][k] = rescaleDg(mismatch1nI37[i][j][k],
					mismatch1nIdH[i][j][k],
					tempf)
				params.mismatch2x3InteriorLoop[i][j][k] = rescaleDg(mismatch23I37[i][j][k],
					mismatch23IdH[i][j][k],
					tempf)
				// if model_details.dangles > 0 {
				mm = rescaleDg(mismatchM37[i][j][k],
					mismatchMdH[i][j][k],
					tempf)
				if mm > 0 {
					params.mismatchMultiLoop[i][j][k] = 0
				} else {
					params.mismatchMultiLoop[i][j][k] = mm
				}

				mm = rescaleDg(mismatchExt37[i][j][k],
					mismatchExtdH[i][j][k],
					tempf)
				if mm > 0 {
					params.mismatchExteriorLoop[i][j][k] = 0
				} else {
					params.mismatchExteriorLoop[i][j][k] = mm
				}
				// } else {
				// 	params.mismatchExt[i][j][k] = 0
				// 	params.mismatchM[i][j][k] = 0
				// }
			}
		}
	}

	/* dangles */
	for i := 0; i <= nbPairs; i++ {
		for j := 0; j < 5; j++ {
			var dd int
			dd = rescaleDg(dangle5_37[i][j],
				dangle5_dH[i][j],
				tempf)
			if dd > 0 {
				params.dangle5[i][j] = 0
			} else {
				params.dangle5[i][j] = dd
			}

			dd = rescaleDg(dangle3_37[i][j],
				dangle3_dH[i][j],
				tempf)
			if dd > 0 {
				params.dangle3[i][j] = 0
			} else {
				params.dangle3[i][j] = dd
			}
		}
	}

	/* interior 1x1 loops */
	for i := 0; i <= nbPairs; i++ {
		for j := 0; j <= nbPairs; j++ {
			for k := 0; k < 5; k++ {
				for l := 0; l < 5; l++ {
					params.interior1x1Loop[i][j][k][l] = rescaleDg(int11_37[i][j][k][l],
						int11_dH[i][j][k][l],
						tempf)
				}
			}
		}
	}

	/* interior 2x1 loops */
	for i := 0; i <= nbPairs; i++ {
		for j := 0; j <= nbPairs; j++ {
			for k := 0; k < 5; k++ {
				for l := 0; l < 5; l++ {
					var m int
					for m = 0; m < 5; m++ {
						params.interior2x1Loop[i][j][k][l][m] = rescaleDg(int21_37[i][j][k][l][m],
							int21_dH[i][j][k][l][m],
							tempf)
					}
				}
			}
		}
	}

	/* interior 2x2 loops */
	for i := 0; i <= nbPairs; i++ {
		for j := 0; j <= nbPairs; j++ {
			for k := 0; k < 5; k++ {
				for l := 0; l < 5; l++ {
					var m, n int
					for m = 0; m < 5; m++ {
						for n = 0; n < 5; n++ {
							params.interior2x2Loop[i][j][k][l][m][n] = rescaleDg(int22_37[i][j][k][l][m][n],
								int22_dH[i][j][k][l][m][n],
								tempf)
						}
					}
				}
			}
		}
	}

	params.tetraloops = Tetraloops
	params.triloops = Triloops
	params.hexaloops = Hexaloops

	return params
}

func rescaleDg(dG, dH, dT int) int {
	return (dH - (dH-dG)*dT)
}

/**
* Logs energy contributions in the same format as ViennaRNA does. This is used
* to compare output to ViennaRNA in test cases.
* Note: 1 is added to the indexes of the closing and enclosed base pairs as
* everyting in ViennaRNA is 1-indexed, but output from the `MinimumFreeEnergy` func
* is 0-indexed.
 */
func logEnergyContributions(energyContribution []EnergyContribution, sequence string) {
	for _, c := range energyContribution {
		switch c.loopType {
		case ExteriorLoop:
			log.Printf("External loop                           : %v\n", c.energy)
		case InteriorLoop:
			var i, j int = c.closingFivePrimeIdx, c.closingThreePrimeIdx
			var k, l int = c.enclosedFivePrimeIdx, c.enclosedThreePrimeIdx
			log.Printf("Interior loop (%v,%v) %v%v; (%v,%v) %v%v: %v\n",
				i+1, j+1,
				string(sequence[i]), string(sequence[j]),
				k+1, l+1,
				string(sequence[k]), string(sequence[l]),
				c.energy)
		case HairpinLoop:
			var i, j int = c.closingFivePrimeIdx, c.closingThreePrimeIdx
			log.Printf("Hairpin loop  (%v,%v) %v%v              : %v\n",
				i+1, j+1,
				string(sequence[i]), string(sequence[j]),
				c.energy)
		case MultiLoop:
			var i, j int = c.closingFivePrimeIdx, c.closingThreePrimeIdx
			log.Printf("Multi   loop (%v,%v) %v%v              : %v\n",
				i+1, j+1,
				string(sequence[i]), string(sequence[j]),
				c.energy)
		}
	}
}
