package poly

import (
	"errors"
	"fmt"
	"log"
	"math"
	"strings"
)

/**
* The MFE calculation functionality has been taken directly from ViennaRNA with the follwing changes:
* - ViennaRNA includes the ability to specify a dangle model (more info available at [src/bin/RNAeval.ggo#L153](https://github.com/ViennaRNA/ViennaRNA/blob/d6fbaf2b33d65944f9249a91ed5ab4b3277f7d06/src/bin/RNAeval.ggo#L153)). This implementation keeps it simple and defaults to the value of -d2.
* - ViennaRNA includes the ability to specify hard and soft constraints (more info available from [their ppt explaining hard and soft constrains](http://benasque.org/2015rna/talks_contr/2713_lorenz_benasque_2015.pdf), [official docs](https://www.tbi.univie.ac.at/RNA/ViennaRNA/doc/html/group__constraints.html), [thier paper](https://almob.biomedcentral.com/articles/10.1186/s13015-016-0070-z)). This implementation keeps it simple and defaults to no hard or soft constraints.
* - ViennaRNA includes the ability to calculate the minimum free energy of co-folded sequences. This implementation keeps it simple and defaults to calculating the mfe of only single sequences.
 */

const (
	/**
	* The number of distinguishable base pairs:
	* CG, GC, GU, UG, AU, UA, & non-standard (see `energy_params.md`)
	 */
	nbPairs            int     = 7
	maxLoop            int     = 30   // The maximum loop length
	DefaultTemperature float64 = 37.0 // Default temperature for free energy evaluation
)

type foldCompound struct {
	length          int           // length of `sequence`
	energy_params   *energyParams // The precomputed free energy contributions for each type of loop
	sequence        string        // The input sequence string
	encodedSequence []int         // Numerical encoding of the sequence (see `encodeSequence()` for more information)
	pairTable       []int         // (see `pairTable()`)
}

// Struct to contain energy contributions of each loop
type EnergyContribution struct {
	closingFivePrimeIdx, closingThreePrimeIdx   int // The closing base pair is the base pair that closes the loop
	loopType                                    int // The type of loop (ExternalLoop, InteriorLoop, HairpinLoop, or MultiLoop)
	energy                                      int // The free energy contribution of the loop in dcal / mol
	enclosedFivePrimeIdx, enclosedThreePrimeIdx int // The base pair (that is enclosed by a closing base pair) which delimits an interior loop (only available for interior loops)
}

const (
	ExternalLoop = iota
	InteriorLoop
	HairpinLoop
	MultiLoop
)

/**
 *  Calculate the free energy of an already folded RNA. Contributions on a per-loop base are logged.
 *
 *  This function allows for detailed energy evaluation of a given sequence/structure pair.
 *
 *  @param sequence         A RNA sequence
 *  @param structure        Secondary structure in dot-bracket notation
 *  @return                 The free energy of the input structure given the input sequence in kcal/mol
 */
func CalculateMFE(sequence, structure string, temperature float64) (float64, []EnergyContribution, error) {
	lenSequence := len(sequence)
	lenStructure := len(structure)

	if lenSequence != lenStructure {
		return 0, nil, fmt.Errorf("length of sequence (%v) != Length of structure (%v)", lenSequence, lenStructure)
	} else if lenStructure == 0 {
		return 0, nil, errors.New("lengths of sequence and structure cannot be 0")
	}

	sequence = strings.ToUpper(sequence)
	// vivek: Should we convert DNA to RNA and should we ensure sequence is only
	// made up of ACGU?

	pairTable, err := pairTable(structure)
	if err != nil {
		return 0, nil, err
	}

	fc := &foldCompound{
		length:          lenSequence,
		energy_params:   scaleEnergyParams(temperature),
		sequence:        sequence,
		encodedSequence: encodeSequence(sequence),
		pairTable:       pairTable,
	}

	energyInt, energyContributions := evaluateFoldCompound(fc)
	energy := float64(energyInt) / 100.0

	return energy, energyContributions, nil
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

	len_sequence := len(sequence)
	encodedSequence := make([]int, len_sequence)

	// encode the sequence based on nucleotideRuneMap
	for i := 0; i < len_sequence; i++ {
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
	interiorLoop                              [maxLoop + 1]int
	mismatchExteriorLoop                      [nbPairs + 1][5][5]int
	mismatchInteriorLoop                      [nbPairs + 1][5][5]int
	mismatch1xnInteriorLoop                   [nbPairs + 1][5][5]int
	mismatch2x3InteriorLoop                   [nbPairs + 1][5][5]int
	mismatchHairpinLoop                       [nbPairs + 1][5][5]int
	mismatchMultiLoop                         [nbPairs + 1][5][5]int
	dangle5                                   [nbPairs + 1][5]int
	dangle3                                   [nbPairs + 1][5]int
	interior1x1Loop                           [nbPairs + 1][nbPairs + 1][5][5]int
	interior2x1Loop                           [nbPairs + 1][nbPairs + 1][5][5][5]int
	interior2x2Loop                           [nbPairs + 1][nbPairs + 1][5][5][5][5]int
	ninio                                     [5]int
	lxc                                       float64
	MLbase, MLclosing, TerminalAU, DuplexInit int
	MLintern                                  [nbPairs + 1]int
	Tetraloops                                string
	Tetraloop                                 [200]int
	Triloops                                  string
	Triloop                                   [40]int
	Hexaloops                                 string
	Hexaloop                                  [40]int
	TripleC, MultipleCA, MultipleCB           int
	basePairEncodedTypeMap                    map[byte]map[byte]int
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
	len_seq := len(structure)

	// vivek: this was a constraint from the original code, but not sure why it's
	// there. I don't think it should be a problem, so I've commented it out for
	// now.
	// if len_seq > SHRT_MAX {
	// 	err := fmt.Sprintf("pairTable: Structure too long to be converted to pair
	//  table (n=%v, max=%v)", len_seq, SHRT_MAX)
	// 	return nil, errors.New(err)
	// }

	pairTable = make([]int, len_seq)

	// the characters of the opening and closing brackets in the structure string
	var openBracket, closeBracket byte = '(', ')'

	// keeps track of the indexes of open brackets. Indexes of open brackets are
	// pushed onto stack and poped off when a closing bracket is encountered
	var openBracketIdxStack []int = make([]int, len_seq)
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
	energy, contribution := externalLoopEnergy(fc)
	energyContributions = append(energyContributions, contribution)

	pairTable := fc.pairTable
	for i := 0; i < fc.length; i++ {
		if pairTable[i] == -1 {
			continue
		}

		energy += stackEnergy(fc, i, &energyContributions)
		i = pairTable[i]
	}

	return energy, energyContributions
}

/*
*	returns the numerical representation of a base pair based on the map
* `basePairEncodedTypeMap`. 7 is a special value that signifies a non-standard
* base pair (see `energy_params.md` for more information).
* As stated in `basePairEncodedTypeMap()`, the numerical representation carries
* no meaning in itself except for where to find the relevent energy
* contributions of the base pair (i, j) in the matrices of the energy
* paramaters (found in `energy_params.go`).
 */
func encodedBasePairType(i, j byte, basePairEncodedTypeMap map[byte]map[byte]int) int {
	var encodedType int = basePairEncodedTypeMap[i][j]

	if encodedType == 0 {
		return 7
	} else {
		return encodedType
	}
}

/**
* Calculate the energy contribution of stabilizing dangling-ends/mismatches
* for all stems branching off the exterior loop
 */
func externalLoopEnergy(fc *foldCompound) (int, EnergyContribution) {
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
		basePairType := encodedBasePairType(fc.sequence[pairFivePrimeIdx],
			fc.sequence[pairThreePrimeIdx],
			fc.energy_params.basePairEncodedTypeMap)

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

		energy += exteriorStemEnergy(basePairType, fivePrimeMismatch, threePrimeMismatch, fc.energy_params)

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
 *  Instead of returning the energy contribution per-se, this function returns
 *  the corresponding Boltzmann factor.
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
		energy += energyParams.TerminalAU
	}

	return energy
}

// TODO: Document
func stackEnergy(fc *foldCompound, pairFivePrimeIdx int, energyContributions *([]EnergyContribution)) int {
	// recursively calculate energy of substructure enclosed by (pairFivePrimeIdx, pairThreePrimeIdx)
	pairTable := fc.pairTable
	energy := 0
	pairThreePrimeIdx := pairTable[pairFivePrimeIdx]

	if fc.energy_params.basePairEncodedTypeMap[fc.sequence[pairFivePrimeIdx]][fc.sequence[pairThreePrimeIdx]] == 0 {
		panic(fmt.Sprintf("bases %v and %v (%v%v) can't pair!",
			pairFivePrimeIdx, pairThreePrimeIdx,
			string(fc.sequence[pairFivePrimeIdx]),
			string(fc.sequence[pairThreePrimeIdx])))
	}

	// iterator from the 5' to 3' direction starting at `pairFivePrimeIdx`
	fivePrimeIter := pairFivePrimeIdx

	// iterator from the 3' to 5' direction starting at `pairThreePrimeIdx`
	threePrimeIter := pairThreePrimeIdx

	for fivePrimeIter < threePrimeIter {
		// process all stacks and interior loops

		// seek to opening pair from 5' end
		fivePrimeIter++
		for pairTable[fivePrimeIter] == -1 {
			fivePrimeIter++
		}

		// seek to closing pair from 3' end
		threePrimeIter--
		for pairTable[threePrimeIter] == -1 {
			threePrimeIter--
		}

		if pairTable[threePrimeIter] != fivePrimeIter || fivePrimeIter > threePrimeIter {
			break
		}

		// vivek: should this be basePairType[fivePrimeIter][threePrimeIter] or is the current order correct?
		// created an issue in the original repo: https://github.com/ViennaRNA/ViennaRNA/issues/125
		if fc.energy_params.basePairEncodedTypeMap[fc.sequence[threePrimeIter]][fc.sequence[fivePrimeIter]] == 0 {
			panic(fmt.Sprintf("bases %d and %d (%c%c) can't pair!",
				fivePrimeIter, threePrimeIter,
				fc.sequence[fivePrimeIter],
				fc.sequence[threePrimeIter]))
		}

		en, contribution := interiorLoopEnergy(fc, pairFivePrimeIdx, pairThreePrimeIdx, fivePrimeIter, threePrimeIter)
		energy += en
		*energyContributions = append(*energyContributions, contribution)

		pairFivePrimeIdx = fivePrimeIter
		pairThreePrimeIdx = threePrimeIter
	} /* end for */

	// fivePrimeIter & threePrimeIter don't pair. Must have found hairpin or multi loop
	if fivePrimeIter > threePrimeIter {
		// hairpin
		en, contribution := hairpinLoopEnergy(fc, pairFivePrimeIdx, pairThreePrimeIdx)
		energy += en
		*energyContributions = append(*energyContributions, contribution)
		return energy
	}

	// (pairFivePrimeIdx, pairThreePrimeIdx) is exterior pair of multi loop
	for fivePrimeIter < pairThreePrimeIdx {
		// add up the contributions of the substructures of the multi loop
		energy += stackEnergy(fc, fivePrimeIter, energyContributions)
		fivePrimeIter = pairTable[fivePrimeIter]
		// search for next base pair in multi loop
		fivePrimeIter++
		for pairTable[fivePrimeIter] == -1 {
			fivePrimeIter++
		}
	}

	en, contribution := multiLoopEnergy(fc, pairFivePrimeIdx, pairTable)
	energy += en
	*energyContributions = append(*energyContributions, contribution)

	return energy
}

/**
 *  Evaluate the free energy contribution of an interior loop with delimiting
 *  base pairs (i,j) and (k,l). See `evaluateInteriorLoop()` for more details.
 */
func interiorLoopEnergy(fc *foldCompound, i, j, k, l int) (int, EnergyContribution) {

	nbUnpairedLeftLoop := k - i - 1
	nbUnpairedRightLoop := j - l - 1

	ijBasePairType := encodedBasePairType(fc.sequence[i], fc.sequence[j],
		fc.energy_params.basePairEncodedTypeMap)
	klBasePairType := encodedBasePairType(fc.sequence[l], fc.sequence[k],
		fc.energy_params.basePairEncodedTypeMap)

	energy := evaluateInteriorLoop(nbUnpairedLeftLoop, nbUnpairedRightLoop,
		ijBasePairType, klBasePairType,
		fc.encodedSequence[i+1], fc.encodedSequence[j-1],
		fc.encodedSequence[k-1], fc.encodedSequence[l+1], fc.energy_params)

	energyContribution := EnergyContribution{
		closingFivePrimeIdx:   i,
		closingThreePrimeIdx:  j,
		enclosedFivePrimeIdx:  k,
		enclosedThreePrimeIdx: l,
		energy:                energy,
		loopType:              InteriorLoop,
	}
	return energy, energyContribution
}

/**
 *  Compute the energy of an interior loop.
 *  This function computes the free energy of an interior loop with the
 *  following structure:
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
func evaluateInteriorLoop(nbUnpairedLeftLoop, nbUnpairedRightLoop int,
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
		// stack
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
				energy += energyParams.TerminalAU
			}

			if enclosedBasePairType > 2 {
				// The encosed base pair is not a GC or CG pair (could be GU, UG, AU, UA, or non-standard)
				energy += energyParams.TerminalAU
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

			// vivek: What does this line mean? Also added above in the 1xn loop
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
	basePairType := encodedBasePairType(fc.sequence[pairFivePrimeIdx],
		fc.sequence[pairThreePrimeIdx], fc.energy_params.basePairEncodedTypeMap)

	energy := evaluateHairpinLoop(nbUnpairedNucleotides, basePairType,
		fc.encodedSequence[pairFivePrimeIdx+1],
		fc.encodedSequence[pairThreePrimeIdx-1],
		fc.sequence[pairFivePrimeIdx-1:], fc.energy_params)

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
		// vivek: this could be a point of failure. Maybe change the above to 7
		// memcpy(tl, sequence, sizeof(char) * 6);
		idx := strings.Index(string(energyParams.Tetraloops), tetraloop)
		if idx != -1 {
			return energyParams.Tetraloop[idx/7]
		}
	} else if size == 6 {
		hexaloop := sequence[:8]
		idx := strings.Index(string(energyParams.Hexaloops), hexaloop)
		if idx != -1 {
			return energyParams.Hexaloop[idx/9]
		}
	} else if size == 3 {
		triloop := sequence[:5]
		idx := strings.Index(string(energyParams.Triloops), triloop)
		if idx != -1 {
			return energyParams.Triloop[idx/6]
		}

		if basePairType > 2 {
			// (X,Y) is not a GC or CG pair (could be GU, UG, AU, UA, or non-standard)
			return energy + energyParams.TerminalAU
		} else {
			// (X,Y) is a GC or CG pair
			return energy
		}
	}

	// size is 5 or greater than 6
	energy += energyParams.mismatchHairpinLoop[basePairType][fivePrimeMismatch][threePrimeMismatch]

	return energy
}

/**
 * pairFivePrimeIdx is the 5'-base of the closing pair
 *
 * since each helix can coaxially stack with at most one of its
 * neighbors we need an auxiliarry variable  cx_energy
 * which contains the best energy given that the last two pairs stack.
 * energy  holds the best energy given the previous two pairs do not
 * stack (i.e. the two current helices may stack)
 * We don't allow the last helix to stack with the first, thus we have to
 * walk around the Loop twice with two starting points and take the minimum
 */
func multiLoopEnergy(fc *foldCompound, closingFivePrimeIdx int, pairTable []int) (int, EnergyContribution) {

	if closingFivePrimeIdx >= pairTable[closingFivePrimeIdx] {
		panic("multiLoopEnergy: pairFivePrimeIdx is not 5' base of a closing pair")
	}

	var closingThreePrimeIdx int
	if closingFivePrimeIdx == 0 {
		closingThreePrimeIdx = fc.length
	} else {
		closingThreePrimeIdx = pairTable[closingFivePrimeIdx]
	}

	// The index of the enclosed base pair at the 5' end
	enclosedFivePrimeIdx := closingFivePrimeIdx + 1

	// seek to the first stem
	for enclosedFivePrimeIdx <= closingThreePrimeIdx && pairTable[enclosedFivePrimeIdx] == -1 {
		enclosedFivePrimeIdx++
	}

	// add inital unpaired nucleotides
	nbUnpairedNucleotides := enclosedFivePrimeIdx - closingFivePrimeIdx - 1

	// MLclosing is energetic penalty imposed when a base pair encloses a
	// multi loop
	energy := fc.energy_params.MLclosing
	for enclosedFivePrimeIdx < closingThreePrimeIdx {
		/* fivePrimeIter must have a pairing partner */
		enclosedThreePrimeIdx := pairTable[enclosedFivePrimeIdx]
		/* get type of the enclosed base pair (enclosedFivePrimeIdx,enclosedThreePrimeIdx) */
		basePairType := encodedBasePairType(fc.sequence[enclosedFivePrimeIdx],
			fc.sequence[enclosedThreePrimeIdx], fc.energy_params.basePairEncodedTypeMap)

		enclosedFivePrimeMismatch := fc.encodedSequence[enclosedFivePrimeIdx-1]
		enclosedThreePrimeMismatch := fc.encodedSequence[enclosedThreePrimeIdx+1]

		energy += multiLoopStemEnergy(basePairType, enclosedFivePrimeMismatch,
			enclosedThreePrimeMismatch, fc.energy_params)

		// seek to the next stem
		enclosedFivePrimeIdx = enclosedThreePrimeIdx + 1

		for enclosedFivePrimeIdx < closingThreePrimeIdx && pairTable[enclosedFivePrimeIdx] == -1 {
			enclosedFivePrimeIdx++
		}

		// add unpaired nucleotides
		nbUnpairedNucleotides += enclosedFivePrimeIdx - enclosedThreePrimeIdx - 1
	}

	if closingFivePrimeIdx > 0 {
		// actual closing pair

		// vivek: Is this a bug? Why are the five and three prime mismatches opposite?
		// Created an issue in the original repo: https://github.com/ViennaRNA/ViennaRNA/issues/126
		basePairType := encodedBasePairType(fc.sequence[closingThreePrimeIdx],
			fc.sequence[closingFivePrimeIdx], fc.energy_params.basePairEncodedTypeMap)
		closingFivePrimeMismatch := fc.encodedSequence[closingThreePrimeIdx-1]
		closingThreePrimeMismatch := fc.encodedSequence[closingFivePrimeIdx+1]

		energy += multiLoopStemEnergy(basePairType, closingFivePrimeMismatch,
			closingThreePrimeMismatch, fc.energy_params)
	} else {
		/* virtual closing pair */
		energy += multiLoopStemEnergy(0, -1, -1, fc.energy_params)
	}

	// add bonus energies for unpaired nucleotides
	energy += nbUnpairedNucleotides * fc.energy_params.MLbase

	energyContribution := EnergyContribution{
		closingFivePrimeIdx:  closingFivePrimeIdx,
		closingThreePrimeIdx: closingThreePrimeIdx,
		energy:               energy,
		loopType:             MultiLoop,
	}
	return energy, energyContribution
}

/**
 *  Compute the energy contribution of a multi-loop stem.
 *
 *  Given a base pair (i,j) (encoded by `basePairType`), compute the
 *  energy contribution including dangling-end/terminal-mismatch contributions.
 *  You can prevent taking 5'-, 3'-dangles or mismatch contributions into
 *  account by passing -1 for `fivePrimeMismatch` and/or `threePrimeMismatch`
 *  respectively.
 *
 *  If either of the adjacent nucleotides i - 1 (`fivePrimeMismatch`) and j + 1
 *  (`threePrimeMismatch`) must not contribute stacking energy, the
 *  corresponding encoding must be -1.
 *
 *  @param  basePairType          The encoded base pair type of (i, j) (see `basePairType()`)
 *  @param  fivePrimeMismatch     The encoded nucleotide directly adjacent (in the 5' direction) to i (may be -1 if index of i is 0)
 *  @param  threePrimeMismatch    The encoded nucleotide directly adjacent (in the 3' direction) to j (may be -1 if index of j is len(sequence) - 1)
 *  @param  energyParams          The pre-computed energy parameters
 *  @return                       The energy contribution of the introduced multi-loop stem
 */
func multiLoopStemEnergy(basePairType, fivePrimeMismatch, threePrimeMismatch int, energyParams *energyParams) int {
	var energy int = 0

	if fivePrimeMismatch >= 0 && threePrimeMismatch >= 0 {
		energy += energyParams.mismatchMultiLoop[basePairType][fivePrimeMismatch][threePrimeMismatch]
	} else if fivePrimeMismatch >= 0 {
		energy += energyParams.dangle5[basePairType][fivePrimeMismatch]
	} else if threePrimeMismatch >= 0 {
		energy += energyParams.dangle3[basePairType][threePrimeMismatch]
	}

	if basePairType > 2 {
		// The closing base pair is not a GC or CG pair (could be GU, UG, AU, UA, or non-standard)
		energy += energyParams.TerminalAU
	}

	energy += energyParams.MLintern[basePairType]

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
		lxc:                    lxc37 * float64(tempf),
		TripleC:                rescaleDg(TripleC37, TripleCdH, tempf),
		MultipleCA:             rescaleDg(MultipleCA37, MultipleCAdH, tempf),
		MultipleCB:             rescaleDg(TerminalAU37, TerminalAUdH, tempf),
		TerminalAU:             rescaleDg(TerminalAU37, TerminalAUdH, tempf),
		MLbase:                 rescaleDg(ML_BASE37, ML_BASEdH, tempf),
		MLclosing:              rescaleDg(ML_closing37, ML_closingdH, tempf),
		basePairEncodedTypeMap: basePairEncodedTypeMap(),
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
		params.Tetraloop[i] = rescaleDg(Tetraloop37[i], TetraloopdH[i], tempf)
	}

	for i := 0; (i * 5) < len(Triloops); i++ {
		params.Triloop[i] = rescaleDg(Triloop37[i], TriloopdH[i], tempf)
	}

	for i := 0; (i * 9) < len(Hexaloops); i++ {
		params.Hexaloop[i] = rescaleDg(Hexaloop37[i], HexaloopdH[i], tempf)
	}

	for i := 0; i <= nbPairs; i++ {
		params.MLintern[i] = rescaleDg(ML_intern37, ML_interndH, tempf)
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

	params.Tetraloops = Tetraloops
	params.Triloops = Triloops
	params.Hexaloops = Hexaloops

	return params
}

func rescaleDg(dG, dH, dT int) int {
	return (dH - (dH-dG)*dT)
}

/**
* Logs energy contributions in the same format as ViennaRNA does. This is used
* to compare output to ViennaRNA in test cases.
* Note: 1 is added to the indexes of the closing and enclosed base pairs as
* everyting in ViennaRNA is 1-indexed, but output from the `CalculateMFE` func
* is 0-indexed.
 */
func logEnergyContributions(energyContribution []EnergyContribution, sequence string) {
	for _, c := range energyContribution {
		switch c.loopType {
		case ExternalLoop:
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
