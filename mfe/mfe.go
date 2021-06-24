package mfe

import (
	"bufio"
	"errors"
	"fmt"
	"log"
	"math"
	"os"
	"regexp"
	"strconv"
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
	// Arbitrary values to denote different types of loops. Used with the
	// `EnergyContribution` struct
	ExteriorLoop = iota
	InteriorLoop
	HairpinLoop
	MultiLoop

	// The number of distinguishable base pairs:
	// CG, GC, GU, UG, AU, UA, & non-standard (see `energy_params.md`)
	nbPairs int = 7
	// The number of distinguishable nucleotides: A, C, G, U
	nbNucleobase int = 4
	// The maximum loop length
	maxLenLoop int = 30
	// DefaultTemperature is the temperature in Celcius for free energy evaluation
	DefaultTemperature float64 = 37.0
	// The temperature at which the energy parameters in `energy_params.go` have
	// been measured at
	energyParamsMeasurementTemperature float64 = 37.0
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
		energyParams:        rawEnergyParamsFromFile("energy_params/rna_turner2004.par").scaleByTemperature(temperature),
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
	var nucleotideRuneEncodedIntMap map[byte]int = map[byte]int{
		'A': 1,
		'C': 2,
		'G': 3,
		'U': 4,
	}

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
Contains all the energy parameters needed for the free energy calculations.

The order of entries to access theses matrices always uses the closing pair or
pairs as the first indices followed by the unpaired bases in 5' to 3' direction.
For example, if we have a 2x2 interior loop:
```
		      5'-GAUA-3'
		      3'-CGCU-5'
```
The closing pairs for the loops are GC and UA (not AU!), and the unpaired bases
are (in 5' to 3' direction, starting at the first pair) A U C G.
Thus, the energy for this sequence is:
```
	pairs:                    GC UA A  U  C  G
						interior2x2Loop[2][6][1][4][2][3]
```
(See `basePairEncodedTypeMap()` and `encodeSequence()` for more details on how
pairs and unpaired nucleotides are encoded)
Note that this sequence is symmetric so the sequence is equivalent to:
```
					5'-UCGC-3'
					3'-AUAG-5'
```
which means the energy of the sequence is equivalent to:
```
	pairs:                    UA GC C  G  A  U
						interior2x2Loop[2][6][1][4][2][3]
```
*/
type energyParams struct {
	/**
	The matrix of free energies for stacked pairs, indexed by the two encoded closing
	pairs. The list should be formatted as symmetric a `7*7` matrix, conforming
	to the order explained above. As an example the stacked pair
	```
						5'-GU-3'
						3'-CA-5'
	```
	corresponds to the entry stackingPair[2][5] (GC=2, AU=5) which should be
	identical to stackingPair[5][2] (AU=5, GC=2).
	size: [nbPairs + 1][nbPairs + 1]int
	*/
	stackingPair [nbPairs + 1][nbPairs + 1]int
	/**
	Free energies of hairpin loops as a function of size. The list should
	contain 31 (maxLenLoop + 1) entries. Since the minimum size of a hairpin loop
	is 3 and we start counting with 0, the first three values should be INF to
	indicate a forbidden value.
	*/
	hairpinLoop [maxLenLoop + 1]int
	/**
	Free energies of bulge loops. Should contain 31 (maxLenLoop + 1) entries, the
	first one being INF.
	*/
	bulge [maxLenLoop + 1]int
	/**
	Free energies of interior loops. Should contain 31 (maxLenLoop + 1) entries,
	the first 4 being INF (since smaller loops are tabulated).

	This field was previous called internal_loop, but has been renamed to
	interior loop to remain consistent with the names of other interior loops
	*/
	interiorLoop [maxLenLoop + 1]int
	/**
		Free energies for the interaction between the closing pair of an interior
	  loop and the two unpaired bases adjacent to the helix. This is a three
	  dimensional array indexed by the type of the closing pair and the two
		unpaired bases. Since we distinguish 5 bases the list contains
	  `7*5*5` entries. The order is such that for example the mismatch

	  ```
	  			       5'-CU-3'
	  			       3'-GC-5'
	  ```
	  corresponds to entry mismatchInteriorLoop[1][4][2] (CG=1, U=4, C=2).

		More information about mismatch energy: https://rna.urmc.rochester.edu/NNDB/turner04/tm.html
	*/
	mismatchInteriorLoop    [nbPairs + 1][nbNucleobase + 1][nbNucleobase + 1]int
	mismatch1xnInteriorLoop [nbPairs + 1][nbNucleobase + 1][nbNucleobase + 1]int
	mismatch2x3InteriorLoop [nbPairs + 1][nbNucleobase + 1][nbNucleobase + 1]int
	mismatchExteriorLoop    [nbPairs + 1][nbNucleobase + 1][nbNucleobase + 1]int
	mismatchHairpinLoop     [nbPairs + 1][nbNucleobase + 1][nbNucleobase + 1]int
	mismatchMultiLoop       [nbPairs + 1][nbNucleobase + 1][nbNucleobase + 1]int
	/**
	Energies for the interaction of an unpaired base on the 5' side and
	adjacent to a helix in multiloops and free ends (the equivalent of mismatch
	energies in interior and hairpin loops). The array is indexed by the type
	of pair closing the helix and the unpaired base and, therefore, forms a `7*5`
	matrix. For example the dangling base in
	```
							5'-GC-3'
							3'- G-5'
	```

	corresponds to entry dangle5[1][3] (CG=1, G=3).

	More information about dangling ends: https://rna.urmc.rochester.edu/NNDB/turner04/de.html
	*/
	dangle5 [nbPairs + 1][nbNucleobase + 1]int
	/**
	Same as above for bases on the 3' side of a helix.
	```
				       5'- A-3'
				       3'-AU-5'
	```
	corresponds to entry dangle3[5][1] (AU=5, A=1).
	*/
	dangle3 [nbPairs + 1][nbNucleobase + 1]int
	/**
	Free energies for symmetric size 2 interior loops. `7*7*5*5` entries.
	Example:
	```
							5'-CUU-3'
							3'-GCA-5'
	```
	corresponds to entry interior1x1Loop[1][5][4][2] (CG=1, AU=5, U=4, C=2),
	which should be identical to interior1x1Loop[5][1][2][4] (AU=5, CG=1, C=2, U=4).
	*/
	interior1x1Loop [nbPairs + 1][nbPairs + 1][nbNucleobase + 1][nbNucleobase + 1]int
	/**
	Free energies for 2x1 interior loops, where 2 is the number of unpaired
	nucelobases on the larger 'side' of the interior loop and 1 is the number of
	unpaired nucelobases on the smaller 'side' of the interior loop.
	`7*7*5*5*5` entries.
	Example:
	```
							5'-CUUU-3'
							3'-GC A-5'
	```
	corresponds to entry interior2x1Loop[1][5][4][4][2] (CG=1, AU=5, U=4, U=4, C=2).
	Note that this matrix is always accessed in the 5' to 3' direction with the
	larger number of unpaired nucleobases first.
	*/
	interior2x1Loop [nbPairs + 1][nbPairs + 1][nbNucleobase + 1][nbNucleobase + 1][nbNucleobase + 1]int
	/**
	Free energies for symmetric size 4 interior loops. `7*7*5*5*5*5` entries.
	Example:
	```
							5'-CUAU-3'
							3'-GCCA-5'
	```
	corresponds to entry interior2x2Loop[1][5][4][1][2][2] (CG=1, AU=5, U=4, A=1, C=2, C=2),
	which should be identical to interior2x2Loop[5][1][2][2][1][4] (AU=5, CG=1, C=2, C=2, A=1, U=4).
	*/
	interior2x2Loop [nbPairs + 1][nbPairs + 1][nbNucleobase + 1][nbNucleobase + 1][nbNucleobase + 1][nbNucleobase + 1]int
	/**
	Parameter for logarithmic loop energy extrapolation. Used to scale energy
	parameters for hairpin, bulge and interior loops when the length of the loop
	is greather than `maxLenLoop`.
	*/
	lxc                                                                                 float64
	multiLoopUnpairedNucelotideBonus, multiLoopClosingPenalty, terminalAUPenalty, ninio int
	multiLoopIntern                                                                     [nbPairs + 1]int
	/**
	Some tetraloops particularly stable tetraloops are assigned an energy
	bonus. For example:
	```
		GAAA    -200
	```
	assigns a bonus energy of -2 kcal/mol to tetraLoop containing
	the sequence GAAA.
	*/
	tetraLoop map[string]int
	triLoop   map[string]int
	hexaLoop  map[string]int
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

func evaluateFoldCompound(fc *foldCompound) (int, []EnergyContribution) {
	energyContributions := make([]EnergyContribution, 0)

	// get the energy contributions due to exterior loops
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
		// seek to end of current loop
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
* See `exteriorStemEnergy()` for information of how exterior loops are evaluated.
* More information available at: https://rna.urmc.rochester.edu/NNDB/turner04/exterior.html
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
		energy += energyParams.terminalAUPenalty
	}

	return energy
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
func stackEnergy(fc *foldCompound, closingFivePrimeIdx int) (int, []EnergyContribution) {
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
 *  base pairs (i,j) and (k,l). See `evaluateStackBulgeInteriorLoop()` for more details.
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

		if nbUnpairedLarger <= maxLenLoop {
			energy = energyParams.bulge[nbUnpairedLarger]
		} else {
			energy = energyParams.bulge[maxLenLoop] + int(energyParams.lxc*math.Log(float64(nbUnpairedLarger)/float64(maxLenLoop)))
		}

		if nbUnpairedLarger == 1 {
			energy += energyParams.stackingPair[closingBasePairType][enclosedBasePairType]
		} else {
			if closingBasePairType > 2 {
				// The closing base pair is not a GC or CG pair (could be GU, UG, AU, UA, or non-standard)
				energy += energyParams.terminalAUPenalty
			}

			if enclosedBasePairType > 2 {
				// The encosed base pair is not a GC or CG pair (could be GU, UG, AU, UA, or non-standard)
				energy += energyParams.terminalAUPenalty
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
				if nbUnpairedLarger+1 <= maxLenLoop {
					energy = energyParams.interiorLoop[nbUnpairedLarger+1]
				} else {
					energy = energyParams.interiorLoop[maxLenLoop] + int(energyParams.lxc*math.Log((float64(nbUnpairedLarger)+1.0)/float64(maxLenLoop)))
				}
				energy += min(maxNinio, (nbUnpairedLarger-nbUnpairedSmaller)*energyParams.ninio)
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
				energy = energyParams.interiorLoop[nbNucleobase+1] + energyParams.ninio
				energy += energyParams.mismatch2x3InteriorLoop[closingBasePairType][closingFivePrimeMismatch][closingThreePrimeMismatch] + energyParams.mismatch2x3InteriorLoop[enclosedBasePairType][enclosedFivePrimeMismatch][enclosedThreePrimeMismatch]
				return energy
			}
		}

		{
			/* generic interior loop (no else here!)*/
			nbUnpairedNucleotides := nbUnpairedLarger + nbUnpairedSmaller

			var energy int
			if nbUnpairedNucleotides <= maxLenLoop {
				energy = energyParams.interiorLoop[nbUnpairedNucleotides]
			} else {
				energy = energyParams.interiorLoop[maxLenLoop] + int(energyParams.lxc*math.Log(float64(nbUnpairedNucleotides)/float64(maxLenLoop)))
			}

			energy += min(maxNinio, (nbUnpairedLarger-nbUnpairedSmaller)*energyParams.ninio)

			energy += energyParams.mismatchInteriorLoop[closingBasePairType][closingFivePrimeMismatch][closingThreePrimeMismatch] + energyParams.mismatchInteriorLoop[enclosedBasePairType][enclosedFivePrimeMismatch][enclosedThreePrimeMismatch]
			return energy
		}
	}
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
func evaluateHairpinLoop(size, basePairType, fivePrimeMismatch, threePrimeMismatch int, sequence string, energyParams *energyParams) int {
	var energy int

	if size <= maxLenLoop {
		energy = energyParams.hairpinLoop[size]
	} else {
		energy = energyParams.hairpinLoop[maxLenLoop] + int(energyParams.lxc*math.Log(float64(size)/float64(maxLenLoop)))
	}

	if size < 3 {
		// should only be the case when folding alignments
		return energy
	}

	if size == 4 {
		tetraLoop := sequence[:6]
		tetraLoopEnergy, present := energyParams.tetraLoop[tetraLoop]
		if present {
			return tetraLoopEnergy
		}
	} else if size == 6 {
		hexaLoop := sequence[:8]
		hexaLoopEnergy, present := energyParams.hexaLoop[hexaLoop]
		if present {
			return hexaLoopEnergy
		}
	} else if size == 3 {
		triLoop := sequence[:5]
		triLoopEnergy, present := energyParams.triLoop[triLoop]
		if present {
			return triLoopEnergy
		}

		if basePairType > 2 {
			// (X,Y) is not a GC or CG pair (could be GU, UG, AU, UA, or non-standard)
			return energy + energyParams.terminalAUPenalty
		}

		// (X,Y) is a GC or CG pair
		return energy
	}

	energy += energyParams.mismatchHairpinLoop[basePairType][fivePrimeMismatch][threePrimeMismatch]

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
 * multiLoopEnergy iterates through the multi-loop and finds each base pair
 * enclosed in the multi-loop. For each enclosed base pair, we add to the total
 * the energy due to that base pair and due to the substructure closed by the
 * enclosed base pair.
 * We also add a bonus energy based on the number of unpaired nucleotides in
 * the multi-loop (though it defaults to 0 right now).
 *
 * More information about multi-loops: https://rna.urmc.rochester.edu/NNDB/turner04/mb.html
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
	multiLoopEnergy += nbUnpairedNucleotides * fc.energyParams.multiLoopUnpairedNucelotideBonus

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
		energy += energyParams.terminalAUPenalty
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

func turner2004Params() rawEnergyParams {
	return rawEnergyParams{
		interior2x2LoopEnergy37C: interior2x2LoopEnergy37,
		interior2x2LoopEnthalpy:  interior2x2LoopEnthalpy,

		interior2x1LoopEnergy37C: interior2x1LoopEnergy37,
		interior2x1LoopEnthalpy:  interior2x1LoopEnthalpy,

		interior1x1LoopEnergy37C: interior1x1LoopEnergy37C,
		interior1x1LoopEnthalpy:  interior1x1LoopEnthalpy,

		mismatchInteriorLoopEnergy37C:    mismatchInteriorLoopEnergy37C,
		mismatchInteriorLoopEnthalpy:     mismatchInteriorLoopEnthalpy,
		mismatchHairpinLoopEnergy37C:     mismatchHairpinLoopEnergy37C,
		mismatchHairpinLoopEnthalpy:      mismatchHairpinLoopEnthalpy,
		mismatchMultiLoopEnergy37C:       mismatchMultiLoopEnergy37C,
		mismatchMultiLoopEnthalpy:        mismatchMultiLoopEnthalpy,
		mismatch1xnInteriorLoopEnergy37C: mismatch1xnInteriorLoopEnergy37C,
		mismatch1xnInteriorLoopEnthalpy:  mismatch1xnInteriorLoopEnthalpy,
		mismatch2x3InteriorLoopEnergy37C: mismatch2x3InteriorLoopEnergy37C,
		mismatch2x3InteriorLoopEnthalpy:  mismatch2x3InteriorLoopEnthalpy,
		mismatchExteriorLoopEnergy37C:    mismatchExteriorLoopEnergy37C,
		mismatchExteriorLoopEnthalpy:     mismatchExteriorLoopEnthalpy,

		dangle5Energy37C:      dangle5Energy37C,
		dangle5Enthalpy:       dangle5Enthalpy,
		dangle3Energy37C:      dangle3Energy37C,
		dangle3Enthalpy:       dangle3Enthalpy,
		stackingPairEnergy37C: stackingPairEnergy37C,
		stackingPairEnthalpy:  stackingPairEnthalpy,

		hairpinLoopEnergy37C:  hairpinLoopEnergy37C[:],
		hairpinLoopEnthalpy:   hairpinLoopEnthalpy[:],
		bulgeEnergy37C:        bulgeEnergy37C[:],
		bulgeEnthalpy:         bulgeEnthalpy[:],
		interiorLoopEnergy37C: interiorLoopEnergy37C[:],
		interiorLoopEnthalpy:  interiorLoopEnthalpy[:],

		multiLoopIntern37C:       multiLoopIntern37C,
		multiLoopInternEnthalpy:  multiLoopInternEnthalpy,
		multiLoopClosing37C:      multiLoopClosing37C,
		multiLoopClosingEnthalpy: multiLoopClosingEnthalpy,
		multiLoopBase37C:         multiLoopBase37C,
		multiLoopBaseEnthalpy:    multiLoopBaseEnthalpy,
		maxNinio:                 maxNinio,
		ninio37C:                 ninio37C,
		ninioEnthalpy:            ninioEnthalpy,
		terminalAU37C:            terminalAU37C,
		terminalAUEnthalpy:       terminalAUEnthalpy,

		tetraLoops:         tetraLoops,
		tetraLoopEnergy37C: tetraLoopEnergy37C,
		tetraLoopEnthalpy:  tetraLoopEnthalpy,

		triLoops:         triLoops,
		triLoopEnergy37C: triLoopEnergy37C,
		triLoopEnthalpy:  triLoopEnthalpy,

		hexaLoops:         hexaLoops,
		hexaLoopEnergy37C: hexaLoopEnergy37C,
		hexaLoopEnthalpy:  hexaLoopEnthalpy,
	}
}

// scales energy paramaters according to the specificed temperatue
func (rawEnergyParams rawEnergyParams) scaleByTemperature(temperature float64) *energyParams {

	// set the non-matix energy parameters
	var params *energyParams = &energyParams{
		lxc:                              rescaleDgFloat64(lxc37, 0, temperature),
		terminalAUPenalty:                rescaleDg(rawEnergyParams.terminalAU37C, rawEnergyParams.terminalAUEnthalpy, temperature),
		multiLoopUnpairedNucelotideBonus: rescaleDg(rawEnergyParams.multiLoopBase37C, rawEnergyParams.multiLoopBaseEnthalpy, temperature),
		multiLoopClosingPenalty:          rescaleDg(rawEnergyParams.multiLoopClosing37C, rawEnergyParams.multiLoopClosingEnthalpy, temperature),
		ninio:                            rescaleDg(rawEnergyParams.ninio37C, rawEnergyParams.ninioEnthalpy, temperature),
	}

	// scale and set hairpin, bulge and interior energy params
	for i := 0; i <= maxLenLoop; i++ {
		params.hairpinLoop[i] = rescaleDg(rawEnergyParams.hairpinLoopEnergy37C[i], rawEnergyParams.hairpinLoopEnthalpy[i], temperature)
		params.bulge[i] = rescaleDg(rawEnergyParams.bulgeEnergy37C[i], rawEnergyParams.bulgeEnthalpy[i], temperature)
		params.interiorLoop[i] = rescaleDg(rawEnergyParams.interiorLoopEnergy37C[i], rawEnergyParams.interiorLoopEnthalpy[i], temperature)
	}

	params.tetraLoop = make(map[string]int)
	for k := range tetraLoops {
		params.tetraLoop[k] = rescaleDg(rawEnergyParams.tetraLoopEnergy37C[k], rawEnergyParams.tetraLoopEnthalpy[k], temperature)
	}

	params.triLoop = make(map[string]int)
	for k := range triLoops {
		params.triLoop[k] = rescaleDg(rawEnergyParams.triLoopEnergy37C[k], rawEnergyParams.triLoopEnthalpy[k], temperature)
	}

	params.hexaLoop = make(map[string]int)
	for k := range hexaLoops {
		params.hexaLoop[k] = rescaleDg(rawEnergyParams.hexaLoopEnergy37C[k], rawEnergyParams.hexaLoopEnthalpy[k], temperature)
	}

	for i := 0; i <= nbPairs; i++ {
		params.multiLoopIntern[i] = rescaleDg(rawEnergyParams.multiLoopIntern37C, rawEnergyParams.multiLoopInternEnthalpy, temperature)
	}

	/* stacks    G(T) = H - [H - G(T0)]*T/T0 */
	for i := 0; i <= nbPairs; i++ {
		for j := 0; j <= nbPairs; j++ {
			params.stackingPair[i][j] = rescaleDg(rawEnergyParams.stackingPairEnergy37C[i][j],
				rawEnergyParams.stackingPairEnthalpy[i][j],
				temperature)
		}
	}

	/* mismatches */
	for i := 0; i <= nbPairs; i++ {
		for j := 0; j <= nbNucleobase; j++ {
			for k := 0; k <= nbNucleobase; k++ {
				var mismatch int
				params.mismatchInteriorLoop[i][j][k] = rescaleDg(rawEnergyParams.mismatchInteriorLoopEnergy37C[i][j][k],
					rawEnergyParams.mismatchInteriorLoopEnthalpy[i][j][k],
					temperature)
				params.mismatchHairpinLoop[i][j][k] = rescaleDg(rawEnergyParams.mismatchHairpinLoopEnergy37C[i][j][k],
					rawEnergyParams.mismatchHairpinLoopEnthalpy[i][j][k],
					temperature)
				params.mismatch1xnInteriorLoop[i][j][k] = rescaleDg(rawEnergyParams.mismatch1xnInteriorLoopEnergy37C[i][j][k],
					rawEnergyParams.mismatch1xnInteriorLoopEnthalpy[i][j][k],
					temperature)
				params.mismatch2x3InteriorLoop[i][j][k] = rescaleDg(rawEnergyParams.mismatch2x3InteriorLoopEnergy37C[i][j][k],
					rawEnergyParams.mismatch2x3InteriorLoopEnthalpy[i][j][k],
					temperature)

				mismatch = rescaleDg(rawEnergyParams.mismatchMultiLoopEnergy37C[i][j][k],
					rawEnergyParams.mismatchMultiLoopEnthalpy[i][j][k],
					temperature)
				params.mismatchMultiLoop[i][j][k] = min(0, mismatch)

				mismatch = rescaleDg(rawEnergyParams.mismatchExteriorLoopEnergy37C[i][j][k],
					rawEnergyParams.mismatchExteriorLoopEnthalpy[i][j][k],
					temperature)
				params.mismatchExteriorLoop[i][j][k] = min(0, mismatch)
			}
		}
	}

	/* dangles */
	for i := 0; i <= nbPairs; i++ {
		for j := 0; j <= nbNucleobase; j++ {
			var dd int
			dd = rescaleDg(rawEnergyParams.dangle5Energy37C[i][j],
				rawEnergyParams.dangle5Enthalpy[i][j],
				temperature)
			params.dangle5[i][j] = min(0, dd)

			dd = rescaleDg(rawEnergyParams.dangle3Energy37C[i][j],
				rawEnergyParams.dangle3Enthalpy[i][j],
				temperature)
			params.dangle3[i][j] = min(0, dd)
		}
	}

	/* interior 1x1 loops */
	for i := 0; i <= nbPairs; i++ {
		for j := 0; j <= nbPairs; j++ {
			for k := 0; k <= nbNucleobase; k++ {
				for l := 0; l <= nbNucleobase; l++ {
					params.interior1x1Loop[i][j][k][l] = rescaleDg(rawEnergyParams.interior1x1LoopEnergy37C[i][j][k][l],
						rawEnergyParams.interior1x1LoopEnthalpy[i][j][k][l],
						temperature)
				}
			}
		}
	}

	/* interior 2x1 loops */
	for i := 0; i < nbPairs; i++ {
		for j := 0; j <= nbPairs; j++ {
			for k := 0; k <= nbNucleobase; k++ {
				for l := 0; l <= nbNucleobase; l++ {
					for m := 0; m <= nbNucleobase; m++ {
						params.interior2x1Loop[i][j][k][l][m] = rescaleDg(rawEnergyParams.interior2x1LoopEnergy37C[i][j][k][l][m],
							rawEnergyParams.interior2x1LoopEnthalpy[i][j][k][l][m],
							temperature)
					}
				}
			}
		}
	}

	/* interior 2x2 loops */
	for i := 0; i < nbPairs; i++ {
		for j := 0; j <= nbPairs; j++ {
			for k := 0; k <= nbNucleobase; k++ {
				for l := 0; l <= nbNucleobase; l++ {
					for m := 0; m <= nbNucleobase; m++ {
						for n := 0; n <= nbNucleobase; n++ {
							params.interior2x2Loop[i][j][k][l][m][n] = rescaleDg(rawEnergyParams.interior2x2LoopEnergy37C[i][j][k][l][m][n],
								rawEnergyParams.interior2x2LoopEnthalpy[i][j][k][l][m][n],
								temperature)
						}
					}
				}
			}
		}
	}

	return params
}

func (rawEnergyParams *rawEnergyParams) setMultiLoopParams(parsedMultiLoopParams []int) {
	rawEnergyParams.multiLoopBase37C = parsedMultiLoopParams[0]
	rawEnergyParams.multiLoopBaseEnthalpy = parsedMultiLoopParams[1]
	rawEnergyParams.multiLoopClosing37C = parsedMultiLoopParams[2]
	rawEnergyParams.multiLoopClosingEnthalpy = parsedMultiLoopParams[3]
	rawEnergyParams.multiLoopIntern37C = parsedMultiLoopParams[4]
	rawEnergyParams.multiLoopInternEnthalpy = parsedMultiLoopParams[5]
}

func (rawEnergyParams *rawEnergyParams) setNinioParams(parsedNinioParams []int) {
	rawEnergyParams.ninio37C = parsedNinioParams[0]
	rawEnergyParams.ninioEnthalpy = parsedNinioParams[1]
	rawEnergyParams.maxNinio = parsedNinioParams[2]
}

/**
Rescale Gibbs free energy according to the equation dG = dH - T * dS
where dG is the change in Gibbs free energy
			dH is the change in enthalpy
			dS is the change in entrophy
			T is the temperature
*/
func rescaleDg(dG, dH int, temperature float64) int {
	// if temperate == energyParamsTemperature then below calculation will always
	// return dG. So we save some computation with this check.
	if temperature == energyParamsMeasurementTemperature {
		return dG
	}

	defaultEnergyParamsTemperatureKelvin := energyParamsMeasurementTemperature + zeroCKelvin
	temperatureKelvin := temperature + zeroCKelvin
	var T float64 = float64(temperatureKelvin / defaultEnergyParamsTemperatureKelvin)

	dGFloat64 := float64(dG)
	dHFloat64 := float64(dH)

	dSFloat64 := dHFloat64 - dGFloat64

	return int(dHFloat64 - dSFloat64*T)
}

// same as rescaleDg, but for floats
func rescaleDgFloat64(dG, dH, temperature float64) float64 {
	// if temperate == energyParamsTemperature then below calculation will always
	// return dG. So we save some computation with this check.
	if temperature == energyParamsMeasurementTemperature {
		return dG
	}

	defaultEnergyParamsTemperatureKelvin := energyParamsMeasurementTemperature + zeroCKelvin
	temperatureKelvin := temperature + zeroCKelvin
	var T float64 = temperatureKelvin / defaultEnergyParamsTemperatureKelvin

	dS := dH - dG
	return dH - dS*T
}

/**
PrintEnergyContributions pretty prints a list of energy contributions returned
from `MinimumFreeEnergy`.
*/
func PrintEnergyContributions(energyContribution []EnergyContribution, sequence string) {
	for _, c := range energyContribution {
		switch c.loopType {
		case ExteriorLoop:
			fmt.Printf("External loop                           : %v\n", c.energy)
		case InteriorLoop:
			var i, j int = c.closingFivePrimeIdx, c.closingThreePrimeIdx
			var k, l int = c.enclosedFivePrimeIdx, c.enclosedThreePrimeIdx
			fmt.Printf("Interior loop (%v,%v) %v%v; (%v,%v) %v%v: %v\n",
				i, j,
				string(sequence[i]), string(sequence[j]),
				k, l,
				string(sequence[k]), string(sequence[l]),
				c.energy)
		case HairpinLoop:
			var i, j int = c.closingFivePrimeIdx, c.closingThreePrimeIdx
			fmt.Printf("Hairpin loop  (%v,%v) %v%v              : %v\n",
				i, j,
				string(sequence[i]), string(sequence[j]),
				c.energy)
		case MultiLoop:
			var i, j int = c.closingFivePrimeIdx, c.closingThreePrimeIdx
			fmt.Printf("Multi   loop (%v,%v) %v%v              : %v\n",
				i, j,
				string(sequence[i]), string(sequence[j]),
				c.energy)
		}
	}
}

type rawEnergyParams struct {
	interior2x2LoopEnergy37C [][][][][][]int
	interior2x2LoopEnthalpy  [][][][][][]int

	interior2x1LoopEnergy37C [][][][][]int
	interior2x1LoopEnthalpy  [][][][][]int

	interior1x1LoopEnergy37C [][][][]int
	interior1x1LoopEnthalpy  [][][][]int

	mismatchInteriorLoopEnergy37C    [][][]int
	mismatchInteriorLoopEnthalpy     [][][]int
	mismatchHairpinLoopEnergy37C     [][][]int
	mismatchHairpinLoopEnthalpy      [][][]int
	mismatchMultiLoopEnergy37C       [][][]int
	mismatchMultiLoopEnthalpy        [][][]int
	mismatch1xnInteriorLoopEnergy37C [][][]int
	mismatch1xnInteriorLoopEnthalpy  [][][]int
	mismatch2x3InteriorLoopEnergy37C [][][]int
	mismatch2x3InteriorLoopEnthalpy  [][][]int
	mismatchExteriorLoopEnergy37C    [][][]int
	mismatchExteriorLoopEnthalpy     [][][]int

	dangle5Energy37C      [][]int
	dangle5Enthalpy       [][]int
	dangle3Energy37C      [][]int
	dangle3Enthalpy       [][]int
	stackingPairEnergy37C [][]int
	stackingPairEnthalpy  [][]int

	hairpinLoopEnergy37C  []int
	hairpinLoopEnthalpy   []int
	bulgeEnergy37C        []int
	bulgeEnthalpy         []int
	interiorLoopEnergy37C []int
	interiorLoopEnthalpy  []int

	multiLoopIntern37C       int
	multiLoopInternEnthalpy  int
	multiLoopClosing37C      int
	multiLoopClosingEnthalpy int
	multiLoopBase37C         int
	multiLoopBaseEnthalpy    int
	maxNinio                 int
	ninio37C                 int
	ninioEnthalpy            int
	terminalAU37C            int
	terminalAUEnthalpy       int

	tetraLoops         map[string]bool
	tetraLoopEnergy37C map[string]int
	tetraLoopEnthalpy  map[string]int

	triLoops         map[string]bool
	triLoopEnergy37C map[string]int
	triLoopEnthalpy  map[string]int

	hexaLoops         map[string]bool
	hexaLoopEnergy37C map[string]int
	hexaLoopEnthalpy  map[string]int
}

func readLine(scanner *bufio.Scanner) (string, bool) {
	// scanner.Split(bufio.ScanLines)
	lineAvailable := scanner.Scan()

	if err := scanner.Err(); err != nil {
		log.Fatal(err)
	}

	if lineAvailable {
		line := scanner.Text()
		return line, lineAvailable
	} else {
		return "", false
	}
}

func rawEnergyParamsFromFile(fileName string) (rawEnergyParams rawEnergyParams) {
	file, err := os.Open(fileName)
	if err != nil {
		log.Fatal(err)
	}
	defer file.Close()

	scanner := bufio.NewScanner(file)

	headerLine, lineAvailable := readLine(scanner)
	if lineAvailable && headerLine != "## RNAfold parameter file v2.0" {
		panic("Missing header line in file.\nMay be this file has not v2.0 format.\n")
	}

	// size_t      line_no;
	// char        *line, ident[256];
	// enum parset type;
	// int         r;

	// line_no = 0

	// if ((!file_content) ||
	//     (!file_content[line_no])) {
	//   return 0
	// }

	// if (strncmp(file_content[line_no++], "## RNAfold parameter file v2.0", 30) != 0) {
	//   vrna_message_warning("Missing header line in file.\n"
	//                        "May be this file has not v2.0 format.\n"
	//                        "Use INTERRUPT-key to stop.");
	// }
	line, lineAvailable := readLine(scanner)
	for lineAvailable {
		var energyParamSet string
		n, _ := fmt.Sscanf(line, "# %255s", &energyParamSet)

		if n == 1 {
			switch energyParamSet {
			case "END":
				break
			case "stack":
				rawEnergyParams.stackingPairEnergy37C = parseParamValues(scanner, nbPairs, nbPairs).([][]int)
				rawEnergyParams.stackingPairEnergy37C = addPreOffset2Dim(rawEnergyParams.stackingPairEnergy37C, 1, 1)
				// fmt.Println(testEq(rawEnergyParams.stackingPairEnergy37C, stackingPairEnergy37C, 2))
			case "stack_enthalpies":
				rawEnergyParams.stackingPairEnthalpy = parseParamValues(scanner, nbPairs, nbPairs).([][]int)
				rawEnergyParams.stackingPairEnthalpy = addPreOffset(rawEnergyParams.stackingPairEnthalpy, 1, 1).([][]int)
				// fmt.Println(testEq(rawEnergyParams.stackingPairEnthalpy, stackingPairEnthalpy, 2))
			case "hairpin":
				rawEnergyParams.hairpinLoopEnergy37C = parseParamValues(scanner, 31).([]int)
				// fmt.Println(testEq(rawEnergyParams.hairpinLoopEnergy37C, hairpinLoopEnergy37C[:], 1))
			case "hairpin_enthalpies":
				rawEnergyParams.hairpinLoopEnthalpy = parseParamValues(scanner, 31).([]int)
				// fmt.Println(testEq(rawEnergyParams.hairpinLoopEnthalpy, hairpinLoopEnthalpy[:], 1))
			case "bulge":
				rawEnergyParams.bulgeEnergy37C = parseParamValues(scanner, 31).([]int)
				// fmt.Println(testEq(rawEnergyParams.bulgeEnergy37C, bulgeEnergy37C[:], 1))
			case "bulge_enthalpies":
				rawEnergyParams.bulgeEnthalpy = parseParamValues(scanner, 31).([]int)
				// fmt.Println(testEq(rawEnergyParams.bulgeEnthalpy, bulgeEnthalpy[:], 1))
			case "interior":
				rawEnergyParams.interiorLoopEnergy37C = parseParamValues(scanner, 31).([]int)
				// fmt.Println(testEq(rawEnergyParams.interiorLoopEnergy37C, interiorLoopEnergy37C[:], 1))
			case "interior_enthalpies":
				rawEnergyParams.interiorLoopEnthalpy = parseParamValues(scanner, 31).([]int)
				// fmt.Println(testEq(rawEnergyParams.interiorLoopEnthalpy, interiorLoopEnthalpy[:], 1))
			case "mismatch_exterior":
				rawEnergyParams.mismatchExteriorLoopEnergy37C = parseParamValues(scanner, nbPairs, nbNucleobase+1, nbNucleobase+1).([][][]int)
				rawEnergyParams.mismatchExteriorLoopEnergy37C = addPreOffset(rawEnergyParams.mismatchExteriorLoopEnergy37C, 1, 0, 0).([][][]int)
				// fmt.Println(testEq(rawEnergyParams.mismatchExteriorLoopEnergy37C, mismatchExteriorLoopEnergy37C, 3))
			case "mismatch_exterior_enthalpies":
				rawEnergyParams.mismatchExteriorLoopEnthalpy = parseParamValues(scanner, nbPairs, nbNucleobase+1, nbNucleobase+1).([][][]int)
				rawEnergyParams.mismatchExteriorLoopEnthalpy = addPreOffset(rawEnergyParams.mismatchExteriorLoopEnthalpy, 1, 0, 0).([][][]int)
				// fmt.Println(testEq(rawEnergyParams.mismatchExteriorLoopEnthalpy, mismatchExteriorLoopEnthalpy, 3))
			case "mismatch_hairpin":
				rawEnergyParams.mismatchHairpinLoopEnergy37C = parseParamValues(scanner, nbPairs, nbNucleobase+1, nbNucleobase+1).([][][]int)
				rawEnergyParams.mismatchHairpinLoopEnergy37C = addPreOffset(rawEnergyParams.mismatchHairpinLoopEnergy37C, 1, 0, 0).([][][]int)
				// fmt.Println(testEq(rawEnergyParams.mismatchHairpinLoopEnergy37C, mismatchHairpinLoopEnergy37C, 3))
			case "mismatch_hairpin_enthalpies":
				rawEnergyParams.mismatchHairpinLoopEnthalpy = parseParamValues(scanner, nbPairs, nbNucleobase+1, nbNucleobase+1).([][][]int)
				rawEnergyParams.mismatchHairpinLoopEnthalpy = addPreOffset(rawEnergyParams.mismatchHairpinLoopEnthalpy, 1, 0, 0).([][][]int)
				// fmt.Println(testEq(rawEnergyParams.mismatchHairpinLoopEnthalpy, mismatchHairpinLoopEnthalpy, 3))
			case "mismatch_interior":
				rawEnergyParams.mismatchInteriorLoopEnergy37C = parseParamValues(scanner, nbPairs, nbNucleobase+1, nbNucleobase+1).([][][]int)
				rawEnergyParams.mismatchInteriorLoopEnergy37C = addPreOffset(rawEnergyParams.mismatchInteriorLoopEnergy37C, 1, 0, 0).([][][]int)
				// fmt.Println(testEq(rawEnergyParams.mismatchInteriorLoopEnergy37C, mismatchInteriorLoopEnergy37C, 3))
			case "mismatch_interior_enthalpies":
				rawEnergyParams.mismatchInteriorLoopEnthalpy = parseParamValues(scanner, nbPairs, nbNucleobase+1, nbNucleobase+1).([][][]int)
				rawEnergyParams.mismatchInteriorLoopEnthalpy = addPreOffset(rawEnergyParams.mismatchInteriorLoopEnthalpy, 1, 0, 0).([][][]int)
				// fmt.Println(testEq(rawEnergyParams.mismatchInteriorLoopEnthalpy, mismatchInteriorLoopEnthalpy, 3))
			case "mismatch_interior_1n":
				rawEnergyParams.mismatch1xnInteriorLoopEnergy37C = parseParamValues(scanner, nbPairs, nbNucleobase+1, nbNucleobase+1).([][][]int)
				rawEnergyParams.mismatch1xnInteriorLoopEnergy37C = addPreOffset(rawEnergyParams.mismatch1xnInteriorLoopEnergy37C, 1, 0, 0).([][][]int)
				// fmt.Println(testEq(rawEnergyParams.mismatch1xnInteriorLoopEnergy37C, mismatch1xnInteriorLoopEnergy37C, 3))
			case "mismatch_interior_1n_enthalpies":
				rawEnergyParams.mismatch1xnInteriorLoopEnthalpy = parseParamValues(scanner, nbPairs, nbNucleobase+1, nbNucleobase+1).([][][]int)
				rawEnergyParams.mismatch1xnInteriorLoopEnthalpy = addPreOffset(rawEnergyParams.mismatch1xnInteriorLoopEnthalpy, 1, 0, 0).([][][]int)
				// fmt.Println(testEq(rawEnergyParams.mismatch1xnInteriorLoopEnthalpy, mismatch1xnInteriorLoopEnthalpy, 3))
			case "mismatch_interior_23":
				rawEnergyParams.mismatch2x3InteriorLoopEnergy37C = parseParamValues(scanner, nbPairs, nbNucleobase+1, nbNucleobase+1).([][][]int)
				rawEnergyParams.mismatch2x3InteriorLoopEnergy37C = addPreOffset(rawEnergyParams.mismatch2x3InteriorLoopEnergy37C, 1, 0, 0).([][][]int)
				// fmt.Println(testEq(rawEnergyParams.mismatch2x3InteriorLoopEnergy37C, mismatch2x3InteriorLoopEnergy37C, 3))
			case "mismatch_interior_23_enthalpies":
				rawEnergyParams.mismatch2x3InteriorLoopEnthalpy = parseParamValues(scanner, nbPairs, nbNucleobase+1, nbNucleobase+1).([][][]int)
				rawEnergyParams.mismatch2x3InteriorLoopEnthalpy = addPreOffset(rawEnergyParams.mismatch2x3InteriorLoopEnthalpy, 1, 0, 0).([][][]int)
				// fmt.Println(testEq(rawEnergyParams.mismatch2x3InteriorLoopEnthalpy, mismatch2x3InteriorLoopEnthalpy, 3))

			case "mismatch_multi":
				rawEnergyParams.mismatchMultiLoopEnergy37C = parseParamValues(scanner, nbPairs, nbNucleobase+1, nbNucleobase+1).([][][]int)
				rawEnergyParams.mismatchMultiLoopEnergy37C = addPreOffset(rawEnergyParams.mismatchMultiLoopEnergy37C, 1, 0, 0).([][][]int)
				// fmt.Println(testEq(rawEnergyParams.mismatchMultiLoopEnergy37C, mismatchMultiLoopEnergy37C, 3))
			case "mismatch_multi_enthalpies":
				rawEnergyParams.mismatchMultiLoopEnthalpy = parseParamValues(scanner, nbPairs, nbNucleobase+1, nbNucleobase+1).([][][]int)
				rawEnergyParams.mismatchMultiLoopEnthalpy = addPreOffset(rawEnergyParams.mismatchMultiLoopEnthalpy, 1, 0, 0).([][][]int)
				// fmt.Println(testEq(rawEnergyParams.mismatchMultiLoopEnthalpy, mismatchMultiLoopEnthalpy, 3))
			case "int11":
				rawEnergyParams.interior1x1LoopEnergy37C = parseParamValues(scanner, nbPairs, nbPairs, nbNucleobase+1, nbNucleobase+1).([][][][]int)
				rawEnergyParams.interior1x1LoopEnergy37C = addPreOffset(rawEnergyParams.interior1x1LoopEnergy37C, 1, 1, 0, 0).([][][][]int)
				// fmt.Println(testEq(rawEnergyParams.interior1x1LoopEnergy37C, interior1x1LoopEnergy37C, 4))
			case "int11_enthalpies":
				rawEnergyParams.interior1x1LoopEnthalpy = parseParamValues(scanner, nbPairs, nbPairs, nbNucleobase+1, nbNucleobase+1).([][][][]int)
				rawEnergyParams.interior1x1LoopEnthalpy = addPreOffset(rawEnergyParams.interior1x1LoopEnthalpy, 1, 1, 0, 0).([][][][]int)
				// fmt.Println(testEq(rawEnergyParams.interior1x1LoopEnthalpy, interior1x1LoopEnthalpy, 4))

			case "int21":
				rawEnergyParams.interior2x1LoopEnergy37C = parseParamValues(scanner, nbPairs, nbPairs, nbNucleobase+1, nbNucleobase+1, nbNucleobase+1).([][][][][]int)
				rawEnergyParams.interior2x1LoopEnergy37C = addPreOffset(rawEnergyParams.interior2x1LoopEnergy37C, 1, 1, 0, 0, 0).([][][][][]int)
				// fmt.Println(testEq(rawEnergyParams.interior2x1LoopEnergy37C, interior2x1LoopEnergy37, 5))
			case "int21_enthalpies":
				rawEnergyParams.interior2x1LoopEnthalpy = parseParamValues(scanner, nbPairs, nbPairs, nbNucleobase+1, nbNucleobase+1, nbNucleobase+1).([][][][][]int)
				rawEnergyParams.interior2x1LoopEnthalpy = addPreOffset(rawEnergyParams.interior2x1LoopEnthalpy, 1, 1, 0, 0, 0).([][][][][]int)
				// fmt.Println(testEq(rawEnergyParams.interior2x1LoopEnthalpy, interior2x1LoopEnthalpy, 5))
			case "int22":
				rawEnergyParams.interior2x2LoopEnergy37C = parseParamValues(scanner, nbPairs-1, nbPairs-1, nbNucleobase, nbNucleobase, nbNucleobase, nbNucleobase).([][][][][][]int)
				rawEnergyParams.interior2x2LoopEnergy37C = addPreOffset(rawEnergyParams.interior2x2LoopEnergy37C, 1, 1, 1, 1, 1, 1).([][][][][][]int)
				rawEnergyParams.interior2x2LoopEnergy37C = addPostOffset(rawEnergyParams.interior2x2LoopEnergy37C, 1, 1, 0, 0, 0, 0).([][][][][][]int)
				rawEnergyParams.interior2x2LoopEnergy37C = updateNST(rawEnergyParams.interior2x2LoopEnergy37C)
				// fmt.Println(testEq(rawEnergyParams.interior2x2LoopEnergy37C, interior2x2LoopEnergy37, 6))
			case "int22_enthalpies":
				rawEnergyParams.interior2x2LoopEnthalpy = parseParamValues(scanner, nbPairs-1, nbPairs-1, nbNucleobase, nbNucleobase, nbNucleobase, nbNucleobase).([][][][][][]int)
				rawEnergyParams.interior2x2LoopEnthalpy = addPreOffset(rawEnergyParams.interior2x2LoopEnthalpy, 1, 1, 1, 1, 1, 1).([][][][][][]int)
				rawEnergyParams.interior2x2LoopEnthalpy = addPostOffset(rawEnergyParams.interior2x2LoopEnthalpy, 1, 1, 0, 0, 0, 0).([][][][][][]int)
				rawEnergyParams.interior2x2LoopEnthalpy = updateNST(rawEnergyParams.interior2x2LoopEnthalpy)
				// fmt.Println(testEq(rawEnergyParams.interior2x2LoopEnthalpy, interior2x2LoopEnthalpy, 6))
			case "dangle5":
				rawEnergyParams.dangle5Energy37C = parseParamValues(scanner, nbPairs, nbNucleobase+1).([][]int)
				rawEnergyParams.dangle5Energy37C = addPreOffset(rawEnergyParams.dangle5Energy37C, 1, 0).([][]int)
				// fmt.Println(testEq(rawEnergyParams.dangle5Energy37C, dangle5Energy37C, 2))
			case "dangle5_enthalpies":
				rawEnergyParams.dangle5Enthalpy = parseParamValues(scanner, nbPairs, nbNucleobase+1).([][]int)
				rawEnergyParams.dangle5Enthalpy = addPreOffset(rawEnergyParams.dangle5Enthalpy, 1, 0).([][]int)
				// fmt.Println(testEq(rawEnergyParams.dangle5Enthalpy, dangle5Enthalpy, 2))
			case "dangle3":
				rawEnergyParams.dangle3Energy37C = parseParamValues(scanner, nbPairs, nbNucleobase+1).([][]int)
				rawEnergyParams.dangle3Energy37C = addPreOffset(rawEnergyParams.dangle3Energy37C, 1, 0).([][]int)
				// fmt.Println(testEq(rawEnergyParams.dangle3Energy37C, dangle3Energy37C, 2))
			case "dangle3_enthalpies":
				rawEnergyParams.dangle3Enthalpy = parseParamValues(scanner, nbPairs, nbNucleobase+1).([][]int)
				rawEnergyParams.dangle3Enthalpy = addPreOffset(rawEnergyParams.dangle3Enthalpy, 1, 0).([][]int)
				// fmt.Println(testEq(rawEnergyParams.dangle3Enthalpy, dangle3Enthalpy, 2))

			case "ML_params":
				mlParams := parseParamValues(scanner, 6).([]int)
				rawEnergyParams.setMultiLoopParams(mlParams)
			case "NINIO":
				ninioParams := parseParamValues(scanner, 3).([]int)
				rawEnergyParams.setNinioParams(ninioParams)
			case "Triloops":
				rawEnergyParams.triLoops,
					rawEnergyParams.triLoopEnergy37C,
					rawEnergyParams.triLoopEnthalpy = parseTriTetraHexaLoopParams(scanner)
			case "Tetraloops":
				rawEnergyParams.tetraLoops,
					rawEnergyParams.tetraLoopEnergy37C,
					rawEnergyParams.tetraLoopEnthalpy = parseTriTetraHexaLoopParams(scanner)
			case "Hexaloops":
				rawEnergyParams.hexaLoops,
					rawEnergyParams.hexaLoopEnergy37C,
					rawEnergyParams.hexaLoopEnthalpy = parseTriTetraHexaLoopParams(scanner)
			case "Misc":
				miscParams := parseParamValues(scanner, 4).([]int)
				rawEnergyParams.terminalAU37C = miscParams[2]
				rawEnergyParams.terminalAUEnthalpy = miscParams[3]
			}
		}
		// else { // skip all other lines
		// 	if len(strings.TrimSpace(removeComments(line))) == 0 {
		// 		fmt.Println(line)
		// 	}
		// }
		line, lineAvailable = readLine(scanner)
	}

	// fmt.Println("and the answer is...")
	// fmt.Println(reflect.DeepEqual(turner2004Params(), rawEnergyParams))
	return
}

func parseTriTetraHexaLoopParams(scanner *bufio.Scanner) (loops map[string]bool, energies map[string]int, enthalpies map[string]int) {
	loops, energies, enthalpies = make(map[string]bool), make(map[string]int), make(map[string]int)
	line, lineAvailable := readLine(scanner)
	for lineAvailable {
		if len(strings.TrimSpace(line)) == 0 {
			// blank line encountered
			return
		} else {
			line = removeComments(line)

			if len(strings.TrimSpace(line)) != 0 {
				values := strings.Fields(line)
				if len(values) != 3 {
					panic(fmt.Sprintf("encountered incorrect number of values. expected 3, got %v", len(values)))
				}
				loop, energy, enthalpy := values[0], parseInt(values[1]), parseInt(values[2])
				loops[loop] = true
				energies[loop] = energy
				enthalpies[loop] = enthalpy
			} // else line only contain comments so continue reading

			line, lineAvailable = readLine(scanner)
		}
	}
	return
}

func printDims(arr interface{}, dims int) {
	switch dims {
	case 0:
		fmt.Println()
	case 1:
		array := arr.([]int)
		fmt.Printf("%v ", len(array))
		printDims(array[0], dims-1)
	case 2:
		array := arr.([][]int)
		fmt.Printf("%v ", len(array))
		printDims(array[0], dims-1)
	case 3:
		array := arr.([][][]int)
		fmt.Printf("%v ", len(array))
		printDims(array[0], dims-1)
	case 4:
		array := arr.([][][][]int)
		fmt.Printf("%v ", len(array))
		printDims(array[0], dims-1)
	case 5:
		array := arr.([][][][][]int)
		fmt.Printf("%v ", len(array))
		printDims(array[0], dims-1)
	case 6:
		array := arr.([][][][][][]int)
		fmt.Printf("%v ", len(array))
		printDims(array[0], dims-1)
	}
}

// returns a slice of length numINF with all values set to INF
func infSlice(numINF int) (ret []int) {
	ret = make([]int, numINF)
	for i := 0; i < numINF; i++ {
		ret[i] = INF
	}
	return
}

// adds numINF `INF`s to the front of a slice
func prependInfsToSlice(slice []int, numINF int) []int {
	return append(infSlice(numINF), slice...)
}

func appendInfsToSlice(slice []int, numINF int) []int {
	return append(slice, infSlice(numINF)...)
}

func addPreOffset(values interface{}, dims ...int) interface{} {
	switch len(dims) {
	case 0:
		panic("invalid number of dims passed to parseParamValues")
	case 1:
		return prependInfsToSlice(values.([]int), dims[0])
	case 2:
		return addPreOffset2Dim(values.([][]int), dims[0], dims[1])
	case 3:
		return addPreOffset3Dim(values.([][][]int), dims[0], dims[1], dims[2])
	case 4:
		return addPreOffset4Dim(values.([][][][]int), dims[0], dims[1], dims[2], dims[3])
	case 5:
		return addPreOffset5Dim(values.([][][][][]int), dims[0], dims[1], dims[2], dims[3], dims[4])
	case 6:
		return addPreOffset6Dim(values.([][][][][][]int), dims[0], dims[1], dims[2], dims[3], dims[4], dims[5])
	}
	return nil
}

func addPreOffset2Dim(values [][]int, dim1Offset, dim2Offset int) (ret [][]int) {
	currLenDim1 := len(values)
	newLenDim1 := currLenDim1 + dim1Offset
	ret = make([][]int, 0, newLenDim1)

	var newLenDim2 int = dim2Offset
	if currLenDim1 > 0 {
		// add check to ensure we don't index an empty slice
		newLenDim2 += len(values[0])
	}
	for i := 0; i < dim1Offset; i++ {
		infSlice := infSlice(newLenDim2)
		ret = append(ret, infSlice)
	}

	for i := 0; i < currLenDim1; i++ {
		ret = append(ret, prependInfsToSlice(values[i], dim2Offset))
	}
	return
}

func addPreOffset3Dim(values [][][]int, dim1Offset, dim2Offset, dim3Offset int) (ret [][][]int) {
	currLenDim1 := len(values)
	newLenDim1 := currLenDim1 + dim1Offset
	ret = make([][][]int, 0, newLenDim1)

	var newLenDim2 int = dim2Offset
	var currLenDim2 int = 0
	if currLenDim1 > 0 {
		// add check to ensure we don't index an empty slice
		currLenDim2 = len(values[0])
		newLenDim2 += currLenDim2
	}
	var newLenDim3 int = dim3Offset
	if currLenDim2 > 0 {
		// add check to ensure we don't index an empty slice
		newLenDim3 += len(values[0][0])
	}
	for i := 0; i < dim1Offset; i++ {
		var infMatrix [][]int
		infMatrix = addPreOffset2Dim(infMatrix, newLenDim2, newLenDim3)
		ret = append(ret, infMatrix)
	}

	for i := 0; i < currLenDim1; i++ {
		ret = append(ret, addPreOffset2Dim(values[i], dim2Offset, dim3Offset))
	}
	return
}

func addPreOffset4Dim(values [][][][]int, dim1Offset, dim2Offset, dim3Offset, dim4Offset int) (ret [][][][]int) {
	currLenDim1 := len(values)
	newLenDim1 := currLenDim1 + dim1Offset
	ret = make([][][][]int, 0, newLenDim1)

	var newLenDim2 int = dim2Offset
	var currLenDim2 int = 0
	if currLenDim1 > 0 {
		// add check to ensure we don't index an empty slice
		currLenDim2 = len(values[0])
		newLenDim2 += currLenDim2
	}
	var newLenDim3 int = dim3Offset
	var currLenDim3 int = 0
	if currLenDim2 > 0 {
		// add check to ensure we don't index an empty slice
		currLenDim3 = len(values[0][0])
		newLenDim3 += currLenDim3
	}
	var newLenDim4 int = dim4Offset
	if currLenDim3 > 0 {
		// add check to ensure we don't index an empty slice
		newLenDim4 += len(values[0][0][0])
	}
	for i := 0; i < dim1Offset; i++ {
		var infMatrix [][][]int
		infMatrix = addPreOffset3Dim(infMatrix, newLenDim2, newLenDim3, newLenDim4)
		ret = append(ret, infMatrix)
	}

	for i := 0; i < len(values); i++ {
		ret = append(ret, addPreOffset3Dim(values[i], dim2Offset, dim3Offset, dim4Offset))
	}
	return
}

func addPreOffset5Dim(values [][][][][]int, dim1Offset, dim2Offset, dim3Offset, dim4Offset, dim5Offset int) (ret [][][][][]int) {
	currLenDim1 := len(values)
	newLenDim1 := currLenDim1 + dim1Offset
	ret = make([][][][][]int, 0, newLenDim1)

	var newLenDim2 int = dim2Offset
	var currLenDim2 int = 0
	if currLenDim1 > 0 {
		// add check to ensure we don't index an empty slice
		currLenDim2 = len(values[0])
		newLenDim2 += currLenDim2
	}
	var newLenDim3 int = dim3Offset
	var currLenDim3 int = 0
	if currLenDim2 > 0 {
		// add check to ensure we don't index an empty slice
		currLenDim3 = len(values[0][0])
		newLenDim3 += currLenDim3
	}
	var newLenDim4 int = dim4Offset
	var currLenDim4 int = 0
	if currLenDim3 > 0 {
		// add check to ensure we don't index an empty slice
		currLenDim4 = len(values[0][0][0])
		newLenDim4 += currLenDim4
	}
	var newLenDim5 int = dim5Offset
	if currLenDim4 > 0 {
		// add check to ensure we don't index an empty slice
		newLenDim5 += len(values[0][0][0][0])
	}
	for i := 0; i < dim1Offset; i++ {
		var infMatrix [][][][]int
		infMatrix = addPreOffset4Dim(infMatrix, newLenDim2, newLenDim3, newLenDim4, newLenDim5)
		ret = append(ret, infMatrix)
	}

	for i := 0; i < len(values); i++ {
		ret = append(ret, addPreOffset4Dim(values[i], dim2Offset, dim3Offset, dim4Offset, dim5Offset))
	}
	return
}

func addPreOffset6Dim(values [][][][][][]int, dim1Offset, dim2Offset, dim3Offset, dim4Offset, dim5Offset, dim6Offset int) (ret [][][][][][]int) {
	currLenDim1 := len(values)
	newLenDim1 := currLenDim1 + dim1Offset
	ret = make([][][][][][]int, 0, newLenDim1)

	var newLenDim2 int = dim2Offset
	var currLenDim2 int = 0
	if currLenDim1 > 0 {
		// add check to ensure we don't index an empty slice
		currLenDim2 = len(values[0])
		newLenDim2 += currLenDim2
	}
	var newLenDim3 int = dim3Offset
	var currLenDim3 int = 0
	if currLenDim2 > 0 {
		// add check to ensure we don't index an empty slice
		currLenDim3 = len(values[0][0])
		newLenDim3 += currLenDim3
	}
	var newLenDim4 int = dim4Offset
	var currLenDim4 int = 0
	if currLenDim3 > 0 {
		// add check to ensure we don't index an empty slice
		currLenDim4 = len(values[0][0][0])
		newLenDim4 += currLenDim4
	}
	var newLenDim5 int = dim5Offset
	var currLenDim5 int = 0
	if currLenDim4 > 0 {
		// add check to ensure we don't index an empty slice
		currLenDim5 = len(values[0][0][0][0])
		newLenDim5 += currLenDim5
	}
	var newLenDim6 int = dim6Offset
	if currLenDim5 > 0 {
		// add check to ensure we don't index an empty slice
		newLenDim6 += len(values[0][0][0][0][0])
	}
	for i := 0; i < dim1Offset; i++ {
		var infMatrix [][][][][]int
		infMatrix = addPreOffset5Dim(infMatrix, newLenDim2, newLenDim3, newLenDim4, newLenDim5, newLenDim6)
		ret = append(ret, infMatrix)
	}

	for i := 0; i < len(values); i++ {
		ret = append(ret, addPreOffset5Dim(values[i], dim2Offset, dim3Offset, dim4Offset, dim5Offset, dim6Offset))
	}
	return
}

// func readLineIntoIntArray(line string, arraySize int) ([]int, error) {
// 	lineSlice := readLineIntoSlice(line)
// 	var slice [4]int
// 	return slice
// }

func parseParamValues(scanner *bufio.Scanner, dims ...int) interface{} {
	switch len(dims) {
	case 0:
		panic("invalid number of dims passed to parseParamValues")
	case 1:
		return parseUntilEnoughValuesIntoSlice(scanner, dims[0])
	case 2:
		return parseParamValuesInto2DimSlice(scanner, dims[0], dims[1])
	case 3:
		return parseParamValuesInto3DimSlice(scanner, dims[0], dims[1], dims[2])
	case 4:
		return parseParamValuesInto4DimSlice(scanner, dims[0], dims[1], dims[2], dims[3])
	case 5:
		return parseParamValuesInto5DimSlice(scanner, dims[0], dims[1], dims[2], dims[3], dims[4])
	case 6:
		return parseParamValuesInto6DimSlice(scanner, dims[0], dims[1], dims[2], dims[3], dims[4], dims[5])
	}
	return nil
}

func parseUntilEnoughValuesIntoSlice(scanner *bufio.Scanner, numValuesToParse int) (ret []int) {
	totalValuesParsed := 0
	ret = make([]int, 0, numValuesToParse)
	for totalValuesParsed < numValuesToParse {
		parsedParamsSlice := parseLineIntoSlice(scanner)
		ret = append(ret, parsedParamsSlice...)
		totalValuesParsed += len(parsedParamsSlice)
	}
	if totalValuesParsed > numValuesToParse {
		panic("parsed too many values")
	}
	return
}

func addPostOffset(values interface{}, dims ...int) interface{} {
	switch len(dims) {
	case 0:
		panic("invalid number of dims passed to parseParamValues")
	case 1:
		return appendInfsToSlice(values.([]int), dims[0])
	case 2:
		return addPostOffset2Dim(values.([][]int), dims[0], dims[1])
	case 3:
		return addPostOffset3Dim(values.([][][]int), dims[0], dims[1], dims[2])
	case 4:
		return addPostOffset4Dim(values.([][][][]int), dims[0], dims[1], dims[2], dims[3])
	case 5:
		return addPostOffset5Dim(values.([][][][][]int), dims[0], dims[1], dims[2], dims[3], dims[4])
	case 6:
		return addPostOffset6Dim(values.([][][][][][]int), dims[0], dims[1], dims[2], dims[3], dims[4], dims[5])
	}
	return nil
}

func addPostOffset2Dim(values [][]int, dim1Offset, dim2Offset int) (ret [][]int) {
	currLenDim1 := len(values)
	newLenDim1 := currLenDim1 + dim1Offset
	ret = make([][]int, 0, newLenDim1)

	var newLenDim2 int = dim2Offset
	if currLenDim1 > 0 {
		// add check to ensure we don't index an empty slice
		newLenDim2 += len(values[0])
	}

	for i := 0; i < currLenDim1; i++ {
		ret = append(ret, appendInfsToSlice(values[i], dim2Offset))
	}

	for i := 0; i < dim1Offset; i++ {
		infSlice := infSlice(newLenDim2)
		ret = append(ret, infSlice)
	}

	return
}

func addPostOffset3Dim(values [][][]int, dim1Offset, dim2Offset, dim3Offset int) (ret [][][]int) {
	currLenDim1 := len(values)
	newLenDim1 := currLenDim1 + dim1Offset
	ret = make([][][]int, 0, newLenDim1)

	var newLenDim2 int = dim2Offset
	var currLenDim2 int = 0
	if currLenDim1 > 0 {
		// add check to ensure we don't index an empty slice
		currLenDim2 = len(values[0])
		newLenDim2 += currLenDim2
	}
	var newLenDim3 int = dim3Offset
	if currLenDim2 > 0 {
		// add check to ensure we don't index an empty slice
		newLenDim3 += len(values[0][0])
	}

	for i := 0; i < currLenDim1; i++ {
		ret = append(ret, addPostOffset2Dim(values[i], dim2Offset, dim3Offset))
	}

	for i := 0; i < dim1Offset; i++ {
		var infMatrix [][]int
		infMatrix = addPostOffset2Dim(infMatrix, newLenDim2, newLenDim3)
		ret = append(ret, infMatrix)
	}
	return
}

func addPostOffset4Dim(values [][][][]int, dim1Offset, dim2Offset, dim3Offset, dim4Offset int) (ret [][][][]int) {
	currLenDim1 := len(values)
	newLenDim1 := currLenDim1 + dim1Offset
	ret = make([][][][]int, 0, newLenDim1)

	var newLenDim2 int = dim2Offset
	var currLenDim2 int = 0
	if currLenDim1 > 0 {
		// add check to ensure we don't index an empty slice
		currLenDim2 = len(values[0])
		newLenDim2 += currLenDim2
	}
	var newLenDim3 int = dim3Offset
	var currLenDim3 int = 0
	if currLenDim2 > 0 {
		// add check to ensure we don't index an empty slice
		currLenDim3 = len(values[0][0])
		newLenDim3 += currLenDim3
	}
	var newLenDim4 int = dim4Offset
	if currLenDim3 > 0 {
		// add check to ensure we don't index an empty slice
		newLenDim4 += len(values[0][0][0])
	}

	for i := 0; i < currLenDim1; i++ {
		ret = append(ret, addPostOffset3Dim(values[i], dim2Offset, dim3Offset, dim4Offset))
	}

	for i := 0; i < dim1Offset; i++ {
		var infMatrix [][][]int
		infMatrix = addPostOffset3Dim(infMatrix, newLenDim2, newLenDim3, newLenDim4)
		ret = append(ret, infMatrix)
	}

	return
}

func addPostOffset5Dim(values [][][][][]int, dim1Offset, dim2Offset, dim3Offset, dim4Offset, dim5Offset int) (ret [][][][][]int) {
	currLenDim1 := len(values)
	newLenDim1 := currLenDim1 + dim1Offset
	ret = make([][][][][]int, 0, newLenDim1)

	var newLenDim2 int = dim2Offset
	var currLenDim2 int = 0
	if currLenDim1 > 0 {
		// add check to ensure we don't index an empty slice
		currLenDim2 = len(values[0])
		newLenDim2 += currLenDim2
	}
	var newLenDim3 int = dim3Offset
	var currLenDim3 int = 0
	if currLenDim2 > 0 {
		// add check to ensure we don't index an empty slice
		currLenDim3 = len(values[0][0])
		newLenDim3 += currLenDim3
	}
	var newLenDim4 int = dim4Offset
	var currLenDim4 int = 0
	if currLenDim3 > 0 {
		// add check to ensure we don't index an empty slice
		currLenDim4 = len(values[0][0][0])
		newLenDim4 += currLenDim4
	}
	var newLenDim5 int = dim5Offset
	if currLenDim4 > 0 {
		// add check to ensure we don't index an empty slice
		newLenDim5 += len(values[0][0][0][0])
	}

	for i := 0; i < currLenDim1; i++ {
		ret = append(ret, addPostOffset4Dim(values[i], dim2Offset, dim3Offset, dim4Offset, dim5Offset))
	}

	for i := 0; i < dim1Offset; i++ {
		var infMatrix [][][][]int
		infMatrix = addPostOffset4Dim(infMatrix, newLenDim2, newLenDim3, newLenDim4, newLenDim5)
		ret = append(ret, infMatrix)
	}

	return
}

func addPostOffset6Dim(values [][][][][][]int, dim1Offset, dim2Offset, dim3Offset, dim4Offset, dim5Offset, dim6Offset int) (ret [][][][][][]int) {
	currLenDim1 := len(values)
	newLenDim1 := currLenDim1 + dim1Offset
	ret = make([][][][][][]int, 0, newLenDim1)

	var newLenDim2 int = dim2Offset
	var currLenDim2 int = 0
	if currLenDim1 > 0 {
		// add check to ensure we don't index an empty slice
		currLenDim2 = len(values[0])
		newLenDim2 += currLenDim2
	}
	var newLenDim3 int = dim3Offset
	var currLenDim3 int = 0
	if currLenDim2 > 0 {
		// add check to ensure we don't index an empty slice
		currLenDim3 = len(values[0][0])
		newLenDim3 += currLenDim3
	}
	var newLenDim4 int = dim4Offset
	var currLenDim4 int = 0
	if currLenDim3 > 0 {
		// add check to ensure we don't index an empty slice
		currLenDim4 = len(values[0][0][0])
		newLenDim4 += currLenDim4
	}
	var newLenDim5 int = dim5Offset
	var currLenDim5 int = 0
	if currLenDim4 > 0 {
		// add check to ensure we don't index an empty slice
		currLenDim5 = len(values[0][0][0][0])
		newLenDim5 += currLenDim5
	}
	var newLenDim6 int = dim6Offset
	if currLenDim5 > 0 {
		// add check to ensure we don't index an empty slice
		newLenDim6 += len(values[0][0][0][0][0])
	}

	for i := 0; i < currLenDim1; i++ {
		ret = append(ret, addPostOffset5Dim(values[i], dim2Offset, dim3Offset, dim4Offset, dim5Offset, dim6Offset))
	}

	for i := 0; i < dim1Offset; i++ {
		var infMatrix [][][][][]int
		infMatrix = addPostOffset5Dim(infMatrix, newLenDim2, newLenDim3, newLenDim4, newLenDim5, newLenDim6)
		ret = append(ret, infMatrix)
	}

	return
}

// func parseUntilEnoughValuesIntoSlice(scanner *bufio.Scanner, numValuesToParse int) (ret []int) {
// 	scanner.Split(bufio.ScanWords)
// 	totalValuesParsed := 0
// 	ret = make([]int, 0, numValuesToParse)
// 	for totalValuesParsed < numValuesToParse {
// 		scanner.Scan()
// 		token := scanner.Text()
// 		if token == "/*" {
// 			readUntilCommentEnd(scanner)
// 		} else {
// 			ret = append(ret, parseToken(token))
// 			totalValuesParsed++
// 		}
// 		// parsedParamsSlice := parseLineIntoSlice(scanner)
// 		// ret = append(ret, parsedParamsSlice...)
// 		// totalValuesParsed += len(parsedParamsSlice)
// 	}
// 	if totalValuesParsed > numValuesToParse {
// 		panic("parsed too many values")
// 	}
// 	return
// }

// func readUntilCommentEnd(scanner *bufio.Scanner) {
// 	scanner.Split(bufio.ScanWords)
// 	token := ""
// 	for token != "*/" {
// 		scanner.Scan()
// 		token = scanner.Text()
// 	}
// }

// func parseToken(token string) int {
// 	if token == "INF" {
// 		return INF
// 	} else {
// 		valueInt64, err := strconv.ParseInt(token, 10, 0)
// 		if err != nil {
// 			panic(err)
// 		}
// 		return int(valueInt64)
// 	}
// }

func parseInt(token string) int {
	valueInt64, err := strconv.ParseInt(token, 10, 0)
	if err != nil {
		panic(err)
	}
	return int(valueInt64)
}

func parseLineIntoSlice(scanner *bufio.Scanner) (paramArray []int) {
	line := readUncommentedLine(scanner)
	values := strings.Fields(line)
	// if len(values) < lenArray {
	// 	panic("too few values to parse into array")
	// }

	for _, value := range values {
		if value == "INF" {
			paramArray = append(paramArray, INF)
		} else {
			paramArray = append(paramArray, parseInt(value))
		}
	}

	return
}

func parseParamValuesInto2DimSlice(scanner *bufio.Scanner, lenDim1, lenDim2 int) (ret [][]int) {
	ret = make([][]int, lenDim1)
	for i := 0; i < lenDim1; i++ {
		ret[i] = parseUntilEnoughValuesIntoSlice(scanner, lenDim2)
	}

	return ret
}

func parseParamValuesInto3DimSlice(scanner *bufio.Scanner, lenDim1, lenDim2, lenDim3 int) (ret [][][]int) {
	ret = make([][][]int, lenDim1)
	for i := 0; i < lenDim1; i++ {
		ret[i] = parseParamValuesInto2DimSlice(scanner, lenDim2, lenDim3)
	}
	return
}

func parseParamValuesInto4DimSlice(scanner *bufio.Scanner, lenDim1, lenDim2, lenDim3, lenDim4 int) (ret [][][][]int) {
	ret = make([][][][]int, lenDim1)
	for i := 0; i < lenDim1; i++ {
		ret[i] = parseParamValuesInto3DimSlice(scanner, lenDim2, lenDim3, lenDim4)
	}
	return
}

func parseParamValuesInto5DimSlice(scanner *bufio.Scanner, lenDim1, lenDim2, lenDim3, lenDim4, lenDim5 int) (ret [][][][][]int) {
	ret = make([][][][][]int, lenDim1)
	for i := 0; i < lenDim1; i++ {
		ret[i] = parseParamValuesInto4DimSlice(scanner, lenDim2, lenDim3, lenDim4, lenDim5)
	}
	return
}

func parseParamValuesInto6DimSlice(scanner *bufio.Scanner, lenDim1, lenDim2, lenDim3, lenDim4, lenDim5, lenDim6 int) (ret [][][][][][]int) {
	ret = make([][][][][][]int, lenDim1)
	for i := 0; i < lenDim1; i++ {
		ret[i] = parseParamValuesInto5DimSlice(scanner, lenDim2, lenDim3, lenDim4, lenDim5, lenDim6)
	}
	return
}

/* update nonstandard nucleotide/basepair involved contributions for int22 */
// array[nbPairs + 1][nbPairs + 1][5][5][5][5]
func updateNST(arrayPointer [][][][][][]int) (array [][][][][][]int) {
	array = arrayPointer
	var max, max2, max3, max4, max5, max6 int
	/* get maxima for one nonstandard nucleotide */
	for i := 1; i < nbPairs; i++ {
		for j := 1; j < nbPairs; j++ {
			for k := 1; k < 5; k++ {
				for l := 1; l < 5; l++ {
					for m := 1; m < 5; m++ {
						max, max2, max3, max4 = -INF, -INF, -INF, -INF /* max of {CGAU} */
						for n := 1; n < 5; n++ {
							max = maxInt(max, array[i][j][k][l][m][n])
							max2 = maxInt(max2, array[i][j][k][l][n][m])
							max3 = maxInt(max3, array[i][j][k][n][l][m])
							max4 = maxInt(max4, array[i][j][n][k][l][m])
						}
						array[i][j][k][l][m][0] = max
						array[i][j][k][l][0][m] = max2
						array[i][j][k][0][l][m] = max3
						array[i][j][0][k][l][m] = max4
					}
				}
			}
		}
	}
	/* get maxima for two nonstandard nucleotides */
	for i := 1; i < nbPairs; i++ {
		for j := 1; j < nbPairs; j++ {
			for k := 1; k < 5; k++ {
				for l := 1; l < 5; l++ {
					max, max2, max3, max4, max5, max6 = -INF, -INF, -INF, -INF, -INF, -INF /* max of {CGAU} */
					for m := 1; m < 5; m++ {
						max = maxInt(max, array[i][j][k][l][m][0])
						max2 = maxInt(max2, array[i][j][k][m][0][l])
						max3 = maxInt(max3, array[i][j][m][0][k][l])
						max4 = maxInt(max4, array[i][j][0][k][l][m])
						max5 = maxInt(max5, array[i][j][0][k][m][l])
						max6 = maxInt(max6, array[i][j][k][0][l][m])
					}
					array[i][j][k][l][0][0] = max
					array[i][j][k][0][0][l] = max2
					array[i][j][0][0][k][l] = max3
					array[i][j][k][0][l][0] = max6
					array[i][j][0][k][0][l] = max5
					array[i][j][0][k][l][0] = max4
				}
			}
		}
	}
	/* get maxima for three nonstandard nucleotides */
	for i := 1; i < nbPairs; i++ {
		for j := 1; j < nbPairs; j++ {
			for k := 1; k < 5; k++ {
				max, max2, max3, max4 = -INF, -INF, -INF, -INF /* max of {CGAU} */
				for l := 1; l < 5; l++ {
					/* should be arbitrary where index l resides in last 3 possible locations */
					max = maxInt(max, array[i][j][k][l][0][0])
					max2 = maxInt(max2, array[i][j][0][k][l][0])
					max3 = maxInt(max3, array[i][j][0][0][k][l])
					max4 = maxInt(max4, array[i][j][0][0][l][k])
				}
				array[i][j][k][0][0][0] = max
				array[i][j][0][k][0][0] = max2
				array[i][j][0][0][k][0] = max3
				array[i][j][0][0][0][k] = max4
			}
		}
	}
	/* get maxima for 4 nonstandard nucleotides */
	for i := 1; i < nbPairs; i++ {
		for j := 1; j < nbPairs; j++ {
			max = -INF /* max of {CGAU} */
			for k := 1; k < 5; k++ {
				max = maxInt(max, array[i][j][k][0][0][0])
			}
			array[i][j][0][0][0][0] = max
		}
	}

	/*
	 * now compute contributions for nonstandard base pairs ...
	 * first, 1 nonstandard bp
	 */
	for i := 1; i < nbPairs; i++ {
		for k := 0; k < 5; k++ {
			for l := 0; l < 5; l++ {
				for m := 0; m < 5; m++ {
					for n := 0; n < 5; n++ {
						max, max2 = -INF, -INF
						for j := 1; j < nbPairs; j++ {
							max = maxInt(max, array[i][j][k][l][m][n])
							max2 = maxInt(max2, array[j][i][k][l][m][n])
						}
						array[i][nbPairs][k][l][m][n] = max
						array[nbPairs][i][k][l][m][n] = max2
					}
				}
			}
		}
	}

	/* now 2 nst base pairs */
	for k := 0; k < 5; k++ {
		for l := 0; l < 5; l++ {
			for m := 0; m < 5; m++ {
				for n := 0; n < 5; n++ {
					max = -INF
					for j := 1; j < nbPairs; j++ {
						max = maxInt(max, array[nbPairs][j][k][l][m][n])
					}
					array[nbPairs][nbPairs][k][l][m][n] = max
				}
			}
		}
	}
	return
}

func maxInt(a, b int) int {
	if a > b {
		return a
	}
	return b
}

// func rd_2dim_slice(content []string, line_no int, array []int, dim [2]int, shift [2]int) {
// 	var post [2]int = int{0, 0}
// 	delta_pre := shift[0] + shift[1]
// 	delta_post := post[0] + post[1]

// 	if delta_pre + delta_post == 0 {
// 		rd_1dim(content, line_no, array, dim[0]*dim[1], 0)
// 		return
// 	}

// 	for i := shift[0]; i < dim[0]-post[0]; i++ {
// 		rd_1dim_slice(content, line_no, array+(i*dim[1]), dim[1], shift[1], post[1])
// 	}
// }

func readUncommentedLine(scanner *bufio.Scanner) (ret string) {
	line, lineAvailable := readLine(scanner)
	for lineAvailable {
		line = removeComments(line)
		if len(strings.TrimSpace(line)) == 0 {
			// line only contain comments so continue reading
			line, lineAvailable = readLine(scanner)
		} else {
			ret = line
			return
		}
	}
	panic("no more lines to read")
}

func testEq1DimSlice(a, b []int) bool {
	if len(a) != len(b) {
		return false
	}
	for i := range a {
		if a[i] != b[i] {
			return false
		}
	}
	return true
}

func testEq2DimSlice(a, b [][]int) bool {
	if len(a) != len(b) {
		return false
	}
	for i := range a {
		if !testEq1DimSlice(a[i], b[i]) {
			return false
		}
	}
	return true
}

func testEq3DimSlice(a, b [][][]int) bool {
	if len(a) != len(b) {
		return false
	}
	for i := range a {
		if !testEq2DimSlice(a[i], b[i]) {
			return false
		}
	}
	return true
}

func testEq4DimSlice(a, b [][][][]int) bool {
	if len(a) != len(b) {
		return false
	}
	for i := range a {
		if !testEq3DimSlice(a[i], b[i]) {
			return false
		}
	}
	return true
}

func testEq5DimSlice(a, b [][][][][]int) bool {
	if len(a) != len(b) {
		return false
	}
	for i := range a {
		if !testEq4DimSlice(a[i], b[i]) {
			return false
		}
	}
	return true
}

func testEq6DimSlice(a, b [][][][][][]int) bool {
	if len(a) != len(b) {
		return false
	}
	for i := range a {
		if !testEq5DimSlice(a[i], b[i]) {
			return false
		}
	}
	return true
}

func testEq(a, b interface{}, dims int) bool {
	switch dims {
	case 1:
		return testEq1DimSlice(a.([]int), b.([]int))
	case 2:
		return testEq2DimSlice(a.([][]int), b.([][]int))
	case 3:
		return testEq3DimSlice(a.([][][]int), b.([][][]int))
	case 4:
		return testEq4DimSlice(a.([][][][]int), b.([][][][]int))
	case 5:
		return testEq5DimSlice(a.([][][][][]int), b.([][][][][]int))
	case 6:
		return testEq6DimSlice(a.([][][][][][]int), b.([][][][][][]int))
	}
	return false
}

func getUncommentedLines(scanner *bufio.Scanner, numLines int) (ret []string) {
	line, lineAvailable := readLine(scanner)
	numUncommentedLines := 0

	for lineAvailable {
		line = removeComments(line)
		if len(strings.TrimSpace(line)) != 0 {
			ret = append(ret, line)
			numUncommentedLines++
		} // skip white space lines and full line comments
	}

	return
}

func removeComments(source string) string {
	commentStartIdx := strings.Index(source, "/*")
	for commentStartIdx != -1 {
		commentEndIdx := strings.Index(source, "*/")
		if commentEndIdx == -1 {
			panic("unclosed comment in parameter file")
		}

		source = source[:commentStartIdx] + source[commentEndIdx+2:]

		commentStartIdx = strings.Index(source, "/*")
	}
	return source
}

// func rd_1dim_slice(char **content, size_t *line_no, int *array, int dim, int shift) {
// 	var cp string
// 	cp = get_array1(content, line_no, array+shift, dim-shift)
// 	if cp {
// 		vrna_message_error("\nrd_1dim: %s", cp)
// 		exit(1)
// 	}
// }

// func get_array1(content []string, line_no int, arr []int, numLines int) string {
//   var i, p, pos, pp, r, last int
// 	var line string
// 	var buf [16]rune

//   i := 0
// 	last := 0
//   for i < numLines {
//     line = content[(*line_no)++]
//     if (!line)
//       vrna_message_error("unexpected end of file in get_array1");

//     ignore_comment(line);
//     pos = 0;
//     while ((i < numLines) && (sscanf(line + pos, "%15s%n", buf, &pp) == 1)) {
//       pos += pp;
//       if (buf[0] == '*') {
//         i++;
//         continue;
//       } else if (buf[0] == 'x') {
//         /* should only be used for loop parameters */
//         if (i == 0)
//           vrna_message_error("can't extrapolate first value");

//         p = arr[last] + (int)(0.5 + lxc37 * log(((double)i) / (double)(last)));
//       } else if (strcmp(buf, "DEF") == 0) {
//         p = DEF;
//       } else if (strcmp(buf, "INF") == 0) {
//         p = INF;
//       } else if (strcmp(buf, "NST") == 0) {
//         p = NST;
//       } else {
//         r = sscanf(buf, "%d", &p);
//         if (r != 1) {
//           return line + pos;
//           vrna_message_error("can't interpret `%s' in get_array1", buf);
//           exit(1);
//         }

//         last = i;
//       }

//       arr[i++] = p;
//     }
//   }

//   return NULL;
// }

// ignore_comment(char *line)
// {
//   /*
//    * excise C style comments
//    * only one comment per line, no multiline comments
//    */
//   char *cp1, *cp2;

//   if ((cp1 = strstr(line, "/*"))) {
//     cp2 = strstr(cp1, "*/");
//     if (cp2 == NULL)
//       vrna_message_error("unclosed comment in parameter file");

//     /* can't use strcpy for overlapping strings */
//     for (cp2 += 2; *cp2 != '\0'; cp2++, cp1++)
//       *cp1 = *cp2;
//     *cp1 = '\0';
//   }

//   return;
// }
