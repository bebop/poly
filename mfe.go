package poly

import (
	"errors"
	"fmt"
	"log"
	"math"
	"strings"
)

/*
* The MFE calculation functionality has been taken directly from ViennaRNA with the follwing changes:
* - ViennaRNA includes the ability to specify a dangle model (more info available at [src/bin/RNAeval.ggo#L153](https://github.com/ViennaRNA/ViennaRNA/blob/d6fbaf2b33d65944f9249a91ed5ab4b3277f7d06/src/bin/RNAeval.ggo#L153)). This implementation keeps it simple and defaults to the value of -d2.
* - ViennaRNA includes the ability to specify hard and soft constraints (more info available from [their ppt explaining hard and soft constrains](http://benasque.org/2015rna/talks_contr/2713_lorenz_benasque_2015.pdf), [official docs](https://www.tbi.univie.ac.at/RNA/ViennaRNA/doc/html/group__constraints.html), [thier paper](https://almob.biomedcentral.com/articles/10.1186/s13015-016-0070-z)). This implementation keeps it simple and defaults to no hard or soft constraints.
* - ViennaRNA includes the ability to calculate the minimum free energy of co-folded sequences. This implementation keeps it simple and defaults to calculating the mfe of only single sequences.
 */

type foldCompound struct {
	length          int           // length of `sequence`
	params          *energyParams // The precomputed free energy contributions for each type of loop
	sequence        string        // The input sequence string
	encodedSequence []int         // Numerical encoding of the sequence (see encodeSequence for the nucleotide->int mapping)
	pairTable       []int         //

}

/**
 *  Calculate the free energy of an already folded RNA. Contributions on a per-loop base are logged.
 *
 *  This function allows for detailed energy evaluation of a given sequence/structure pair.
 *
 *  @param sequence         A RNA sequence
 *  @param structure        Secondary structure in dot-bracket notation
 *  @return                 The free energy of the input structure given the input sequence in kcal/mol
 */
func CalculateMFE(sequence, structure string) (float64, error) {
	lenSequence := len(sequence)
	lenStructure := len(structure)

	if lenSequence != lenStructure {
		return 0, fmt.Errorf("length of sequence (%v) != Length of structure (%v)", lenSequence, lenStructure)
	} else if lenStructure == 0 {
		return 0, errors.New("lengths of sequence and structure cannot be 0")
	}

	// convert DNA to RNA
	// sequence = strings.ToUpper(strings.ReplaceAll(sequence, "T", "U"))

	// fmt.Println(regexp.MatchString("[^ATCG]+", sequence))

	pairTable, err := pairTable(structure)
	if err != nil {
		return 0, err
	}

	fc := &foldCompound{
		length:          lenSequence,
		params:          get_scaled_params(),
		sequence:        sequence,
		encodedSequence: encodeSequence(sequence),
		pairTable:       pairTable,
	}

	energy := float64(evalFoldCompound(fc)) / 100.0

	return energy, nil
}

func rescaleDg(dG, dH, dT int) int {
	return (dH - (dH-dG)*dT)
}

const (
	nbPairs            int     = 7    /** The number of distinguishable base pairs */
	maxLoop            int     = 30   /** The maximum loop length */
	defaultTemperature float64 = 37.0 // Default temperature for structure prediction and free energy evaluation
)

/* vivek: No clue what these values correspond to. Why isn't the matrix symmetric?
* Created issue in original repo to get more info (https://github.com/ViennaRNA/ViennaRNA/issues/124)
 */
func defaultBasePairEncodedTypeMap() map[byte]map[byte]int {
	//   			  _  A  C  G  U
	// /* _ */ {0, 0, 0, 0, 0},
	// /* A */ {0, 0, 0, 0, 5},
	// /* C */ {0, 0, 0, 1, 0},
	// /* G */ {0, 0, 2, 0, 3},
	// /* U */ {0, 6, 0, 4, 0},

	nucelotideAEncodedTypeMap := map[byte]int{'U': 5}
	nucelotideCEncodedTypeMap := map[byte]int{'G': 1}
	nucelotideGEncodedTypeMap := map[byte]int{'C': 2, 'U': 3}
	nucelotideUEncodedTypeMap := map[byte]int{'A': 6, 'G': 4}

	pair := map[byte]map[byte]int{
		'A': nucelotideAEncodedTypeMap,
		'C': nucelotideCEncodedTypeMap,
		'G': nucelotideGEncodedTypeMap,
		'U': nucelotideUEncodedTypeMap,
	}

	return pair
}

// contains all the energy information needed for the free energy calculations
type energyParams struct {
	stack                                     [nbPairs + 1][nbPairs + 1]int
	hairpin                                   [31]int
	bulge                                     [maxLoop + 1]int
	internal_loop                             [maxLoop + 1]int
	mismatchExt                               [nbPairs + 1][5][5]int
	mismatchI                                 [nbPairs + 1][5][5]int
	mismatch1nI                               [nbPairs + 1][5][5]int
	mismatch23I                               [nbPairs + 1][5][5]int
	mismatchH                                 [nbPairs + 1][5][5]int
	mismatchM                                 [nbPairs + 1][5][5]int
	dangle5                                   [nbPairs + 1][5]int
	dangle3                                   [nbPairs + 1][5]int
	int11                                     [nbPairs + 1][nbPairs + 1][5][5]int
	int21                                     [nbPairs + 1][nbPairs + 1][5][5][5]int
	int22                                     [nbPairs + 1][nbPairs + 1][5][5][5][5]int
	ninio                                     [5]int
	lxc                                       float64
	MLbase, MLclosing, TerminalAU, DuplexInit int
	MLintern                                  [nbPairs + 1]int
	Tetraloop_E                               [200]int
	Tetraloops                                string
	Triloop_E                                 [40]int
	Triloops                                  string
	Hexaloop_E                                [40]int
	Hexaloops                                 string
	TripleC, MultipleCA, MultipleCB           int
	/* Integer representation of a base pair
	* vivek: This field was previously in the vrna_md_t struct. However since it
	* was the only field accessed from that struct, I added it here and removed
	* that struct
	 */
	basePairEncodedTypeMap map[byte]map[byte]int
}

/* vivek:
* returns a slice `pairTable` where `pairTable[i]`` returns the index of the
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

	// if len_seq > SHRT_MAX {
	// 	err := fmt.Sprintf("pairTable: Structure too long to be converted to pair table (n=%v, max=%v)", len_seq, SHRT_MAX)
	// 	return nil, errors.New(err)
	// }

	pairTable = make([]int, len_seq)

	// the characters of the opening and closing brackets
	var open_bracket, close_bracket byte = '(', ')'

	// keeps track of the indexes of open brackets. indexes of open brackets are pushed onto stack
	// and poped off when a closing bracket is encountered

	var open_bracket_idx_stack []int = make([]int, len_seq)
	var stack_idx int = 0

	// iterate through structure_char_slice and create pairTable
	for i := 0; i < len(structure); i++ {
		if structure[i] == open_bracket {
			// push index of open bracket onto stack
			open_bracket_idx_stack[stack_idx] = i
			stack_idx++
		} else if structure[i] == close_bracket {
			stack_idx--

			if stack_idx < 0 {
				return nil, fmt.Errorf("%v\nunbalanced brackets '%v%v' found while extracting base pairs", structure,
					open_bracket, close_bracket)
			}

			open_bracket_idx := open_bracket_idx_stack[stack_idx]
			// current index of one-indexed sequence
			pairTable[i] = open_bracket_idx
			pairTable[open_bracket_idx] = i
		} else {
			pairTable[i] = -1
		}
	}

	if stack_idx != 0 {
		return nil, fmt.Errorf("%v\nunbalanced brackets '%v%v' found while extracting base pairs", structure,
			open_bracket, close_bracket)
	}

	return pairTable, nil
}

/* vivek:
*
*
 */
func evalFoldCompound(fc *foldCompound) int {
	pairTable := fc.pairTable
	var i, length, energy int

	length = fc.length

	energy = externalLoopEnergy(fc)
	logExternalLoopEnergy(energy)

	for i = 0; i < length; i++ {
		if pairTable[i] == -1 {
			continue
		}

		energy += stack_energy(fc, i)
		i = pairTable[i]
	}

	return energy
}

func logExternalLoopEnergy(energy int) {
	log.Printf("External loop                           : %v\n", energy)
}

// encodes a sequence into its numerical representaiton
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

/*
* Calculate the energy contribution of
* stabilizing dangling-ends/mismatches
* for all stems branching off the exterior
* loop
 */
func externalLoopEnergy(fc *foldCompound) int {
	pairTable := fc.pairTable
	var energy int = 0
	length := fc.length
	pairOpenIdx := 0

	/* seek to opening base of first stem (pair) */
	for pairOpenIdx < length && pairTable[pairOpenIdx] == -1 {
		pairOpenIdx++
	}

	for pairOpenIdx < length {
		/* pairOpenIdx must have a pairing partner */
		pairCloseIdx := pairTable[pairOpenIdx]

		/* get type of base pair (pairOpenIdx, pairCloseIdx) */
		basePairType := encodeBasePairType(fc.sequence[pairOpenIdx],
			fc.sequence[pairCloseIdx],
			fc.params.basePairEncodedTypeMap)

		// vivek: cryptic variable names. 5 and 3 have to do with the 5' and 3' ends
		// of a DNA/RNA strand
		var fivePrimeMismatch, threePrimeMismatch int
		if pairOpenIdx > 0 {
			fivePrimeMismatch = fc.encodedSequence[pairOpenIdx-1]
		} else {
			fivePrimeMismatch = -1
		}

		if pairCloseIdx < length-1 {
			threePrimeMismatch = fc.encodedSequence[pairCloseIdx+1]
		} else {
			threePrimeMismatch = -1
		}

		energy += exteriorStemEnergy(basePairType, fivePrimeMismatch, threePrimeMismatch, fc.params)

		/* seek to the next stem */
		pairOpenIdx = pairCloseIdx + 1
		for pairOpenIdx < length && pairTable[pairOpenIdx] == -1 {
			pairOpenIdx++
		}
	}

	return energy
}

/**
 *  @brief  Evaluate a stem branching off the exterior loop
 *
 *  Given a base pair (i,j) encoded by type, compute the energy contribution
 *  including dangling-end/terminal-mismatch contributions. Instead of returning the
 *  energy contribution per-se, this function returns the corresponding Boltzmann factor.
 *  If either of the adjacent nucleotides i - 1 (`fivePrimeAdjacent`) and j + 1 (`threePrimeAdjacent`) must not
 *  contribute stacking energy, the corresponding encoding must be -1.
 *
 *  @param  basePairType  The base pair type encoding (see `basePairType()`)
 *  @param  fivePrimeAdjacent   The encoded nucleotide directly adjacent at the 5' side of the base pair (may be -1 if index of i is 0)
 *  @param  threePrimeAdjacent   The encoded nucleotide directly adjacent at the 3' side of the base pair (may be -1 if index of j is len(sequence) - 1)
 *  @param  p     The pre-computed energy parameters
 *  @return       The energy contribution of the introduced exterior-loop stem
 */
func exteriorStemEnergy(encodedBasePairType int, fivePrimeAdjacent int, threePrimeAdjacent int, energyParams *energyParams) int {
	var energy int = 0

	if fivePrimeAdjacent >= 0 && threePrimeAdjacent >= 0 {
		energy += energyParams.mismatchExt[encodedBasePairType][fivePrimeAdjacent][threePrimeAdjacent]
	} else if fivePrimeAdjacent >= 0 {
		// `j` is the last nucleotide of the sequence
		energy += energyParams.dangle5[encodedBasePairType][fivePrimeAdjacent]
	} else if threePrimeAdjacent >= 0 {
		// `i` is the first nucleotide of the sequence
		energy += energyParams.dangle3[encodedBasePairType][threePrimeAdjacent]
	}

	if encodedBasePairType > 2 {
		energy += energyParams.TerminalAU
	}

	return energy
}

/* vivek:
*	what is i? what is stack_energy?
 */
func stack_energy(fc *foldCompound, pairOpenIdx int) int {
	/* recursively calculate energy of substructure enclosed by (pairOpenIdx, pairCloseIdx) */
	pairTable := fc.pairTable
	energy := 0
	pairCloseIdx := pairTable[pairOpenIdx]
	sequence := fc.sequence

	if fc.params.basePairEncodedTypeMap[fc.sequence[pairOpenIdx]][fc.sequence[pairCloseIdx]] == 0 {
		panic(fmt.Sprintf("bases %v and %v (%v%v) can't pair!",
			pairOpenIdx, pairCloseIdx,
			string(sequence[pairOpenIdx]),
			string(sequence[pairCloseIdx])))
	}

	n5p_iterator := pairOpenIdx
	n3p_iterator := pairCloseIdx

	for n5p_iterator < n3p_iterator {
		/* process all stacks and interior loops */

		// seek to opening bracket from 5' end
		n5p_iterator++
		for pairTable[n5p_iterator] == -1 {
			n5p_iterator++
		}

		// seek to closing bracket from 3' end
		n3p_iterator--
		for pairTable[n3p_iterator] == -1 {
			n3p_iterator--
		}

		if (pairTable[n3p_iterator] != n5p_iterator) || (n5p_iterator > n3p_iterator) {
			break
		}

		// vivek: should this be basePairType[n5p_iterator][n3p_iterator] or is the current order correct?
		if fc.params.basePairEncodedTypeMap[fc.sequence[n3p_iterator]][fc.sequence[n5p_iterator]] == 0 {
			panic(fmt.Sprintf("bases %d and %d (%c%c) can't pair!",
				n5p_iterator, n3p_iterator,
				sequence[n5p_iterator],
				sequence[n3p_iterator]))
		}

		ee := eval_int_loop(fc, pairOpenIdx, pairCloseIdx, n5p_iterator, n3p_iterator)
		energy += ee

		vrna_cstr_print_eval_int_loop(pairOpenIdx, pairCloseIdx,
			sequence[pairOpenIdx], sequence[pairCloseIdx],
			n5p_iterator, n3p_iterator,
			sequence[n5p_iterator], sequence[n3p_iterator], ee)

		pairOpenIdx = n5p_iterator
		pairCloseIdx = n3p_iterator
	} /* end while */

	/* n5p_iterator, n3p_iterator don't pair must have found hairpin or multiloop */

	if n5p_iterator > n3p_iterator {
		/* hairpin */

		ee := vrna_eval_hp_loop(fc, pairOpenIdx, pairCloseIdx)
		energy += ee

		vrna_cstr_print_eval_hp_loop(pairOpenIdx, pairCloseIdx, sequence[pairOpenIdx], sequence[pairCloseIdx], ee)

		return energy
	}

	/* (pairOpenIdx, pairCloseIdx) is exterior pair of multiloop */
	for n5p_iterator < pairCloseIdx {
		/* add up the contributions of the substructures of the ML */
		energy += stack_energy(fc, n5p_iterator)
		n5p_iterator = pairTable[n5p_iterator]
		/* search for next base pair in multiloop */
		n5p_iterator++
		for pairTable[n5p_iterator] == -1 {
			n5p_iterator++
		}
	}

	var ee int = 0
	var err error

	ee, err = energy_of_ml_pairTable(fc, pairOpenIdx, pairTable)
	if err != nil {
		panic(err)
	}

	energy += ee

	vrna_cstr_print_eval_mb_loop(pairOpenIdx, pairCloseIdx, sequence[pairOpenIdx], sequence[pairCloseIdx], ee)

	return energy
}

func vrna_cstr_print_eval_int_loop(i, j int, si, sj byte, k, l int, sk, sl byte, energy int) {
	log.Printf("Interior loop (%v,%v) %v%v; (%v,%v) %v%v: %v\n",
		i+1, j+1,
		string(si), string(sj),
		k+1, l+1,
		string(sk), string(sl),
		energy)
}

func vrna_cstr_print_eval_hp_loop(i, j int, si, sj byte, energy int) {
	log.Printf("Hairpin loop  (%v,%v) %v%v              : %v\n",
		i+1, j+1,
		string(si), string(sj),
		energy)
}

func vrna_cstr_print_eval_mb_loop(i, j int, si, sj byte, energy int) {
	log.Printf("Multi   loop (%v,%v) %v%v              : %v\n",
		i+1, j+1,
		string(si), string(sj),
		energy)
}

/**
 *  @brief Evaluate the free energy contribution of an interior loop with delimiting
 *  base pairs (i,j) and (k,l)
 *
 */
func eval_int_loop(fc *foldCompound, i, j, k, l int) int {
	var u1, u2 int
	energy := 0

	u1 = k - i - 1
	u2 = j - l - 1

	ij_basePairType := encodeBasePairType(fc.sequence[i], fc.sequence[j], fc.params.basePairEncodedTypeMap)
	kl_basePairType := encodeBasePairType(fc.sequence[l], fc.sequence[k], fc.params.basePairEncodedTypeMap)

	energy = E_IntLoop(u1, u2, ij_basePairType, kl_basePairType, fc.encodedSequence[i+1], fc.encodedSequence[j-1], fc.encodedSequence[k-1], fc.encodedSequence[l+1], fc.params)

	return energy
}

/* vivek:
*	Couldn't find any documentation about what `md.pair` is.
* Created an issue to get more info: https://github.com/ViennaRNA/ViennaRNA/issues/124
* retruns type of base pair (p,q)
 */
func encodeBasePairType(i, j byte, basePairEncodedTypeMap map[byte]map[byte]int) int {
	var encodedType int = basePairEncodedTypeMap[i][j]

	if encodedType == 0 {
		return 7
	} else {
		return encodedType
	}
}

/**
 *  Compute the Energy of an interior-loop
 *  This function computes the free energy of an interior-loop with the
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
 *  This general structure depicts an interior-loop that is closed by the base pair (X,Y).
 *  The enclosed base pair is (V,U) which leaves the unpaired bases a_1-a_n and b_1-b_n
 *  that constitute the loop.
 *  The base pair (X,Y) will be refered to as `closingBasePair`, and the base pair
 *  (U,V) will be refered to as `enclosedBasePair`.
 *  In this example, the length of the interior-loop is `n+m`
 *  where n or m may be 0 resulting in a bulge-loop or base pair stack.
 *  The mismatching nucleotides for the closing pair (X,Y) are:
 *  5'-mismatch: a_1
 *  3'-mismatch: b_m
 *  and for the enclosed base pair (V,U):
 *  5'-mismatch: b_1
 *  3'-mismatch: a_n
 *  @note Base pairs are always denoted in 5'->3' direction. Thus the `enclosedBasePair`
 *  must be 'turned arround' when evaluating the free energy of the interior-loop
 *
 *  @param  nbUnpairedLeftLoop      The size of the 'left'-loop (number of unpaired nucleotides)
 *  @param  nbUnpairedRightLoop     The size of the 'right'-loop (number of unpaired nucleotides)
 *  @param  closingBasePairType    	The encoded type of the closing base pair of the interior loop
 *  @param  enclosedBasePairType    The encoded type of the enclosed base pair
 *  @param  cbpFivePrimeMismatch     The 5'-mismatching nucleotide of the closing pair
 *  @param  cbpThreePrimeMismatch     The 3'-mismatching nucleotide of the closing pair
 *  @param  ebpThreePrimeMismatch     The 3'-mismatching nucleotide of the enclosed pair
 *  @param  ebpFivePrimeMismatch     The 5'-mismatching nucleotide of the enclosed pair
 *  @param  P       The datastructure containing scaled energy parameters
 *  @return The Free energy of the Interior-loop in dcal/mol
 */
func E_IntLoop(nbUnpairedLeftLoop, nbUnpairedRightLoop int,
	closingBasePairType, enclosedBasePairType int,
	cbpFivePrimeMismatch, cbpThreePrimeMismatch, ebpThreePrimeMismatch, ebpFivePrimeMismatch int, P *energyParams) int {
	/* compute energy of degree 2 loop (stack bulge or interior) */
	var nl, ns, u, energy int

	if nbUnpairedLeftLoop > nbUnpairedRightLoop {
		nl = nbUnpairedLeftLoop
		ns = nbUnpairedRightLoop
	} else {
		nl = nbUnpairedRightLoop
		ns = nbUnpairedLeftLoop
	}

	if nl == 0 {
		return P.stack[closingBasePairType][enclosedBasePairType] /* stack */
	}

	if ns == 0 {
		/* bulge */
		if nl <= maxLoop {
			energy = P.bulge[nl]
		} else {
			energy = P.bulge[30] + int(P.lxc*math.Log(float64(nl)/30.0))
		}

		if nl == 1 {
			energy += P.stack[closingBasePairType][enclosedBasePairType]
		} else {
			if closingBasePairType > 2 {
				energy += P.TerminalAU
			}

			if enclosedBasePairType > 2 {
				energy += P.TerminalAU
			}
		}

		return energy
	} else {
		/* interior loop */
		if ns == 1 {
			if nl == 1 {
				/* 1x1 loop */
				return P.int11[closingBasePairType][enclosedBasePairType][cbpFivePrimeMismatch][cbpThreePrimeMismatch]
			}

			if nl == 2 {
				/* 2x1 loop */
				if nbUnpairedLeftLoop == 1 {
					energy = P.int21[closingBasePairType][enclosedBasePairType][cbpFivePrimeMismatch][ebpFivePrimeMismatch][cbpThreePrimeMismatch]
				} else {
					energy = P.int21[enclosedBasePairType][closingBasePairType][ebpFivePrimeMismatch][cbpFivePrimeMismatch][ebpThreePrimeMismatch]
				}
				return energy
			} else {
				/* 1xn loop */
				if nl+1 <= maxLoop {
					energy = P.internal_loop[nl+1]
				} else {
					energy = P.internal_loop[30] + int(P.lxc*math.Log((float64(nl)+1.0)/30.0))
				}
				energy += min(MAX_NINIO, (nl-ns)*P.ninio[2])
				energy += P.mismatch1nI[closingBasePairType][cbpFivePrimeMismatch][cbpThreePrimeMismatch] + P.mismatch1nI[enclosedBasePairType][ebpFivePrimeMismatch][ebpThreePrimeMismatch]
				return energy
			}
		} else if ns == 2 {
			if nl == 2 {
				/* 2x2 loop */
				return P.int22[closingBasePairType][enclosedBasePairType][cbpFivePrimeMismatch][ebpThreePrimeMismatch][ebpFivePrimeMismatch][cbpThreePrimeMismatch]
			} else if nl == 3 {
				/* 2x3 loop */
				energy = P.internal_loop[5] + P.ninio[2]
				energy += P.mismatch23I[closingBasePairType][cbpFivePrimeMismatch][cbpThreePrimeMismatch] + P.mismatch23I[enclosedBasePairType][ebpFivePrimeMismatch][ebpThreePrimeMismatch]
				return energy
			}
		}

		{
			/* generic interior loop (no else here!)*/
			u = nl + ns
			if u <= maxLoop {
				energy = P.internal_loop[u]
			} else {
				energy = P.internal_loop[30] + int(P.lxc*math.Log(float64(u)/30.0))
			}

			energy += min(MAX_NINIO, (nl-ns)*P.ninio[2])

			energy += P.mismatchI[closingBasePairType][cbpFivePrimeMismatch][cbpThreePrimeMismatch] + P.mismatchI[enclosedBasePairType][ebpFivePrimeMismatch][ebpThreePrimeMismatch]
		}
	}

	return energy
}

/**
 *  @brief Evaluate free energy of a hairpin loop
 *
 *  @ingroup eval
 *
 *  @note This function is polymorphic! The provided #foldCompound may be of type
 *  #VRNA_FC_TYPE_SINGLE or #VRNA_FC_TYPE_COMPARATIVE
 *
 *  @param  fc  The #foldCompound for the particular energy evaluation
 *  @param  i   5'-position of the base pair
 *  @param  j   3'-position of the base pair
 *  @returns    Free energy of the hairpin loop closed by @f$ (i,j) @f$ in deka-kal/mol
 */
func vrna_eval_hp_loop(fc *foldCompound, pairOpenIdx, pairCloseIdx int) int {
	var u, e int
	var P *energyParams

	P = fc.params

	u = pairCloseIdx - pairOpenIdx - 1
	basePairType := encodeBasePairType(fc.sequence[pairOpenIdx], fc.sequence[pairCloseIdx], fc.params.basePairEncodedTypeMap)

	e = E_Hairpin(u, basePairType, fc.encodedSequence[pairOpenIdx+1], fc.encodedSequence[pairCloseIdx-1], fc.sequence[pairOpenIdx-1:], P)
	// fmt.Printf("hairpin e: %v\n", e)
	return e
}

/**
 *  @brief Compute the Energy of a hairpin-loop
 *
 *  To evaluate the free energy of a hairpin-loop, several parameters have to be known.
 *  A general hairpin-loop has this structure:<BR>
 *  <PRE>
 *        a3 a4
 *      a2     a5
 *      a1     a6
 *        X - Y
 *        |   |
 *        5'  3'
 *  </PRE>
 *  where X-Y marks the closing pair [e.g. a <B>(G,C)</B> pair]. The length of this loop is 6 as there are
 *  six unpaired nucleotides (a1-a6) enclosed by (X,Y). The 5' mismatching nucleotide is
 *  a1 while the 3' mismatch is a6. The nucleotide sequence of this loop is &quot;a1.a2.a3.a4.a5.a6&quot; <BR>
 *  @note The parameter sequence should contain the sequence of the loop in capital letters of the nucleic acid
 *  alphabet if the loop size is below 7. This is useful for unusually stable tri-, tetra- and hexa-loops
 *  which are treated differently (based on experimental data) if they are tabulated.
 *  @see scale_parameters()
 *  @see energyParams
 *  @warning Not (really) thread safe! A threadsafe implementation will replace this function in a future release!\n
 *  Energy evaluation may change due to updates in global variable "tetra_loop"
 *
 *  @param  size  The size of the loop (number of unpaired nucleotides)
 *  @param  type  The pair type of the base pair closing the hairpin
 *  @param  si1   The 5'-mismatching nucleotide
 *  @param  sj1   The 3'-mismatching nucleotide
 *  @param  string  The sequence of the loop (May be @p NULL, otherwise mst be at least @f$size + 2@f$ long)
 *  @param  P     The datastructure containing scaled energy parameters
 *  @return The Free energy of the Hairpin-loop in dcal/mol
 */
func E_Hairpin(size int, type_1 int, si1, sj1 int, sequence string, P *energyParams) int {
	var energy int

	if size <= 30 {
		energy = P.hairpin[size]
	} else {
		energy = P.hairpin[30] + int(P.lxc*math.Log(float64(size)/30.0))
	}

	if size < 3 {
		return energy /* should only be the case when folding alignments */
	}

	// if P.model_details.special_hp == 1 {
	var tl string
	var idx int
	if size == 4 {
		tl = sequence[:6]
		// vivek: this could be a point of failure. Maybe change the above to 7
		// memcpy(tl, sequence, sizeof(char) * 6);
		idx = strings.Index(string(P.Tetraloops), tl)
		if idx != -1 {
			return P.Tetraloop_E[idx/7]
		}
	} else if size == 6 {
		tl = sequence[:8]
		idx = strings.Index(string(P.Hexaloops), tl)
		if idx != -1 {
			return P.Hexaloop_E[idx/9]
		}
	} else if size == 3 {
		tl = sequence[:5]
		idx = strings.Index(string(P.Triloops), tl)
		if idx != -1 {
			return P.Triloop_E[idx/6]
		}

		if type_1 > 2 {
			return energy + P.TerminalAU
		} else {
			return energy
		}
	}
	// }

	energy += P.mismatchH[type_1][si1][sj1]

	return energy
}

/**
 *** pairOpenIdx is the 5'-base of the closing pair
 ***
 *** since each helix can coaxially stack with at most one of its
 *** neighbors we need an auxiliarry variable  cx_energy
 *** which contains the best energy given that the last two pairs stack.
 *** energy  holds the best energy given the previous two pairs do not
 *** stack (i.e. the two current helices may stack)
 *** We don't allow the last helix to stack with the first, thus we have to
 *** walk around the Loop twice with two starting points and take the minimum
 ***/
func energy_of_ml_pairTable(vc *foldCompound, pairOpenIdx int, pairTable []int) (int, error) {

	var energy int
	var pairCloseIdx, n5p_iterator_close_idx, num_unpaired_nucs, fivePrimeMismatch, threePrimeMismatch int

	length := int(vc.length)

	bonus := 0

	if pairOpenIdx >= pairTable[pairOpenIdx] {
		return INF, errors.New("energy_of_ml_pairTable: pairOpenIdx is not 5' base of a closing pair")
	}

	if pairOpenIdx == 0 {
		pairCloseIdx = length
	} else {
		pairCloseIdx = pairTable[pairOpenIdx]
	}

	/* init the variables */
	energy = 0
	num_unpaired_nucs = 0 /* the total number of unpaired nucleotides */
	n5p_iterator := pairOpenIdx + 1

	for n5p_iterator <= pairCloseIdx && pairTable[n5p_iterator] == -1 {
		n5p_iterator++
	}

	/* add bonus energies for first stretch of unpaired nucleotides */
	num_unpaired_nucs += n5p_iterator - pairOpenIdx - 1

	for n5p_iterator < pairCloseIdx {
		/* p must have a pairing partner */
		n5p_iterator_close_idx = pairTable[n5p_iterator]
		/* get type of base pair (p,q) */
		basePairType := encodeBasePairType(vc.sequence[n5p_iterator], vc.sequence[n5p_iterator_close_idx], vc.params.basePairEncodedTypeMap)

		fivePrimeMismatch = vc.encodedSequence[n5p_iterator-1]
		threePrimeMismatch = vc.encodedSequence[n5p_iterator_close_idx+1]

		energy += E_MLstem(basePairType, fivePrimeMismatch, threePrimeMismatch, vc.params)

		/* seek to the next stem */
		n5p_iterator = n5p_iterator_close_idx + 1

		for n5p_iterator < pairCloseIdx && pairTable[n5p_iterator] == -1 {
			n5p_iterator++
		}
		num_unpaired_nucs += n5p_iterator - n5p_iterator_close_idx - 1 /* add unpaired nucleotides */
	}

	if pairOpenIdx > 0 {
		/* actual closing pair */
		basePairType := encodeBasePairType(vc.sequence[pairCloseIdx], vc.sequence[pairOpenIdx], vc.params.basePairEncodedTypeMap)
		fivePrimeMismatch = vc.encodedSequence[pairCloseIdx-1]
		threePrimeMismatch = vc.encodedSequence[pairOpenIdx+1]

		energy += E_MLstem(basePairType, fivePrimeMismatch, threePrimeMismatch, vc.params)
	} else {
		/* virtual closing pair */
		energy += E_MLstem(0, -1, -1, vc.params)
	}

	energy += vc.params.MLclosing

	energy += num_unpaired_nucs * vc.params.MLbase

	return energy + bonus, nil
}

/**
 *  @def E_MLstem(A,B,C,D)
 *  <H2>Compute the Energy contribution of a Multiloop stem</H2>
 *  This definition is a wrapper for the E_Stem() funtion.
 *  It is substituted by an E_Stem() funtion call with argument
 *  extLoop=0, so the energy contribution returned reflects a
 *  stem introduced in a multiloop.<BR>
 *  As for the parameters B (si1) and C (sj1) of the substituted
 *  E_Stem() function, you can inhibit to take 5'-, 3'-dangles
 *  or mismatch contributions to be taken into account by passing
 *  -1 to these parameters.
 *
 *  @see    E_Stem()
 *  @param  A The pair type of the stem-closing pair
 *  @param  B The 5'-mismatching nucleotide
 *  @param  C The 3'-mismatching nucleotide
 *  @param  D The datastructure containing scaled energy parameters
 *  @return   The energy contribution of the introduced multiloop stem
 */
func E_MLstem(type_1 int, si1, sj1 int, P *energyParams) int {
	var energy int = 0

	if si1 >= 0 && sj1 >= 0 {
		energy += P.mismatchM[type_1][si1][sj1]
	} else if si1 >= 0 {
		energy += P.dangle5[type_1][si1]
	} else if sj1 >= 0 {
		energy += P.dangle3[type_1][sj1]
	}

	if type_1 > 2 {
		energy += P.TerminalAU
	}

	energy += P.MLintern[type_1]

	return energy
}

// Returns the minimum of two ints
func min(a, b int) int {
	if a < b {
		return a
	}
	return b
}

func get_scaled_params() *energyParams {
	var tempf int = int((defaultTemperature + K0) / Tmeasure)

	var params *energyParams = &energyParams{
		lxc:                    lxc37 * float64(tempf),
		TripleC:                rescaleDg(TripleC37, TripleCdH, tempf),
		MultipleCA:             rescaleDg(MultipleCA37, MultipleCAdH, tempf),
		MultipleCB:             rescaleDg(TerminalAU37, TerminalAUdH, tempf),
		TerminalAU:             rescaleDg(TerminalAU37, TerminalAUdH, tempf),
		MLbase:                 rescaleDg(ML_BASE37, ML_BASEdH, tempf),
		MLclosing:              rescaleDg(ML_closing37, ML_closingdH, tempf),
		basePairEncodedTypeMap: defaultBasePairEncodedTypeMap(),
	}

	params.ninio[2] = rescaleDg(ninio37, niniodH, tempf)

	for i := 0; i < 31; i++ {
		params.hairpin[i] = rescaleDg(hairpin37[i], hairpindH[i], tempf)
	}

	var i int
	for i = 0; i <= min(30, maxLoop); i++ {
		params.bulge[i] = rescaleDg(bulge37[i], bulgedH[i], tempf)
		params.internal_loop[i] = rescaleDg(internal_loop37[i], internal_loopdH[i], tempf)
	}

	for ; i <= maxLoop; i++ {
		params.bulge[i] = params.bulge[30] +
			int(params.lxc*math.Log(float64(i)/30.0))
		params.internal_loop[i] = params.internal_loop[30] +
			int(params.lxc*math.Log(float64(i)/30.0))
	}

	// This could be 281 or only the amount of chars
	for i := 0; (i * 7) < len(Tetraloops); i++ {
		params.Tetraloop_E[i] = rescaleDg(Tetraloop37[i], TetraloopdH[i], tempf)
	}

	for i := 0; (i * 5) < len(Triloops); i++ {
		params.Triloop_E[i] = rescaleDg(Triloop37[i], TriloopdH[i], tempf)
	}

	for i := 0; (i * 9) < len(Hexaloops); i++ {
		params.Hexaloop_E[i] = rescaleDg(Hexaloop37[i], HexaloopdH[i], tempf)
	}

	for i := 0; i <= nbPairs; i++ {
		params.MLintern[i] = rescaleDg(ML_intern37, ML_interndH, tempf)
	}

	/* stacks    G(T) = H - [H - G(T0)]*T/T0 */
	for i := 0; i <= nbPairs; i++ {
		for j := 0; j <= nbPairs; j++ {
			params.stack[i][j] = rescaleDg(stack37[i][j],
				stackdH[i][j],
				tempf)
		}
	}

	/* mismatches */
	for i := 0; i <= nbPairs; i++ {
		for j := 0; j < 5; j++ {
			for k := 0; k < 5; k++ {
				var mm int
				params.mismatchI[i][j][k] = rescaleDg(mismatchI37[i][j][k],
					mismatchIdH[i][j][k],
					tempf)
				params.mismatchH[i][j][k] = rescaleDg(mismatchH37[i][j][k],
					mismatchHdH[i][j][k],
					tempf)
				params.mismatch1nI[i][j][k] = rescaleDg(mismatch1nI37[i][j][k],
					mismatch1nIdH[i][j][k],
					tempf)
				params.mismatch23I[i][j][k] = rescaleDg(mismatch23I37[i][j][k],
					mismatch23IdH[i][j][k],
					tempf)
				// if model_details.dangles > 0 {
				mm = rescaleDg(mismatchM37[i][j][k],
					mismatchMdH[i][j][k],
					tempf)
				if mm > 0 {
					params.mismatchM[i][j][k] = 0
				} else {
					params.mismatchM[i][j][k] = mm
				}

				mm = rescaleDg(mismatchExt37[i][j][k],
					mismatchExtdH[i][j][k],
					tempf)
				if mm > 0 {
					params.mismatchExt[i][j][k] = 0
				} else {
					params.mismatchExt[i][j][k] = mm
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
					params.int11[i][j][k][l] = rescaleDg(int11_37[i][j][k][l],
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
						params.int21[i][j][k][l][m] = rescaleDg(int21_37[i][j][k][l][m],
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
							params.int22[i][j][k][l][m][n] = rescaleDg(int22_37[i][j][k][l][m][n],
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
