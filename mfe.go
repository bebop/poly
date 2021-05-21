package poly

import (
	"errors"
	"fmt"
	"log"
	"math"
	"strings"
	"unicode"
)

/*
* The MFE calculation functionality has been taken directly from ViennaRNA with the follwing changes:
* - ViennaRNA includes the ability to specify a dangle model (more info available at [src/bin/RNAeval.ggo#L153](https://github.com/ViennaRNA/ViennaRNA/blob/d6fbaf2b33d65944f9249a91ed5ab4b3277f7d06/src/bin/RNAeval.ggo#L153)). This implementation keeps it simple and defaults to the value of -d2.
* - ViennaRNA includes the ability to specify hard and soft constraints (more info available from [their ppt explaining hard and soft constrains](http://benasque.org/2015rna/talks_contr/2713_lorenz_benasque_2015.pdf), [official docs](https://www.tbi.univie.ac.at/RNA/ViennaRNA/doc/html/group__constraints.html), [thier paper](https://almob.biomedcentral.com/articles/10.1186/s13015-016-0070-z)). This implementation keeps it simple and defaults to no hard or soft constraints.
* - ViennaRNA includes the ability to calculate the minimum free energy of co-folded sequences. This implementation keeps it simple and defaults to calculating the mfe of only single sequences.
 */

/**
 *  @brief  The most basic data structure required by many functions throughout the RNAlib
 *
 */
type vrna_fold_compound_t struct {
	length int /**<  @brief  The length of the sequence (or sequence alignment). vivek: In ViennaRNA, length is stored as unsigned int (a primitive C type which ranges from 0 to 4,294,967,295). uint32 is the equivalent type in Go. */

	params            *vrna_param_t /**<  @brief  The precomputed free energy contributions for each type of loop */
	sequence          string        /* The input sequence string */
	sequence_encoding []int         /**<  @brief  Numerical encoding of the input sequence
	*    @see    vrna_sequence_encode()
	 */
}

func min(a, b int) int {
	if a < b {
		return a
	}
	return b
}

func CalculateMfe(seq, structure string) (float64, error) {
	if len(seq) != len(structure) {
		return 0, fmt.Errorf("length of sequence (%v) != Length of structure (%v)", len(seq), len(structure))
	} else if len(seq) == 0 {
		return 0, errors.New("lengths of sequence and structure cannot be 0")
	}

	fc := &vrna_fold_compound_t{
		length:            len(seq),
		params:            vrna_params(),
		sequence:          seq,
		sequence_encoding: encodeSequence(seq),
	}

	energy, err := vrna_eval_structure_cstr(fc, structure)
	if err != nil {
		return 0, err
	}
	return energy, nil
}

func vrna_params() *vrna_param_t {
	return get_scaled_params()
}

func rescaleDg(dG, dH, dT int) int {
	return (dH - (dH-dG)*dT)
}

var (
	DEFAULT_TEMP float64 = 37.0
	NB_PAIRS     int     = 7 /** The number of distinguishable base pairs */
)

const (
	MAXALPHA int = 20 /** @brief Maximal length of alphabet */
	MAXLOOP  int = 30 /** The maximum loop length */
)

type vrna_md_t struct {
	temperature float64                         /**<  @brief  The temperature used to scale the thermodynamic parameters */
	pair        [MAXALPHA + 1][MAXALPHA + 1]int /**<  @brief  Integer representation of a base pair */
}

var (
	// Default temperature for structure prediction and free energy evaluation in $^\circ C$
	VRNA_MODEL_DEFAULT_TEMPERATURE float64 = 37.0
)

const NBPAIRS = 7 /** The number of distinguishable base pairs */

func DefaultModelDetails() *vrna_md_t {
	return &vrna_md_t{
		temperature: VRNA_MODEL_DEFAULT_TEMPERATURE,
		// vivek: Why isn't this a symmetric matrix?
		// I couldn't find anything about this field in their docs
		pair: [MAXALPHA + 1][MAXALPHA + 1]int{
			//			 _  A  C  G  U  X  K  I
			/* _ */ {0, 0, 0, 0, 0, 0, 0, 0},
			/* A */ {0, 0, 0, 0, 5, 0, 0, 5},
			/* C */ {0, 0, 0, 1, 0, 0, 0, 0},
			/* G */ {0, 0, 2, 0, 3, 0, 0, 0},
			/* U */ {0, 6, 0, 4, 0, 0, 0, 6},
			/* X */ {0, 0, 0, 0, 0, 0, 2, 0},
			/* K */ {0, 0, 0, 0, 0, 1, 0, 0},
			/* I */ {0, 6, 0, 0, 5, 0, 0, 0},
		},
	}
}

type vrna_param_t struct {
	stack                                     [NBPAIRS + 1][NBPAIRS + 1]int
	hairpin                                   [31]int
	bulge                                     [MAXLOOP + 1]int
	internal_loop                             [MAXLOOP + 1]int
	mismatchExt                               [NBPAIRS + 1][5][5]int
	mismatchI                                 [NBPAIRS + 1][5][5]int
	mismatch1nI                               [NBPAIRS + 1][5][5]int
	mismatch23I                               [NBPAIRS + 1][5][5]int
	mismatchH                                 [NBPAIRS + 1][5][5]int
	mismatchM                                 [NBPAIRS + 1][5][5]int
	dangle5                                   [NBPAIRS + 1][5]int
	dangle3                                   [NBPAIRS + 1][5]int
	int11                                     [NBPAIRS + 1][NBPAIRS + 1][5][5]int
	int21                                     [NBPAIRS + 1][NBPAIRS + 1][5][5][5]int
	int22                                     [NBPAIRS + 1][NBPAIRS + 1][5][5][5][5]int
	ninio                                     [5]int
	lxc                                       float64
	MLbase, MLclosing, TerminalAU, DuplexInit int
	MLintern                                  [NBPAIRS + 1]int
	Tetraloop_E                               [200]int
	Tetraloops                                string
	Triloop_E                                 [40]int
	Triloops                                  string
	Hexaloop_E                                [40]int
	Hexaloops                                 string
	TripleC, MultipleCA, MultipleCB           int
	temperature                               float64    /**<  @brief  Temperature used for loop contribution scaling */
	model_details                             *vrna_md_t /**<  @brief  Model details to be used in the recursions */
}

func get_scaled_params() *vrna_param_t {
	var model_details *vrna_md_t = DefaultModelDetails()
	// float64
	var i, j, k, l int
	var tempf float64 = (model_details.temperature + K0) / Tmeasure

	var params *vrna_param_t = &vrna_param_t{
		model_details: model_details,
		temperature:   model_details.temperature,
		lxc:           lxc37 * tempf,
		TripleC:       rescaleDg(TripleC37, TripleCdH, int(tempf)),
		MultipleCA:    rescaleDg(MultipleCA37, MultipleCAdH, int(tempf)),
		MultipleCB:    rescaleDg(TerminalAU37, TerminalAUdH, int(tempf)),
		TerminalAU:    rescaleDg(TerminalAU37, TerminalAUdH, int(tempf)),
		// DuplexInit:            rescaleDg(DuplexInit37, DuplexInitdH, int(tempf)),
		MLbase:    rescaleDg(ML_BASE37, ML_BASEdH, int(tempf)),
		MLclosing: rescaleDg(ML_closing37, ML_closingdH, int(tempf)),
	}

	params.ninio[2] = rescaleDg(ninio37, niniodH, int(tempf))

	for i = 0; i < 31; i++ {
		params.hairpin[i] = rescaleDg(hairpin37[i], hairpindH[i], int(tempf))
	}

	for i = 0; i <= min(30, MAXLOOP); i++ {
		params.bulge[i] = rescaleDg(bulge37[i], bulgedH[i], int(tempf))
		params.internal_loop[i] = rescaleDg(internal_loop37[i], internal_loopdH[i], int(tempf))
	}

	for ; i <= MAXLOOP; i++ {
		params.bulge[i] = params.bulge[30] +
			int(params.lxc*math.Log(float64(i)/30.0))
		params.internal_loop[i] = params.internal_loop[30] +
			int(params.lxc*math.Log(float64(i)/30.0))
	}

	// This could be 281 or only the amount of chars
	for i = 0; (i * 7) < len(Tetraloops); i++ {
		params.Tetraloop_E[i] = rescaleDg(Tetraloop37[i], TetraloopdH[i], int(tempf))
	}

	for i = 0; (i * 5) < len(Triloops); i++ {
		params.Triloop_E[i] = rescaleDg(Triloop37[i], TriloopdH[i], int(tempf))
	}

	for i = 0; (i * 9) < len(Hexaloops); i++ {
		params.Hexaloop_E[i] = rescaleDg(Hexaloop37[i], HexaloopdH[i], int(tempf))
	}

	for i = 0; i <= NBPAIRS; i++ {
		params.MLintern[i] = rescaleDg(ML_intern37, ML_interndH, int(tempf))
	}

	/* stacks    G(T) = H - [H - G(T0)]*T/T0 */
	for i = 0; i <= NBPAIRS; i++ {
		for j = 0; j <= NBPAIRS; j++ {
			params.stack[i][j] = rescaleDg(stack37[i][j],
				stackdH[i][j],
				int(tempf))
		}
	}

	/* mismatches */
	for i = 0; i <= NBPAIRS; i++ {
		for j = 0; j < 5; j++ {
			for k = 0; k < 5; k++ {
				var mm int
				params.mismatchI[i][j][k] = rescaleDg(mismatchI37[i][j][k],
					mismatchIdH[i][j][k],
					int(tempf))
				params.mismatchH[i][j][k] = rescaleDg(mismatchH37[i][j][k],
					mismatchHdH[i][j][k],
					int(tempf))
				params.mismatch1nI[i][j][k] = rescaleDg(mismatch1nI37[i][j][k],
					mismatch1nIdH[i][j][k],
					int(tempf))
				params.mismatch23I[i][j][k] = rescaleDg(mismatch23I37[i][j][k],
					mismatch23IdH[i][j][k],
					int(tempf))
				// if model_details.dangles > 0 {
				mm = rescaleDg(mismatchM37[i][j][k],
					mismatchMdH[i][j][k],
					int(tempf))
				if mm > 0 {
					params.mismatchM[i][j][k] = 0
				} else {
					params.mismatchM[i][j][k] = mm
				}

				mm = rescaleDg(mismatchExt37[i][j][k],
					mismatchExtdH[i][j][k],
					int(tempf))
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
	for i = 0; i <= NBPAIRS; i++ {
		for j = 0; j < 5; j++ {
			var dd int
			dd = rescaleDg(dangle5_37[i][j],
				dangle5_dH[i][j],
				int(tempf))
			if dd > 0 {
				params.dangle5[i][j] = 0
			} else {
				params.dangle5[i][j] = dd
			}

			dd = rescaleDg(dangle3_37[i][j],
				dangle3_dH[i][j],
				int(tempf))
			if dd > 0 {
				params.dangle3[i][j] = 0
			} else {
				params.dangle3[i][j] = dd
			}
		}
	}

	/* interior 1x1 loops */
	for i = 0; i <= NBPAIRS; i++ {
		for j = 0; j <= NBPAIRS; j++ {
			for k = 0; k < 5; k++ {
				for l = 0; l < 5; l++ {
					params.int11[i][j][k][l] = rescaleDg(int11_37[i][j][k][l],
						int11_dH[i][j][k][l],
						int(tempf))
				}
			}
		}
	}

	/* interior 2x1 loops */
	for i = 0; i <= NBPAIRS; i++ {
		for j = 0; j <= NBPAIRS; j++ {
			for k = 0; k < 5; k++ {
				for l = 0; l < 5; l++ {
					var m int
					for m = 0; m < 5; m++ {
						params.int21[i][j][k][l][m] = rescaleDg(int21_37[i][j][k][l][m],
							int21_dH[i][j][k][l][m],
							int(tempf))
					}
				}
			}
		}
	}

	/* interior 2x2 loops */
	for i = 0; i <= NBPAIRS; i++ {
		for j = 0; j <= NBPAIRS; j++ {
			for k = 0; k < 5; k++ {
				for l = 0; l < 5; l++ {
					var m, n int
					for m = 0; m < 5; m++ {
						for n = 0; n < 5; n++ {
							params.int22[i][j][k][l][m][n] = rescaleDg(int22_37[i][j][k][l][m][n],
								int22_dH[i][j][k][l][m][n],
								int(tempf))
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

func vrna_eval_structure_cstr(fc *vrna_fold_compound_t, structure string) (float64, error) {
	// var pt []int
	pair_table, err := vrna_pair_table_from_string(structure)
	if err != nil {
		return 0, err
	}
	// log.Printf("vrna_eval_structure_cstr pair_table: %v", pair_table)

	// var en float64
	energy := float64(eval_pair_table(fc, pair_table)) / 100.0
	// en, err := wrap_eval_structure(fc, structure, pair_table)
	// if err != nil {
	// 	return 0, err
	// }

	return energy, nil
}

var SHRT_MAX int = 32767

/* vivek:
* returns a slice pair_table where pair_table[sequence[i]] returns the index of the
* of the paired opening/closing bracket of i if i is an opening/closing bracket,
* else 0.
* Note: If indexes in pair_table are 1-indexed as 0 is used as the default value
* Examples -
* Input:    ". .  (  (  (  ( . . . ) ) ) ) . . .  (  ( . . . . . . . .  )  ) . ."
* Output:[30 0 0 13 12 11 10 0 0 0 6 5 4 3 0 0 0 28 27 0 0 0 0 0 0 0 0 18 17 0 0 0]
*
* Input:    "(  .  (  (  (  ( . . . ) ) ) ) . . .  (  ( . . . . . . . .  )  ) . )"
* Output:[30 30 0 13 12 11 10 0 0 0 6 5 4 3 0 0 0 28 27 0 0 0 0 0 0 0 0 18 17 0 1 0]
 */
func vrna_pair_table_from_string(structure string) ([]int, error) {
	var pair_table []int
	len_seq := len(structure)

	if len_seq > SHRT_MAX {
		err := fmt.Sprintf("vrna_pair_table_from_string: Structure too long to be converted to pair table (n=%v, max=%v)", len_seq, SHRT_MAX)
		return nil, errors.New(err)
	}

	pair_table = make([]int, len_seq+2)
	// pair_table[0] = len_seq

	// a slice of the characters in the structure
	var structure_char_slice []rune = []rune(structure)

	// the characters of the opening and closing brackets
	var open_bracket, close_bracket rune = '(', ')'

	// keeps track of the indexes of open brackets. indexes of open brackets are pushed onto stack
	// and poped off when a closing bracket is encountered

	var open_bracket_idx_stack []int = make([]int, len_seq)
	var stack_idx int = 0

	// iterate through structure_char_slice and create pair_table
	for i := 0; i < len(structure_char_slice); i++ {
		curr_idx := i + 1
		if structure_char_slice[i] == open_bracket {
			// pop index of open bracket onto stack
			open_bracket_idx_stack[stack_idx] = curr_idx
			stack_idx++
		} else if structure_char_slice[i] == close_bracket {
			stack_idx--

			if stack_idx < 0 {
				return nil, fmt.Errorf("%v\nunbalanced brackets '%v%v' found while extracting base pairs", structure,
					open_bracket, close_bracket)
			}

			open_bracket_idx := open_bracket_idx_stack[stack_idx]
			// current index of one-indexed sequence
			pair_table[curr_idx] = open_bracket_idx
			pair_table[open_bracket_idx] = curr_idx
		} else {
			pair_table[curr_idx] = -1
		}
	}

	if stack_idx != 0 {
		return nil, fmt.Errorf("%v\nunbalanced brackets '%v%v' found while extracting base pairs", structure,
			open_bracket, close_bracket)
	}

	return pair_table, nil
}

/* vivek:
*
*
 */
func eval_pair_table(fc *vrna_fold_compound_t, pair_table []int) int {
	var i, length, energy int

	length = int(fc.length)

	energy = energy_of_extLoop_pair_table(fc, 0, pair_table)

	vrna_cstr_print_eval_ext_loop(energy)
	// vivek: assume input is only A, T, C, G, or U
	for i = 1; i <= length; i++ {
		if pair_table[i] == -1 {
			continue
		}

		energy += stack_energy(fc, i, pair_table)
		i = pair_table[i]
	}

	return energy
}

func vrna_cstr_print_eval_ext_loop(energy int) {
	log.Printf("External loop                           : %v\n", energy)
}

/* vivek:
* encodes `sequence` based on `encodeNucelotide` into `S[1:len(sequence)-1]`
* `S[0]` is the length of `sequence`
* `S[len(sequence) + 1]` (or `S[-1]`) == `S[1]` which makes the sequence circular
* not sure why the orig repo makes the sequence cirular
* Contains the numerical encoding of the pair type for each pair (i,j) used

* manipulates the encoded sequence of `sequence` (see doc of `vrna_seq_encode_simple` for more
* details of the encoding process) by setting `S[0]` to `S[len(sequence)]`
* If a sequence is axxx...xxxb where a and b are the first and last nucleotides
* of the sequence, and AXXX...XXXB is the encoded `sequence`, `encodeSequence`
* returns BAXXX...XXXBA.
* thinking out loud:
* * could be done to make the sequence circular
* * may be needed in functions that require the previous and next values to
* compute a value

* Encodes a sequence into its numerical representation (see `encodeNucelotide`
* for the rune to int mapping) into `S[1:len(sequence)-1]`.

 */
func encodeSequence(sequence string) []int {
	var l int = len(sequence)
	var S []int = make([]int, l)

	for i := 0; i < l; i++ { /* make numerical encoding of sequence */
		S[i] = encodeNucelotide(([]rune(sequence))[i])
	}

	// S[0] = S[l]

	return S
}

// vivek: I don't know if this variable name is a joke or not...
// var Law_and_Order []rune = []rune("_ACGUTXKI")
var Law_and_Order string = "_ACGUTXKI"

/* vivek:
* Maps:
*	A -> 1
* C -> 2
* G -> 3
* T / U -> 4
 */
func encodeNucelotide(c rune) int {
	/* return numerical representation of nucleotide used e.g. in vrna_md_t.pair[][] */ //
	var code int = -1

	c = unicode.ToUpper(c)

	pos := strings.Index(Law_and_Order, string(c))
	if pos == -1 {
		code = 0
	} else {
		code = pos
	}

	if code > 5 {
		code = 0
	}

	if code > 4 {
		code-- /* make T and U equivalent */
	}

	return code
}

/*
* Calculate the energy contribution of
* stabilizing dangling-ends/mismatches
* for all stems branching off the exterior
* loop
* vivek: What is `i`?
 */
func energy_of_extLoop_pair_table(fc *vrna_fold_compound_t, i int, pair_table []int) int {
	var energy int = 0
	length := int(fc.length)
	pair_table_idx := 1

	/* seek to opening base of first stem (pair) */
	for pair_table_idx <= length && (pair_table[pair_table_idx] == -1) {
		pair_table_idx++
	}

	for pair_table_idx < length {
		/* pair_table_idx  must have a pairing partner */
		pair_close_idx := int(pair_table[pair_table_idx])

		/* get type of base pair (p,q) */
		pair_type := vrna_get_pair_type_md(fc.sequence_encoding[pair_table_idx-1],
			fc.sequence_encoding[pair_close_idx-1],
			fc.params.model_details)

		// vivek: cryptic variable names. 5 and 3 have to do with the 5' and 3' ends
		// of a DNA/RNA strand
		var mm5, mm3 int
		if pair_table_idx > 1 {
			mm5 = fc.sequence_encoding[pair_table_idx-1-1]
		} else {
			mm5 = -1
		}

		if pair_close_idx < length {
			mm3 = fc.sequence_encoding[pair_close_idx]
		} else {
			mm3 = -1
		}

		energy += vrna_E_ext_stem(pair_type, mm5, mm3, fc.params)

		/* seek to the next stem */
		pair_table_idx = pair_close_idx + 1
		for pair_table_idx <= length && pair_table[pair_table_idx] == -1 {
			pair_table_idx++
		}

		// vivek: is this a bug? should it be >= instead?
		// Is it likely that pair_table_idx == i given we seek to next stem above?
		if pair_table_idx == i {
			break /* cut was in loop */
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
 *  If either of the adjacent nucleotides i - 1 (`n5d`) and j + 1 (`n3d`) must not
 *  contribute stacking energy, the corresponding encoding must be -1.
 *
 *  @see vrna_E_exp_stem()
 *
 *  @param  type  The base pair encoding
 *  @param  n5d   The encoded nucleotide directly adjacent at the 5' side of the base pair (may be -1)
 *  @param  n3d   The encoded nucleotide directly adjacent at the 3' side of the base pair (may be -1)
 *  @param  p     The pre-computed energy parameters
 *  @return       The energy contribution of the introduced exterior-loop stem
 *
 * vivek: `n5d` and `n3d` are always positive unless the opening bracket is the
 * first nucelotide of the sequence (in which case n5d < 0), or the closing bracket
 * is the last nucelotide of the sequence (in which case n3d < 0).
 */
func vrna_E_ext_stem(pair_type int, n5d int, n3d int, p *vrna_param_t) int {
	var energy int = 0

	if n5d >= 0 && n3d >= 0 {
		energy += p.mismatchExt[pair_type][n5d][n3d]
	} else if n5d >= 0 {
		// only occurs in the case opening bracket of pair is first nucelotide of sequence
		energy += p.dangle5[pair_type][n5d]
	} else if n3d >= 0 {
		// only occurs in the case closing bracket of pair is last nucelotide of sequence
		energy += p.dangle3[pair_type][n3d]
	}

	if pair_type > 2 {
		energy += p.TerminalAU
	}

	return energy
}

/* vivek:
*	what is i? what is stack_energy?
 */
func stack_energy(fc *vrna_fold_compound_t, pair_open_idx int, pair_table []int) int {
	/* recursively calculate energy of substructure enclosed by (pair_open_idx, pair_close_idx) */
	energy := 0
	pair_close_idx := pair_table[pair_open_idx]
	sequence := fc.sequence

	if fc.params.model_details.pair[fc.sequence_encoding[pair_open_idx-1]][fc.sequence_encoding[pair_close_idx-1]] == 0 {
		panic(fmt.Sprintf("bases %v and %v (%v%v) can't pair!",
			pair_open_idx, pair_close_idx,
			string(sequence[pair_open_idx-1]),
			string(sequence[pair_close_idx-1])))
	}

	n5p_iterator := pair_open_idx
	n3p_iterator := pair_close_idx

	for n5p_iterator < n3p_iterator {
		/* process all stacks and interior loops */

		// seek to opening bracket from 5' end
		n5p_iterator++
		for pair_table[n5p_iterator] == -1 {
			n5p_iterator++
		}

		// seek to closing bracket from 3' end
		n3p_iterator--
		for pair_table[n3p_iterator] == -1 {
			n3p_iterator--
		}

		if (pair_table[n3p_iterator] != n5p_iterator) || (n5p_iterator > n3p_iterator) {
			break
		}

		// vivek: should this be pair[n5p_iterator][n3p_iterator] or is the current order correct?
		if fc.params.model_details.pair[fc.sequence_encoding[n3p_iterator-1]][fc.sequence_encoding[n5p_iterator-1]] == 0 {
			panic(fmt.Sprintf("bases %d and %d (%c%c) can't pair!",
				n5p_iterator, n3p_iterator,
				sequence[n5p_iterator-1],
				sequence[n3p_iterator-1]))
		}

		ee := eval_int_loop(fc, pair_open_idx, pair_close_idx, n5p_iterator, n3p_iterator)
		energy += ee

		vrna_cstr_print_eval_int_loop(pair_open_idx, pair_close_idx,
			sequence[pair_open_idx-1], sequence[pair_close_idx-1],
			n5p_iterator, n3p_iterator,
			sequence[n5p_iterator-1], sequence[n3p_iterator-1], ee)

		pair_open_idx = n5p_iterator
		pair_close_idx = n3p_iterator
	} /* end while */

	/* n5p_iterator, n3p_iterator don't pair must have found hairpin or multiloop */

	if n5p_iterator > n3p_iterator {
		/* hairpin */

		ee := vrna_eval_hp_loop(fc, pair_open_idx, pair_close_idx)
		energy += ee

		vrna_cstr_print_eval_hp_loop(pair_open_idx, pair_close_idx, sequence[pair_open_idx-1], sequence[pair_close_idx-1], ee)

		return energy
	}

	/* (pair_open_idx, pair_close_idx) is exterior pair of multiloop */
	for n5p_iterator < pair_close_idx {
		/* add up the contributions of the substructures of the ML */
		energy += stack_energy(fc, n5p_iterator, pair_table)
		n5p_iterator = pair_table[n5p_iterator]
		/* search for next base pair in multiloop */
		n5p_iterator++
		for pair_table[n5p_iterator] == -1 {
			n5p_iterator++
		}
	}

	var ee int = 0
	var err error

	ee, err = energy_of_ml_pair_table(fc, pair_open_idx, pair_table)
	if err != nil {
		panic(err)
	}

	energy += ee

	vrna_cstr_print_eval_mb_loop(pair_open_idx, pair_close_idx, sequence[pair_open_idx-1], sequence[pair_close_idx-1], ee)

	return energy
}

func vrna_cstr_print_eval_int_loop(i, j int, si, sj byte, k, l int, sk, sl byte, energy int) {
	log.Printf("Interior loop (%v,%v) %v%v; (%v,%v) %v%v: %v\n",
		i, j,
		string(si), string(sj),
		k, l,
		string(sk), string(sl),
		energy)
}

func vrna_cstr_print_eval_hp_loop(i, j int, si, sj byte, energy int) {
	log.Printf("Hairpin loop  (%v,%v) %v%v              : %v\n",
		i, j,
		string(si), string(sj),
		energy)
}

func vrna_cstr_print_eval_mb_loop(i, j int, si, sj byte, energy int) {
	log.Printf("Multi   loop (%v,%v) %v%v              : %v\n",
		i, j,
		string(si), string(sj),
		energy)
}

/**
 *  @brief Evaluate the free energy contribution of an interior loop with delimiting
 *  base pairs (i,j) and (k,l)
 *
 */
func eval_int_loop(fc *vrna_fold_compound_t, i, j, k, l int) int {
	var u1, u2 int
	energy := 0

	u1 = k - i - 1
	u2 = j - l - 1

	ij_pair_type := vrna_get_pair_type_md(fc.sequence_encoding[i-1], fc.sequence_encoding[j-1], fc.params.model_details)
	kl_pair_type := vrna_get_pair_type_md(fc.sequence_encoding[l-1], fc.sequence_encoding[k-1], fc.params.model_details)

	energy = E_IntLoop(u1, u2, ij_pair_type, kl_pair_type, fc.sequence_encoding[i], fc.sequence_encoding[j-1-1], fc.sequence_encoding[k-1-1], fc.sequence_encoding[l], fc.params)

	return energy
}

/* vivek:
*	Couldn't find any documentation about what `md.pair` is.
* Created an issue to get more info: https://github.com/ViennaRNA/ViennaRNA/issues/124
* retruns type of base pair (p,q)
 */
func vrna_get_pair_type_md(i, j int, md *vrna_md_t) int {
	var tt int = md.pair[i][j]

	if tt == 0 {
		return 7
	} else {
		return tt
	}
}

/**
 *  <H2>Compute the Energy of an interior-loop</H2>
 *  This function computes the free energy @f$\Delta G@f$ of an interior-loop with the
 *  following structure: <BR>
 *  <PRE>
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
 *  </PRE>
 *  This general structure depicts an interior-loop that is closed by the base pair (X,Y).
 *  The enclosed base pair is (V,U) which leaves the unpaired bases a_1-a_n and b_1-b_n
 *  that constitute the loop. In this example, the length of the interior-loop is @f$(n+m)@f$
 *  where n or m may be 0 resulting in a bulge-loop or base pair stack.
 *  The mismatching nucleotides for the closing pair (X,Y) are:
 *  5'-mismatch: a_1
 *  3'-mismatch: b_m
 *  and for the enclosed base pair (V,U):
 *  5'-mismatch: b_1
 *  3'-mismatch: a_n
 *  @note Base pairs are always denoted in 5'->3' direction. Thus the enclosed base pair
 *  must be 'turned arround' when evaluating the free energy of the interior-loop
 *
 *  @param  n1      The size of the 'left'-loop (number of unpaired nucleotides)
 *  @param  n2      The size of the 'right'-loop (number of unpaired nucleotides)
 *  @param  type    The pair type of the base pair closing the interior loop
 *  @param  type_2  The pair type of the enclosed base pair
 *  @param  si1     The 5'-mismatching nucleotide of the closing pair
 *  @param  sj1     The 3'-mismatching nucleotide of the closing pair
 *  @param  sp1     The 3'-mismatching nucleotide of the enclosed pair
 *  @param  sq1     The 5'-mismatching nucleotide of the enclosed pair
 *  @param  P       The datastructure containing scaled energy parameters
 *  @return The Free energy of the Interior-loop in dcal/mol
 */
func E_IntLoop(n1, n2 int, type_1, type_2 int, si1, sj1, sp1, sq1 int, P *vrna_param_t) int {
	/* compute energy of degree 2 loop (stack bulge or interior) */
	var nl, ns, u, energy int

	if n1 > n2 {
		nl = n1
		ns = n2
	} else {
		nl = n2
		ns = n1
	}

	if nl == 0 {
		return P.stack[type_1][type_2] /* stack */
	}

	if ns == 0 {
		/* bulge */
		if nl <= MAXLOOP {
			energy = P.bulge[nl]
		} else {
			energy = P.bulge[30] + int(P.lxc*math.Log(float64(nl)/30.0))
		}

		if nl == 1 {
			energy += P.stack[type_1][type_2]
		} else {
			if type_1 > 2 {
				energy += P.TerminalAU
			}

			if type_2 > 2 {
				energy += P.TerminalAU
			}
		}

		return energy
	} else {
		/* interior loop */
		if ns == 1 {
			if nl == 1 {
				/* 1x1 loop */
				return P.int11[type_1][type_2][si1][sj1]
			}

			if nl == 2 {
				/* 2x1 loop */
				if n1 == 1 {
					energy = P.int21[type_1][type_2][si1][sq1][sj1]
				} else {
					energy = P.int21[type_2][type_1][sq1][si1][sp1]
				}
				return energy
			} else {
				/* 1xn loop */
				if nl+1 <= MAXLOOP {
					energy = P.internal_loop[nl+1]
				} else {
					energy = P.internal_loop[30] + int(P.lxc*math.Log((float64(nl)+1.0)/30.0))
				}
				energy += min(MAX_NINIO, (nl-ns)*P.ninio[2])
				energy += P.mismatch1nI[type_1][si1][sj1] + P.mismatch1nI[type_2][sq1][sp1]
				return energy
			}
		} else if ns == 2 {
			if nl == 2 {
				/* 2x2 loop */
				return P.int22[type_1][type_2][si1][sp1][sq1][sj1]
			} else if nl == 3 {
				/* 2x3 loop */
				energy = P.internal_loop[5] + P.ninio[2]
				energy += P.mismatch23I[type_1][si1][sj1] + P.mismatch23I[type_2][sq1][sp1]
				return energy
			}
		}

		{
			/* generic interior loop (no else here!)*/
			u = nl + ns
			if u <= MAXLOOP {
				energy = P.internal_loop[u]
			} else {
				energy = P.internal_loop[30] + int(P.lxc*math.Log(float64(u)/30.0))
			}

			energy += min(MAX_NINIO, (nl-ns)*P.ninio[2])

			energy += P.mismatchI[type_1][si1][sj1] + P.mismatchI[type_2][sq1][sp1]
		}
	}

	return energy
}

/**
 *  @brief Evaluate free energy of a hairpin loop
 *
 *  @ingroup eval
 *
 *  @note This function is polymorphic! The provided #vrna_fold_compound_t may be of type
 *  #VRNA_FC_TYPE_SINGLE or #VRNA_FC_TYPE_COMPARATIVE
 *
 *  @param  fc  The #vrna_fold_compound_t for the particular energy evaluation
 *  @param  i   5'-position of the base pair
 *  @param  j   3'-position of the base pair
 *  @returns    Free energy of the hairpin loop closed by @f$ (i,j) @f$ in deka-kal/mol
 */
func vrna_eval_hp_loop(fc *vrna_fold_compound_t, pair_open_idx, pair_close_idx int) int {
	var u, e int
	var P *vrna_param_t
	var md *vrna_md_t

	P = fc.params
	md = P.model_details

	u = pair_close_idx - pair_open_idx - 1
	pair_type := vrna_get_pair_type_md(fc.sequence_encoding[pair_open_idx-1], fc.sequence_encoding[pair_close_idx-1], md)

	e = E_Hairpin(u, pair_type, fc.sequence_encoding[pair_open_idx], fc.sequence_encoding[pair_close_idx-1-1], fc.sequence[pair_open_idx-1:], P)
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
 *  @see vrna_param_t
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
func E_Hairpin(size int, type_1 int, si1, sj1 int, sequence string, P *vrna_param_t) int {
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
 *** pair_open_idx is the 5'-base of the closing pair
 ***
 *** since each helix can coaxially stack with at most one of its
 *** neighbors we need an auxiliarry variable  cx_energy
 *** which contains the best energy given that the last two pairs stack.
 *** energy  holds the best energy given the previous two pairs do not
 *** stack (i.e. the two current helices may stack)
 *** We don't allow the last helix to stack with the first, thus we have to
 *** walk around the Loop twice with two starting points and take the minimum
 ***/
func energy_of_ml_pair_table(vc *vrna_fold_compound_t, pair_open_idx int, pair_table []int) (int, error) {

	var energy int
	var pair_close_idx, n5p_iterator_close_idx, num_unpaired_nucs, mm5, mm3 int

	length := int(vc.length)

	bonus := 0

	if pair_open_idx >= pair_table[pair_open_idx] {
		return INF, errors.New("energy_of_ml_pair_table: pair_open_idx is not 5' base of a closing pair")
	}

	if pair_open_idx == 0 {
		pair_close_idx = length + 1
	} else {
		pair_close_idx = pair_table[pair_open_idx]
	}

	/* init the variables */
	energy = 0
	num_unpaired_nucs = 0 /* the total number of unpaired nucleotides */
	n5p_iterator := pair_open_idx + 1

	for n5p_iterator <= pair_close_idx && pair_table[n5p_iterator] == -1 {
		n5p_iterator++
	}

	/* add bonus energies for first stretch of unpaired nucleotides */
	num_unpaired_nucs += n5p_iterator - pair_open_idx - 1

	for n5p_iterator < pair_close_idx {
		/* p must have a pairing partner */
		n5p_iterator_close_idx = pair_table[n5p_iterator]
		/* get type of base pair (p,q) */
		pair_type := vrna_get_pair_type_md(vc.sequence_encoding[n5p_iterator-1], vc.sequence_encoding[n5p_iterator_close_idx-1], vc.params.model_details)

		mm5 = vc.sequence_encoding[n5p_iterator-1-1]
		mm3 = vc.sequence_encoding[n5p_iterator_close_idx]

		energy += E_MLstem(pair_type, mm5, mm3, vc.params)

		/* seek to the next stem */
		n5p_iterator = n5p_iterator_close_idx + 1

		for n5p_iterator < pair_close_idx && pair_table[n5p_iterator] == -1 {
			n5p_iterator++
		}
		num_unpaired_nucs += n5p_iterator - n5p_iterator_close_idx - 1 /* add unpaired nucleotides */
	}

	if pair_open_idx > 0 {
		/* actual closing pair */
		pair_type := vrna_get_pair_type_md(vc.sequence_encoding[pair_close_idx-1], vc.sequence_encoding[pair_open_idx-1], vc.params.model_details)
		mm5 = vc.sequence_encoding[pair_close_idx-1-1]
		mm3 = vc.sequence_encoding[pair_open_idx]

		energy += E_MLstem(pair_type, mm5, mm3, vc.params)
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
func E_MLstem(type_1 int, si1, sj1 int, P *vrna_param_t) int {
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
