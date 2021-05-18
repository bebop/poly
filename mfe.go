package poly

import (
	"errors"
	"fmt"
	"math"
	"strings"
	"unicode"
)

/* Design decisions behind ViennaRNA (The software package from which the functions in this file have been derived from):
* About hard constraints (the following points have been taken directly from the paper: https://almob.biomedcentral.com/articles/10.1186/s13015-016-0070-z):
* Thermodynamics-based pseudo-knot free RNA secondary structure prediction is by no means perfect.
* The same is true of SCFG-based (stochastic context-free grammars) approaches.
* It is therefore of key interest to guide the RNA secondary structure predictions by incorporating experimental evidence beyond the parameters of the standard energy model.
* This view was emphasized already by proposing a constraint programming framework for RNA folding.
( Thermodynamic folding software (ViennaRNA Package) includes the possibility to constrain the set of allowed base pairs or to force individual nucleotides to be paired.
* Early approaches towards incorporating additional information, e.g. from chemical or enzymatic probing data or known chemical modification, used large energy penalties to effectively prohibit certain conformations.
* We refer to such all-or-none decisions as hard constraints.
*/

/**
 *  @brief  The most basic data structure required by many functions throughout the RNAlib
 *
 *  @note   Please read the documentation of this data structure carefully! Some attributes are only available for
 *  specific types this data structure can adopt.
 *
 *  @warning  Reading/Writing from/to attributes that are not within the scope of the current type usually result
 *  in undefined behavior!
 *
 *  @see  #vrna_fold_compound_t.type, vrna_fold_compound(), vrna_fold_compound_comparative(), vrna_fold_compound_free(),
 *        #VRNA_FC_TYPE_SINGLE, #VRNA_FC_TYPE_COMPARATIVE
 */
type vrna_fold_compound_t struct {
	length   uint32 /**<  @brief  The length of the sequence (or sequence alignment). vivek: In ViennaRNA, length is stored as unsigned int (a primitive C type which ranges from 0 to 4,294,967,295). uint32 is the equivalent type in Go. */
	cutpoint int    /*  @brief  The position of the (cofold) cutpoint within the provided sequence. If there is no cutpoint, this field will be set to -1 */

	strand_number []uint32 /**<  @brief  The strand number a particular nucleotide is associated with */
	strand_order  []uint32 /**<  @brief  The strand order, i.e. permutation of current concatenated sequence */
	strand_start  []uint32 /**<  @brief  The start position of a particular strand within the current concatenated sequence */
	strand_end    []uint32 /**<  @brief  The end (last) position of a particular strand within the current concatenated sequence */

	strands     uint32       /* vivek: The number of strands in this folding compound. Since we only calculate MFE for one compound, this is generally one.*/
	nucleotides []vrna_seq_t /* */
	// alignment   *vrna_msa_t

	// hc *vrna_hc_t /**<  @brief  The hard constraints data structure used for structure prediction */

	// matrices *vrna_mx_mfe_t /**<  @brief  The MFE DP matrices */
	// exp_matrices *vrna_mx_pf_t  /**<  @brief  The PF DP matrices  */

	params *vrna_param_t /**<  @brief  The precomputed free energy contributions for each type of loop */
	// exp_params *vrna_exp_param_t    /**<  @brief  The precomputed free energy contributions as Boltzmann factors  */

	iindx []int /**<  @brief  DP matrix accessor  */
	jindx []int /**<  @brief  DP matrix accessor  */

	/**
	 *  @}
	 *
	 *  @name Secondary Structure Decomposition (grammar) related data fields
	 *  @{
	 */

	/* data structure to adjust additional structural domains, such as G-quadruplexes */
	// domains_struc *vrna_sd_t             /**<  @brief  Additional structured domains */

	/* data structure to adjust additional contributions to unpaired stretches, e.g. due to protein binding */
	// domains_up     *vrna_ud_                /**<  @brief  Additional unstructured domains */

	/**
	 *  @}
	 */

	/**
	 *  @name Data fields available for single/hybrid structure prediction
	 *  @{
	 */
	sequence string /**<  @brief  The input sequence string
	 *    @warning   Only available if @verbatim type==VRNA_FC_TYPE_SINGLE @endverbatim
	 */
	sequence_encoding []int /**<  @brief  Numerical encoding of the input sequence
	 *    @see    vrna_sequence_encode()
	 *    @warning   Only available if @verbatim type==VRNA_FC_TYPE_SINGLE @endverbatim
	 */
	sequence_encoding2 []int
	ptype              string /**<  @brief  Pair type array
	 *
	 *    Contains the numerical encoding of the pair type for each pair (i,j) used
	 *    in MFE, Partition function and Evaluation computations.
	 *    @note This array is always indexed via jindx, in contrast to previously
	 *    different indexing between mfe and pf variants!
	 *    @warning   Only available if @verbatim type==VRNA_FC_TYPE_SINGLE @endverbatim
	 *    @see    vrna_idx_col_wise(), vrna_ptypes()
	 */

	sc *vrna_sc_t /**<  @brief  The soft constraints for usage in structure prediction and evaluation
	 *    @warning   Only available if @verbatim type==VRNA_FC_TYPE_SINGLE @endverbatim
	 */

	/**
	 *  @}
	 */

	// /**
	//  *  @name Data fields for consensus structure prediction
	//  *  @{
	//  */
	//     sequences []string        /**<  @brief  The aligned sequences
	//                                        *    @note   The end of the alignment is indicated by a NULL pointer in the second dimension
	//                                        *    @warning   Only available if @verbatim type==VRNA_FC_TYPE_COMPARATIVE @endverbatim
	//                                        */
	//     n_seq uint32              /**<  @brief  The number of sequences in the alignment
	//                                        *    @warning   Only available if @verbatim type==VRNA_FC_TYPE_COMPARATIVE @endverbatim
	//                                        */
	//     cons_seq string          /**<  @brief  The consensus sequence of the aligned sequences
	//                                        *    @warning   Only available if @verbatim type==VRNA_FC_TYPE_COMPARATIVE @endverbatim
	//                                        */
	//     S_cons []int           /**<  @brief  Numerical encoding of the consensus sequence
	//                                        *    @warning   Only available if @verbatim type==VRNA_FC_TYPE_COMPARATIVE @endverbatim
	//                                        */
	//     S   [][]int            /**<  @brief  Numerical encoding of the sequences in the alignment
	//                                        *    @warning   Only available if @verbatim type==VRNA_FC_TYPE_COMPARATIVE @endverbatim
	//                                        */
	//     S5  [][]int             /**<  @brief    S5[s][i] holds next base 5' of i in sequence s
	//                                        *    @warning  Only available if @verbatim type==VRNA_FC_TYPE_COMPARATIVE @endverbatim
	//                                        */
	//     S3  [][]int             /**<  @brief    Sl[s][i] holds next base 3' of i in sequence s
	//                                        *    @warning  Only available if @verbatim type==VRNA_FC_TYPE_COMPARATIVE @endverbatim
	//                                        */
	// 		Ss []string
	// 		a2s [][]uint32
	//     pscore []int              /**<  @brief  Precomputed array of pair types expressed as pairing scores
	//                                          *    @warning   Only available if @verbatim type==VRNA_FC_TYPE_COMPARATIVE @endverbatim
	//                                          */
	//     pscore_local [][]int       /**<  @brief  Precomputed array of pair types expressed as pairing scores
	//                                          *    @warning   Only available if @verbatim type==VRNA_FC_TYPE_COMPARATIVE @endverbatim
	//                                          */
	//     pscore_pf_compat []int    /**<  @brief  Precomputed array of pair types expressed as pairing scores indexed via iindx
	//                                          *    @deprecated  This attribute will vanish in the future!
	//                                          *    @warning   Only available if @verbatim type==VRNA_FC_TYPE_COMPARATIVE @endverbatim
	//                                          */
	//     scs []*vrna_sc_t                /**<  @brief  A set of soft constraints (for each sequence in the alignment)
	//                                          *    @warning   Only available if @verbatim type==VRNA_FC_TYPE_COMPARATIVE @endverbatim
	//                                          */
	// 		oldAliEn int

	/**
	 *  @}
	 */

	/**
	 *  @name Additional data fields for Distance Class Partitioning
	 *
	 *  These data fields are typically populated with meaningful data only if used in the context of Distance Class Partitioning
	 *  @{
	 */
	maxD1         uint32 /**<  @brief  Maximum allowed base pair distance to first reference */
	maxD2         uint32 /**<  @brief  Maximum allowed base pair distance to second reference */
	reference_pt1 []int  /**<  @brief  A pairtable of the first reference structure */
	reference_pt2 []int  /**<  @brief  A pairtable of the second reference structure */

	referenceBPs1 []uint32 /**<  @brief  Matrix containing number of basepairs of reference structure1 in interval [i,j] */
	referenceBPs2 []uint32 /**<  @brief  Matrix containing number of basepairs of reference structure2 in interval [i,j] */
	bpdist        []uint32 /**<  @brief  Matrix containing base pair distance of reference structure 1 and 2 on interval [i,j] */

	mm1 []uint32 /**<  @brief  Maximum matching matrix, reference struct 1 disallowed */
	mm2 []uint32 /**<  @brief  Maximum matching matrix, reference struct 2 disallowed */

	/**
	 *  @}
	 */

	/**
	 *  @name Additional data fields for local folding
	 *
	 *  These data fields are typically populated with meaningful data only if used in the context of local folding
	 *  @{
	 */
	window_size int      /**<  @brief  window size for local folding sliding window approach */
	ptype_local []string /**<  @brief  Pair type array (for local folding) */

	/**
	 *  @}
	 */
}

func default_vrna_fold_compound_t() *vrna_fold_compound_t {
	return &vrna_fold_compound_t{
		length:        0,
		strands:       0,
		cutpoint:      -1,
		strand_number: nil,
		strand_order:  nil,
		strand_start:  nil,
		strand_end:    nil,
		nucleotides:   nil,
		alignment:     nil,

		// hc:       nil,
		// matrices: nil,
		// exp_matrices: nil,
		params: nil,
		// exp_params:   nil,
		iindx: nil,
		jindx: nil,

		// stat_cb:      nil,
		// auxdata:      nil,
		// free_auxdata: nil,

		// domains_struc: nil,
		// domains_up:    nil,
		// aux_grammar:   nil,

		sequence:           "",
		sequence_encoding:  make([]int, 0),
		sequence_encoding2: nil,
		ptype:              "",
		sc:                 nil,

		// axD1:          0,
		// axD2:          0,
		reference_pt1: nil,
		reference_pt2: nil,
		referenceBPs1: nil,
		referenceBPs2: nil,
		bpdist:        nil,
		mm1:           nil,
		mm2:           nil,

		window_size: -1,
		ptype_local: nil,
	}
}

/**
 *  @brief  Data structure representing a nucleotide sequence
 */
type vrna_seq_t struct {
	sequence string /**< @brief The string representation of the sequence */
	encoding []int  /**< @brief The integer representation of the sequence */
	// encoding5 []int  /* vivek: no documentation of this variable given in the orig repo, and not used for mfe calculation so commented out */
	// encoding3 []int  /* vivek: no documentation of this variable given in the orig repo, and not used for mfe calculation so commented out */
	length uint32 /**< @brief The length of the sequence */
}

// type vrna_msa_t struct {
// 	n_seq        uint32
// 	sequences    []vrna_seq_t
// 	gapfree_seq  []string
// 	gapfree_size []int    /* for MAF alignment coordinates */
// 	genome_size  []uint32 /* for MAF alignment coordinates */
// 	start        []uint32 /* for MAF alignment coordinates */
// 	orientation  string   /* for MAF alignment coordinates */
// 	a2s          [][]uint32
// }

// type vrna_hc_t struct {
// 	n uint32

// 	state rune

// 	mx []rune

// 	matrix_local [][]rune

// 	up_ext []int /**<  @brief  A linear array that holds the number of allowed
// 	 *            unpaired nucleotides in an exterior loop
// 	 */
// 	up_hp []int /**<  @brief  A linear array that holds the number of allowed
// 	 *            unpaired nucleotides in a hairpin loop
// 	 */
// 	up_int []int /**<  @brief  A linear array that holds the number of allowed
// 	 *            unpaired nucleotides in an interior loop
// 	 */
// 	up_ml []int /**<  @brief  A linear array that holds the number of allowed
// 	 *            unpaired nucleotides in a multi branched loop
// 	 */
// }

/**
 *  @brief  Minimum Free Energy (MFE) Dynamic Programming (DP) matrices data structure required within the #vrna_fold_compound_t
 */
// type vrna_mx_mfe_t struct {
// 	/** @name Common fields for MFE matrices
// 	 *  @{
// 	 */

// 	length uint32 /**<  @brief  Length of the sequence, therefore an indicator of the size of the DP matrices */
// 	/**
// 	 *  @}
// 	 */

// 	/** @name Default DP matrices
// 	 *  @note These data fields are available if
// 	 *        @code vrna_mx_mfe_t.type == VRNA_MX_DEFAULT @endcode
// 	 * @{
// 	 */
// 	c   []int /**<  @brief  Energy array, given that i-j pair */
// 	f5  []int /**<  @brief  Energy of 5' end */
// 	f3  []int /**<  @brief  Energy of 3' end */
// 	fc  []int /**<  @brief  Energy from i to cutpoint (and vice versa if i>cut) */
// 	fML []int /**<  @brief  Multi-loop auxiliary energy array */
// 	fM1 []int /**<  @brief  Second ML array, only for unique multibrnach loop decomposition */
// 	fM2 []int /**<  @brief  Energy for a multibranch loop region with exactly two stems, extending to 3' end */
// 	ggg []int /**<  @brief  Energies of g-quadruplexes */
// 	Fc  int   /**<  @brief  Minimum Free Energy of entire circular RNA */
// 	FcH int
// 	FcI int
// 	FcM int
// 	/**
// 	 * @}
// 	 */

// 	/** @name Local Folding DP matrices using window approach
// 	 *  @note These data fields are available if
// 	 *        @code vrna_mx_mfe_t.type == VRNA_MX_WINDOW @endcode
// 	 * @{
// 	 */
// 	c_local   [][]int /**<  @brief  Energy array, given that i-j pair */
// 	f3_local  []int   /**<  @brief  Energy of 5' end */
// 	fML_local [][]int /**<  @brief  Multi-loop auxiliary energy array */
// 	ggg_local [][]int /**<  @brief  Energies of g-quadruplexes */
// 	/**
// 	 * @}
// 	 */

// 	/** @name Distance Class DP matrices
// 	 *  @note These data fields are available if
// 	 *        @code vrna_mx_mfe_t.type == VRNA_MX_2DFOLD @endcode
// 	 * @{
// 	 */
// 	E_F5     [][][]int
// 	l_min_F5 [][]int
// 	l_max_F5 [][]int
// 	k_min_F5 []int
// 	k_max_F5 []int

// 	E_F3     [][][]int
// 	l_min_F3 [][]int
// 	l_max_F3 [][]int
// 	k_min_F3 []int
// 	k_max_F3 []int

// 	E_C     [][][]int
// 	l_min_C [][]int
// 	l_max_C [][]int
// 	k_min_C []int
// 	k_max_C []int

// 	E_M     [][][]int
// 	l_min_M [][]int
// 	l_max_M [][]int
// 	k_min_M []int
// 	k_max_M []int

// 	E_M1     [][][]int
// 	l_min_M1 [][]int
// 	l_max_M1 [][]int
// 	k_min_M1 []int
// 	k_max_M1 []int

// 	E_M2     [][][]int
// 	l_min_M2 [][]int
// 	l_max_M2 [][]int
// 	k_min_M2 []int
// 	k_max_M2 []int

// 	E_Fc     [][]int
// 	l_min_Fc []int
// 	l_max_Fc []int
// 	k_min_Fc int
// 	k_max_Fc int

// 	E_FcH     [][]int
// 	l_min_FcH []int
// 	l_max_FcH []int
// 	k_min_FcH int
// 	k_max_FcH int

// 	E_FcI     [][]int
// 	l_min_FcI []int
// 	l_max_FcI []int
// 	k_min_FcI int
// 	k_max_FcI int

// 	E_FcM     [][]int
// 	l_min_FcM []int
// 	l_max_FcM []int
// 	k_min_FcM int
// 	k_max_FcM int

// 	/* auxilary arrays for remaining set of coarse graining (k,l) > (k_max, l_max) */
// 	E_F5_rem []int
// 	E_F3_rem []int
// 	E_C_rem  []int
// 	E_M_rem  []int
// 	E_M1_rem []int
// 	E_M2_rem []int

// 	E_Fc_rem  int
// 	E_FcH_rem int
// 	E_FcI_rem int
// 	E_FcM_rem int

// 	// #ifdef COUNT_STATES
// 	//   unsigned long ***N_F5;
// 	//   unsigned long ***N_C;
// 	//   unsigned long ***N_M;
// 	//   unsigned long ***N_M1;
// 	// #endif

// 	/**
// 	 * @}
// 	 */
// }

/**
 *  @brief  A base pair constraint
 */
type vrna_sc_bp_storage_t struct {
	interval_start, interval_end uint32
	e                            int
}

/**
 *  @brief  The soft constraints data structure
 *
 *  @ingroup soft_constraints
 */
type vrna_sc_t struct {
	n uint32

	state rune

	energy_up     [][]int     /**<  @brief Energy contribution for stretches of unpaired nucleotides */
	exp_energy_up [][]float64 /**<  @brief Boltzmann Factors of the energy contributions for unpaired sequence stretches */

	up_storage []int /**<  @brief  Storage container for energy contributions per unpaired nucleotide */
	// bp_storage [][]vrna_sc_bp_storage_t /**<  @brief  Storage container for energy contributions per base pair */

	energy_bp     []int     /**<  @brief Energy contribution for base pairs */
	exp_energy_bp []float64 /**<  @brief Boltzmann Factors of the energy contribution for base pairs */

	energy_bp_local     [][]int     /**<  @brief Energy contribution for base pairs (sliding window approach) */
	exp_energy_bp_local [][]float64 /**<  @brief Boltzmann Factors of the energy contribution for base pairs (sliding window approach) */

	energy_stack     []int     /**<  @brief Pseudo Energy contribution per base pair involved in a stack */
	exp_energy_stack []float64 /**<  @brief Boltzmann weighted pseudo energy contribution per nucleotide involved in a stack */

	// /* generic soft contraints below */
	// vrna_callback_sc_energy     *f;     /**<  @brief  A function pointer used for pseudo
	//                                      *            energy contribution in MFE calculations
	//                                      *    @see    vrna_sc_add_f()
	//                                      */

	// vrna_callback_sc_backtrack  *bt;    /**<  @brief  A function pointer used to obtain backtraced
	//                                      *            base pairs in loop regions that were altered
	//                                      *            by soft constrained pseudo energy contributions
	//                                      *    @see    vrna_sc_add_bt()
	//                                      */

	// vrna_callback_sc_exp_energy *exp_f; /**<  @brief  A function pointer used for pseudo energy
	//                                      *            contribution boltzmann factors in PF
	//                                      *            calculations
	//                                      *    @see    vrna_sc_add_exp_f()
	//                                      */

	// void                        *data;  /**<  @brief  A pointer to the data object provided for
	//                                      *            for pseudo energy contribution functions of the
	//                                      *            generic soft constraints feature
	//                                      */
	// vrna_callback_free_auxdata  *free_data;
}

/**
 *  @brief  A base pair constraint
 */

// /**
//  *  @brief  Partition function (PF) Dynamic Programming (DP) matrices data structure required within the #vrna_fold_compound_t
//  */
// type vrna_mx_pf_t struct {
//   /** @name Common fields for DP matrices
//    *  @{
//    */

//   length uint32
//   scale []float64
//   expMLbase []float64

//   /**
//    *  @}
//    */

//   /** @name Default PF matrices
//    *  @note These data fields are available if
//    *        @code vrna_mx_pf_t.type == VRNA_MX_DEFAULT @endcode
//    *  @{
//    */
//   q []float64
//   qb []float64
//   qm []float64
//   qm1 []float64
//   probs []float64
//   q1k []float64
//   qln []float64
//   G []float64

//   qo float64
//   qm2 []float64
//   qho float64
//   qio float64
//   qmo float64

//   /**
//    *  @}
//    */

// /** @name Local Folding DP matrices using window approach
//    *  @note These data fields are available if
//    *        @code vrna_mx_mfe_t.type == VRNA_MX_WINDOW @endcode
//    * @{
//    */
// 	 q_local [][]float64
// 	 qb_local [][]float64
// 	 qm_local [][]float64
// 	 pR [][]float64
// 	 qm2_local [][]float64
// 	 QI5 [][]float64
// 	 q2l [][]float64
// 	 qmb [][]float64
// 	 G_local [][]float64
// 	 /**
// 		*  @}
// 		*/

// 	 /** @name Distance Class DP matrices
// 		*  @note These data fields are available if
// 		*        @code vrna_mx_pf_t.type == VRNA_MX_2DFOLD @endcode
// 		*  @{
// 		*/
// 	 Q [][][]float64
// 	 l_min_Q [][]int
// 	 l_max_Q [][]int
// 	 int *k_min_Q;
// 	 int *k_max_Q;

// 	 Q_B [][][]float64
// 	 l_min_Q_B [][]int
// 	 l_max_Q_B [][]int
// 	 int *k_min_Q_B;
// 	 int *k_max_Q_B;

// 	 Q_M [][][]float64
// 	 l_min_Q_M [][]int
// 	 l_max_Q_M [][]int
// 	 int *k_min_Q_M;
// 	 int *k_max_Q_M;

// 	 Q_M1 [][][]float64
// 	 l_min_Q_M1 [][]int
// 	 l_max_Q_M1 [][]int
// 	 int *k_min_Q_M1;
// 	 int *k_max_Q_M1;

// 	 Q_M2 [][][]float64
// 	 l_min_Q_M2 [][]int
// 	 l_max_Q_M2 [][]int
// 	 int *k_min_Q_M2;
// 	 int *k_max_Q_M2;

// 	 Q_c [][]float64
// 	 int *l_min_Q_c;
// 	 int *l_max_Q_c;
// 	 int k_min_Q_c;
// 	 int k_max_Q_c;

// 	 Q_cH [][]float64
// 	 int *l_min_Q_cH;
// 	 int *l_max_Q_cH;
// 	 int k_min_Q_cH;
// 	 int k_max_Q_cH;

// 	 Q_cI [][]float64
// 	 int *l_min_Q_cI;
// 	 int *l_max_Q_cI;
// 	 int k_min_Q_cI;
// 	 int k_max_Q_cI;

// 	 Q_cM [][]float64
// 	 int *l_min_Q_cM;
// 	 int *l_max_Q_cM;
// 	 int k_min_Q_cM;
// 	 int k_max_Q_cM;

// 	 /* auxilary arrays for remaining set of coarse graining (k,l) > (k_max, l_max) */
// 	 FLT_OR_DBL *Q_rem;
// 	 FLT_OR_DBL *Q_B_rem;
// 	 FLT_OR_DBL *Q_M_rem;
// 	 FLT_OR_DBL *Q_M1_rem;
// 	 FLT_OR_DBL *Q_M2_rem;

// 	 FLT_OR_DBL Q_c_rem;
// 	 FLT_OR_DBL Q_cH_rem;
// 	 FLT_OR_DBL Q_cI_rem;
// 	 FLT_OR_DBL Q_cM_rem;
// 	 /**
// 		*  @}
// 		*/

//  #ifndef VRNA_DISABLE_C11_FEATURES
// 	 /* C11 support for unnamed unions/structs */
//  };
//  };
//  #endif
//  };

func CalculateMfe(seq, structure string) (float64, error) {
	if len(seq) != len(structure) {
		return 0, fmt.Errorf("length of sequence (%v) != Length of structure (%v)", len(seq), len(structure))
	} else if len(seq) == 0 {
		return 0, errors.New("lengths of sequence and structure cannot be 0")
	}

	// fc is everything with NIL
	// md is the default model details
	// add_params(fc, &md, options);
	/*
	 * ALWAYS provide regular energy parameters
	 * remove previous parameters if present and they differ from current model
	 */

	fc := default_vrna_fold_compound_t()
	fc.params = vrna_params()
	fc.sequence = seq
	fc.length = uint32(len(seq))
	// fc := &vrna_fold_compound_t{params: vrna_params(), sequence: seq, structure: structure}

	// vrna_params_prepare(fc, options);
	sanitize_bp_span(fc)
	set_fold_compound(fc)

	// tmp is structure

	energy, err := vrna_eval_structure_cstr(fc, structure)
	if err != nil {
		return 0, err
	}
	return energy, nil
}

func vrna_params() *vrna_param_t {
	return get_scaled_params()
}

func RESCALE_dG_int(dG, dH, dT int) int {
	return (dH - (dH-dG)*dT)
}

func RESCALE_dG_float64(dG, dH, dT float64) float64 {
	return (dH - (dH-dG)*dT)
}

var (
	DEFAULT_TEMP float64 = 37.0
	NB_PAIRS     int     = 7 /** The number of distinguishable base pairs */
)

const (
	MAXALPHA int = 7  /** @brief Maximal length of alphabet */
	MAXLOOP  int = 30 /** The maximum loop length */

)

type vrna_md_t struct {
	temperature float64 /**<  @brief  The temperature used to scale the thermodynamic parameters */
	betaScale   float64 /**<  @brief  A scaling factor for the thermodynamic temperature of the Boltzmann factors */
	pf_smooth   int     /**<  @brief  A flat specifying whether energies in Boltzmann factors need to be smoothed */
	dangles     int     /**<  @brief  Specifies the dangle model used in any energy evaluation (0,1,2 or 3)
	 *
	 *    If set to 0 no stabilizing energies are assigned to bases adjacent to
	 *    helices in free ends and multiloops (so called dangling ends). Normally
	 *    (dangles = 1) dangling end energies are assigned only to unpaired
	 *    bases and a base cannot participate simultaneously in two dangling ends. In
	 *    the partition function algorithm vrna_pf() these checks are neglected.
	 *    To provide comparability between free energy minimization and partition function
	 *    algorithms, the default setting is 2.
	 *    This treatment of dangling ends gives more favorable energies to helices
	 *    directly adjacent to one another, which can be beneficial since such
	 *    helices often do engage in stabilizing interactions through co-axial
	 *    stacking.\n
	 *    If set to 3 co-axial stacking is explicitly included for
	 *    adjacent helices in multiloops. The option affects only mfe folding
	 *    and energy evaluation (vrna_mfe() and vrna_eval_structure()), as
	 *    well as suboptimal folding (vrna_subopt()) via re-evaluation of energies.
	 *    Co-axial stacking with one intervening mismatch is not considered so far.
	 *    @note   Some function do not implement all dangle model but only a subset of
	 *            (0,1,2,3). In particular, partition function algorithms can only handle
	 *            0 and 2. Read the documentation of the particular recurrences or
	 *            energy evaluation function for information about the provided dangle
	 *            model.
	 */
	special_hp     int  /**<  @brief  Include special hairpin contributions for tri, tetra and hexaloops */
	noLP           int  /**<  @brief  Only consider canonical structures, i.e. no 'lonely' base pairs */
	noGU           int  /**<  @brief  Do not allow GU pairs */
	noGUclosure    int  /**<  @brief  Do not allow loops to be closed by GU pair */
	logML          int  /**<  @brief  Use logarithmic scaling for multiloops */
	circ           int  /**<  @brief  Assume RNA to be circular instead of linear */
	gquad          int  /**<  @brief  Include G-quadruplexes in structure prediction */
	uniq_ML        int  /**<  @brief  Flag to ensure unique multi-branch loop decomposition during folding */
	energy_set     int  /**<  @brief  Specifies the energy set that defines set of compatible base pairs */
	backtrack      int  /**<  @brief  Specifies whether or not secondary structures should be backtraced */
	backtrack_type rune /**<  @brief  Specifies in which matrix to backtrack */
	compute_bpp    int  /**<  @brief  Specifies whether or not backward recursions for base pair probability (bpp) computation will be performed */
	// nonstandards   [64]rune /**<  @brief  contains allowed non standard bases */
	max_bp_span int /**<  @brief  maximum allowed base pair span */

	min_loop_size int /**<  @brief  Minimum size of hairpin loops
	 *    @note The default value for this field is #TURN, however, it may
	 *    be 0 in cofolding context.
	 */
	window_size int     /**<  @brief  Size of the sliding window for locally optimal structure prediction */
	oldAliEn    int     /**<  @brief  Use old alifold energy model */
	ribo        int     /**<  @brief  Use ribosum scoring table in alifold energy model */
	cv_fact     float64 /**<  @brief  Co-variance scaling factor for consensus structure prediction */
	nc_fact     float64 /**<  @brief  Scaling factor to weight co-variance contributions of non-canonical pairs */
	sfact       float64 /**<  @brief  Scaling factor for partition function scaling */
	rtype       [8]int  /**<  @brief  Reverse base pair type array */
	// vivek: values below are MAXALPHA + 1. Just evaluating to int because it gives errors
	alias [MAXALPHA + 1]int               /**<  @brief  alias of an integer nucleotide representation */
	pair  [MAXALPHA + 1][MAXALPHA + 1]int /**<  @brief  Integer representation of a base pair */
}

var (
	// Default temperature for structure prediction and free energy evaluation in $^\circ C$
	VRNA_MODEL_DEFAULT_TEMPERATURE float64 = 37.0
	VRNA_MODEL_DEFAULT_PF_SMOOTH   int     = 1
	// Default dangling end model
	VRNA_MODEL_DEFAULT_DANGLES int = 2
	// Default model behavior for lookup of special tri-, tetra-, and hexa-loops
	VRNA_MODEL_DEFAULT_SPECIAL_HP int = 1
	// Default model behavior for so-called 'lonely pairs'
	VRNA_MODEL_DEFAULT_NO_LP int = 1
	// Default model behavior for G-U base pairs
	VRNA_MODEL_DEFAULT_NO_GU int = 0
	// Default model behavior for G-U base pairs closing a loop
	VRNA_MODEL_DEFAULT_NO_GU_CLOSURE int = 0
	// Default model behavior on how to evaluate the energy contribution of multi-branch loops
	VRNA_MODEL_DEFAULT_LOG_ML int = 0
	// Default model behavior to treat a molecule as a circular RNA (DNA)
	VRNA_MODEL_DEFAULT_CIRC int = 0
	// Default model behavior regarding the treatment of G-Quadruplexes
	VRNA_MODEL_DEFAULT_GQUAD int = 0
	// Default behavior of the model regarding unique multi-branch loop decomposition
	VRNA_MODEL_DEFAULT_UNIQ_ML int = 0
	// Default model behavior on which energy set to use
	VRNA_MODEL_DEFAULT_ENERGY_SET int = 0
	// Default model behavior with regards to backtracking of structures
	VRNA_MODEL_DEFAULT_BACKTRACK int = 1
	// Default model behavior on what type of backtracking to perform
	VRNA_MODEL_DEFAULT_BACKTRACK_TYPE rune = 'F'
	// Default model behavior with regards to computing base pair probabilities
	VRNA_MODEL_DEFAULT_COMPUTE_BPP int = 1
	// Default model behavior for the allowed maximum base pair span
	VRNA_MODEL_DEFAULT_MAX_BP_SPAN int = -1
	// The minimum loop length
	TURN int = 3
	// Default model behavior for the sliding window approach
	VRNA_MODEL_DEFAULT_WINDOW_SIZE int = -1
	// Default model behavior for consensus structure energy evaluation
	VRNA_MODEL_DEFAULT_ALI_OLD_EN int = 0
	// Default model behavior for consensus structure co-variance contribution assessment
	VRNA_MODEL_DEFAULT_ALI_RIBO int = 0
	// Default model behavior for weighting the co-variance score in consensus structure prediction
	VRNA_MODEL_DEFAULT_ALI_CV_FACT float64 = 1.0
	// Default model behavior for weighting the nucleotide conservation? in consensus structure prediction
	VRNA_MODEL_DEFAULT_ALI_NC_FACT float64 = 1.0
)

const NBPAIRS = 7 /** The number of distinguishable base pairs */

func DefaultModelDetails() *vrna_md_t {
	return &vrna_md_t{
		temperature:    VRNA_MODEL_DEFAULT_TEMPERATURE,
		betaScale:      1.,
		pf_smooth:      VRNA_MODEL_DEFAULT_PF_SMOOTH,
		dangles:        VRNA_MODEL_DEFAULT_DANGLES,
		special_hp:     VRNA_MODEL_DEFAULT_SPECIAL_HP,
		noLP:           VRNA_MODEL_DEFAULT_NO_LP,
		noGU:           VRNA_MODEL_DEFAULT_NO_GU,
		noGUclosure:    VRNA_MODEL_DEFAULT_NO_GU_CLOSURE,
		logML:          VRNA_MODEL_DEFAULT_LOG_ML,
		circ:           VRNA_MODEL_DEFAULT_CIRC,
		gquad:          VRNA_MODEL_DEFAULT_GQUAD,
		uniq_ML:        VRNA_MODEL_DEFAULT_UNIQ_ML,
		energy_set:     VRNA_MODEL_DEFAULT_ENERGY_SET,
		backtrack:      VRNA_MODEL_DEFAULT_BACKTRACK,
		backtrack_type: VRNA_MODEL_DEFAULT_BACKTRACK_TYPE,
		compute_bpp:    VRNA_MODEL_DEFAULT_COMPUTE_BPP,
		// nonstandards:   []rune{0},
		max_bp_span:   VRNA_MODEL_DEFAULT_MAX_BP_SPAN,
		min_loop_size: TURN,
		window_size:   VRNA_MODEL_DEFAULT_WINDOW_SIZE,
		oldAliEn:      VRNA_MODEL_DEFAULT_ALI_OLD_EN,
		ribo:          VRNA_MODEL_DEFAULT_ALI_RIBO,
		cv_fact:       VRNA_MODEL_DEFAULT_ALI_CV_FACT,
		nc_fact:       VRNA_MODEL_DEFAULT_ALI_NC_FACT,
		sfact:         1.07,
		rtype:         [8]int{0, 2, 1, 4, 3, 6, 5, 7},
		alias:         [MAXALPHA + 1]int{0, 1, 2, 3, 4, 3, 2, 0},
		pair: [MAXALPHA + 1][MAXALPHA + 1]int{
			{0, 0, 0, 0, 0, 0, 0, 0},
			{0, 0, 0, 0, 5, 0, 0, 5},
			{0, 0, 0, 1, 0, 0, 0, 0},
			{0, 0, 2, 0, 3, 0, 0, 0},
			{0, 6, 0, 4, 0, 0, 0, 6},
			{0, 0, 0, 0, 0, 0, 2, 0},
			{0, 0, 0, 0, 0, 1, 0, 0},
			{0, 6, 0, 0, 5, 0, 0, 0},
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
	// Tetraloops                                [1401]rune
	Tetraloops string
	Triloop_E  [40]int
	// Triloops                                  [241]rune
	Triloops   string
	Hexaloop_E [40]int
	Hexaloops  string
	// Hexaloops                                 [1801]rune
	TripleC, MultipleCA, MultipleCB           int
	gquad                                     [VRNA_GQUAD_MAX_STACK_SIZE + 1][3*VRNA_GQUAD_MAX_LINKER_LENGTH + 1]int
	gquadLayerMismatch, gquadLayerMismatchMax int
	temperature                               float64    /**<  @brief  Temperature used for loop contribution scaling */
	model_details                             *vrna_md_t /**<  @brief  Model details to be used in the recursions */
}

func get_scaled_params() *vrna_param_t {
	var model_details *vrna_md_t = DefaultModelDetails()
	// float64
	var i, j, k, l int
	var tempf float64 = (model_details.temperature + K0) / Tmeasure

	var params *vrna_param_t = &vrna_param_t{
		model_details:         model_details,
		temperature:           model_details.temperature,
		lxc:                   lxc37 * tempf,
		TripleC:               RESCALE_dG_int(TripleC37, TripleCdH, int(tempf)),
		MultipleCA:            RESCALE_dG_int(MultipleCA37, MultipleCAdH, int(tempf)),
		MultipleCB:            RESCALE_dG_int(TerminalAU37, TerminalAUdH, int(tempf)),
		TerminalAU:            RESCALE_dG_int(TerminalAU37, TerminalAUdH, int(tempf)),
		DuplexInit:            RESCALE_dG_int(DuplexInit37, DuplexInitdH, int(tempf)),
		MLbase:                RESCALE_dG_int(ML_BASE37, ML_BASEdH, int(tempf)),
		MLclosing:             RESCALE_dG_int(ML_closing37, ML_closingdH, int(tempf)),
		gquadLayerMismatch:    RESCALE_dG_int(GQuadLayerMismatch37, GQuadLayerMismatchH, int(tempf)),
		gquadLayerMismatchMax: GQuadLayerMismatchMax,
	}

	params.ninio[2] = RESCALE_dG_int(ninio37, niniodH, int(tempf))

	for i = VRNA_GQUAD_MIN_STACK_SIZE; i <= VRNA_GQUAD_MAX_STACK_SIZE; i++ {
		for j = 3 * VRNA_GQUAD_MIN_LINKER_LENGTH; j <= 3*VRNA_GQUAD_MAX_LINKER_LENGTH; j++ {
			var GQuadAlpha_T float64 = float64(RESCALE_dG_int(GQuadAlpha37, GQuadAlphadH, int(tempf)))
			var GQuadBeta_T float64 = float64(RESCALE_dG_int(GQuadBeta37, GQuadBetadH, int(tempf)))
			params.gquad[i][j] = int(GQuadAlpha_T)*(i-1) + int(float64(GQuadBeta_T)*math.Log(float64(j)-2.0))
		}
	}

	for i = 0; i < 31; i++ {
		params.hairpin[i] = RESCALE_dG_int(hairpin37[i], hairpindH[i], int(tempf))
	}

	for i = 0; i <= Min(30, MAXLOOP); i++ {
		params.bulge[i] = RESCALE_dG_int(bulge37[i], bulgedH[i], int(tempf))
		params.internal_loop[i] = RESCALE_dG_int(internal_loop37[i], internal_loopdH[i], int(tempf))
	}

	for ; i <= MAXLOOP; i++ {
		params.bulge[i] = params.bulge[30] +
			int(params.lxc*math.Log(float64(i)/30.0))
		params.internal_loop[i] = params.internal_loop[30] +
			int(params.lxc*math.Log(float64(i)/30.0))
	}

	// This could be 281 or only the amount of chars
	for i = 0; (i * 7) < len(Tetraloops); i++ {
		params.Tetraloop_E[i] = RESCALE_dG_int(Tetraloop37[i], TetraloopdH[i], int(tempf))
	}

	for i = 0; (i * 5) < len(Triloops); i++ {
		params.Triloop_E[i] = RESCALE_dG_int(Triloop37[i], TriloopdH[i], int(tempf))
	}

	for i = 0; (i * 9) < len(Hexaloops); i++ {
		params.Hexaloop_E[i] = RESCALE_dG_int(Hexaloop37[i], HexaloopdH[i], int(tempf))
	}

	for i = 0; i <= NBPAIRS; i++ {
		params.MLintern[i] = RESCALE_dG_int(ML_intern37, ML_interndH, int(tempf))
	}

	/* stacks    G(T) = H - [H - G(T0)]*T/T0 */
	for i = 0; i <= NBPAIRS; i++ {
		for j = 0; j <= NBPAIRS; j++ {
			params.stack[i][j] = RESCALE_dG_int(stack37[i][j],
				stackdH[i][j],
				int(tempf))
		}
	}

	/* mismatches */
	for i = 0; i <= NBPAIRS; i++ {
		for j = 0; j < 5; j++ {
			for k = 0; k < 5; k++ {
				var mm int
				params.mismatchI[i][j][k] = RESCALE_dG_int(mismatchI37[i][j][k],
					mismatchIdH[i][j][k],
					int(tempf))
				params.mismatchH[i][j][k] = RESCALE_dG_int(mismatchH37[i][j][k],
					mismatchHdH[i][j][k],
					int(tempf))
				params.mismatch1nI[i][j][k] = RESCALE_dG_int(mismatch1nI37[i][j][k],
					mismatch1nIdH[i][j][k],
					int(tempf))
				params.mismatch23I[i][j][k] = RESCALE_dG_int(mismatch23I37[i][j][k],
					mismatch23IdH[i][j][k],
					int(tempf))
				if model_details.dangles > 0 {
					mm = RESCALE_dG_int(mismatchM37[i][j][k],
						mismatchMdH[i][j][k],
						int(tempf))
					if mm > 0 {
						params.mismatchM[i][j][k] = 0
					} else {
						params.mismatchM[i][j][k] = mm
					}

					mm = RESCALE_dG_int(mismatchExt37[i][j][k],
						mismatchExtdH[i][j][k],
						int(tempf))
					if mm > 0 {
						params.mismatchExt[i][j][k] = 0
					} else {
						params.mismatchExt[i][j][k] = mm
					}
				} else {
					params.mismatchExt[i][j][k] = 0
					params.mismatchM[i][j][k] = 0
				}
			}
		}
	}

	/* dangles */
	for i = 0; i <= NBPAIRS; i++ {
		for j = 0; j < 5; j++ {
			var dd int
			dd = RESCALE_dG_int(dangle5_37[i][j],
				dangle5_dH[i][j],
				int(tempf))
			if dd > 0 {
				params.dangle5[i][j] = 0
			} else {
				params.dangle5[i][j] = dd
			}

			dd = RESCALE_dG_int(dangle3_37[i][j],
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
					params.int11[i][j][k][l] = RESCALE_dG_int(int11_37[i][j][k][l],
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
						params.int21[i][j][k][l][m] = RESCALE_dG_int(int21_37[i][j][k][l][m],
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
							params.int22[i][j][k][l][m][n] = RESCALE_dG_int(int22_37[i][j][k][l][m][n],
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
	pt, err := vrna_ptable_from_string(structure)
	if err != nil {
		return 0, err
	}

	// var en float64
	en, err := wrap_eval_structure(fc, structure, pt)
	if err != nil {
		return 0, err
	}

	return en, nil
}

var SHRT_MAX uint32 = 32767

func vrna_ptable_from_string(structure string) ([]int, error) {
	var pt []int
	var n uint32

	n = uint32(len(structure))

	if n > SHRT_MAX {
		err := fmt.Sprintf("vrna_ptable_from_string: Structure too long to be converted to pair table (n=%v, max=%v)", n, SHRT_MAX)
		return nil, errors.New(err)
	}

	pt = make([]int, n+2)
	pt[0] = int(n)

	pt, err := extract_pairs(pt, structure, []rune("()"))
	if err != nil {
		return nil, err
	}

	return pt, nil
}

/* requires that pt[0] already contains the length of the string! */
func extract_pairs(pt []int, structure string, pair []rune) ([]int, error) {
	var ptr []rune = []rune(structure)
	var open, close rune
	var stack []int
	var i, j, n uint32
	var hx, ptr_idx int

	n = uint32(pt[0])
	stack = make([]int, n+1)

	open = pair[0]
	close = pair[1]

	for hx, i, ptr_idx = 0, 1, 0; i <= n && ptr_idx < len(ptr); ptr_idx++ {
		if ptr[ptr_idx] == open {
			stack[hx] = int(i)
			hx++
		} else if ptr[ptr_idx] == close {
			hx--
			j = uint32(stack[hx])

			if hx < 0 {
				// vrna_message_warning("%s\nunbalanced brackets '%2s' found while extracting base pairs",
				//                      );
				return nil, errors.New(fmt.Sprintf("%v\nunbalanced brackets '%v' found while extracting base pairs", structure,
					pair))
				// free(stack);
				// return 0;
			}

			pt[i] = int(j)
			pt[j] = int(i)
		}
		i++
	}

	// free(stack);

	if hx != 0 {
		return nil, errors.New(fmt.Sprintf("%v\nunbalanced brackets '%v' found while extracting base pairs", structure,
			pair))
		// return 0;
	}

	return pt, nil
	// return 1; /* success */
}

func wrap_eval_structure(fc *vrna_fold_compound_t, structure string, pt []int) (float64, error) {
	var res, gq int
	// var l [3]int
	var energy float64

	energy = float64(INF) / 100.0

	gq = fc.params.model_details.gquad
	fc.params.model_details.gquad = 0

	// can add support for circular strands with this
	// if (vc.params.model_details.circ)
	//   res = eval_circ_pt(vc, pt, output_stream, verbosity);
	// else
	res = eval_pt(fc, pt)
	fc.params.model_details.gquad = gq

	// if gq == 1 && (parse_gquad(structure, &L, l) > 0) {
	// 	if (verbosity > 0)
	// 		vrna_cstr_print_eval_sd_corr(output_stream);

	// 	res += en_corr_of_loop_gquad(vc, 1, vc.length, structure, pt, output_stream, verbosity);
	// }
	energy = float64(res) / 100.0
	return energy, nil
}

func eval_pt(fc *vrna_fold_compound_t, pt []int) int {
	// var sn []uint32
	var i, length, energy int

	length = int(fc.length)

	// vrna_sc_prepare(fc)
	// if fc.params.model_details.backtrack_type == 'M' {
	// 	energy = energy_of_ml_pt(vc, 0, pt)
	// } else {

	// }
	// panic(fmt.Sprintf("%v %v", fc.length, len(pt)))
	energy = energy_of_extLoop_pt(fc, 0, pt)
	// panic(fmt.Sprintf("fc: %v \n\n\n pt: %v", fc, pt))
	// energy = 0

	for i = 1; i <= length; i++ {
		if pt[i] == 0 {
			continue
		}

		energy += stack_energy(fc, i, pt)
		i = pt[i]
	}
	for i = 1; fc.strand_number[i] != fc.strand_number[length]; i++ {
		if fc.strand_number[i] != fc.strand_number[pt[i]] {
			energy += fc.params.DuplexInit
			break
		}
	}

	return energy
}

// func vrna_params_prepare(fc *vrna_fold_compound_t) {
// // VRNA_OPTION_MFE | VRNA_OPTION_EVAL_ONLY

//     var md_p *vrna_md_t

//     /*
//      *  every vrna_fold_compound_t must have a vrna_paramt_t structure attached
//      *  to it that holds the current model details. So we just use this here as
//      *  the reference model
//      */
//     md_p = &(fc.params.model_details)

//     if (options & VRNA_OPTION_PF)
//     {
//       /* remove previous parameters if present and they differ from reference model */
//       if (fc.exp_params)
//       {
//         if (memcmp(md_p, &(fc.exp_params.model_details), sizeof(vrna_md_t)) != 0)
//         {
//           free(fc.exp_params);
//           fc.exp_params = NULL;
//         }
//       }

//       if (!fc.exp_params)
//         fc.exp_params = (fc.type == VRNA_FC_TYPE_SINGLE) ? vrna_exp_params(md_p) : vrna_exp_params_comparative(fc.n_seq, md_p);
//     }

// }

func sanitize_bp_span(fc *vrna_fold_compound_t) {
	var md *vrna_md_t
	md = fc.params.model_details

	/* non-local fold mode */
	md.window_size = int(fc.length)

	if md.max_bp_span <= 0 || md.max_bp_span > md.window_size {
		md.max_bp_span = md.window_size
	}
}

func set_fold_compound(fc *vrna_fold_compound_t) {
	// VRNA_OPTION_MFE | VRNA_OPTION_EVAL_ONLY
	var sequence string
	// var sequences []string
	// var length, s int
	var md_p *vrna_md_t

	md_p = fc.params.model_details

	sequence = fc.sequence

	fc.sequence = ""
	fc.length = 0

	/* split input sequences at default delimiter '&' */
	// sequences = vrna_strsplit(sequence, NULL);

	vrna_sequence_add(fc, sequence)

	if fc.strands > 1 {
		fc.cutpoint = int(fc.nucleotides[0].length) + 1

		if md_p.min_loop_size == TURN {
			md_p.min_loop_size = 0 /* is it safe to set this here? */
		}
	}

	// if !(options & VRNA_OPTION_EVAL_ONLY) {
	// 	fc.ptype = (aux & WITH_PTYPE) ? vrna_ptypes(fc.sequence_encoding2, md_p) : NULL;
	// 	/* backward compatibility ptypes */
	// 	fc.ptype_pf_compat =
	// 			(aux & WITH_PTYPE_COMPAT) ? get_ptypes(fc.sequence_encoding2, md_p, 1) : NULL;
	// }

	vrna_sequence_prepare(fc)

	fc.iindx = vrna_idx_row_wise(fc.length)
	fc.jindx = vrna_idx_col_wise(fc.length)

}

// vivek: big potential point of failure fiven the cryptic memcpy-ing
func vrna_sequence_add(fc *vrna_fold_compound_t, sequence string) {
	// VRNA_SEQUENCE_RNA
	var add_length uint32

	add_length = uint32(len(sequence))

	/* add the sequence to the nucleotides container */
	// fc.nucleotides = (vrna_seq_t *)vrna_realloc(vc.nucleotides,
	// 																							sizeof(vrna_seq_t) *
	// 																							(vc.strands + 1));
	fc.nucleotides = make([]vrna_seq_t, fc.strands+1)

	set_sequence(&fc.nucleotides[fc.strands],
		sequence,
		fc.params.model_details)

	/* increase strands counter */
	fc.strands++

	/* add new sequence to initial order of all strands */
	// fc.sequence = (char *)vrna_realloc(vc.sequence,
	// 																		sizeof(char) *
	// 																		(vc.length + add_length + 1));
	// vivek: should we add an ampersand between sequences?
	fc.sequence = fc.sequence + sequence

	/* add encoding for new strand */
	// fc.sequence_encoding = make([]int, fc.length+add_length+2)

	// fc.sequence_encoding = fc.sequence_encoding + '\0' + fc.nucleotides[fc.strands - 1].encoding

	fc.sequence_encoding = append(fc.sequence_encoding, fc.nucleotides[fc.strands-1].encoding[1:]...)

	/* restore circular encoding */
	fc.sequence_encoding = append(fc.sequence_encoding, fc.sequence_encoding[1])
	fmt.Println(len(fc.sequence_encoding))
	fc.sequence_encoding[0] = fc.sequence_encoding[len(fc.sequence_encoding)-1]

	/* add encoding2 (simple encoding) for new strand */
	// fc.sequence_encoding2 = make([]int, fc.length+add_length+2)

	enc := vrna_seq_encode_simple(fc.nucleotides[fc.strands-1].sequence,
		fc.params.model_details)
	fc.sequence_encoding2 = append(fc.sequence_encoding2, enc[1:]...)
	// panic(fmt.Sprintf("%v", fc.sequence_encoding2))
	// fc.sequence_encoding2[fc.length+1:] = enc[1 : 1+add_length]
	// memcpy(vc.sequence_encoding2 + vc.length + 1,
	// 				enc + 1,
	// 				add_length * sizeof(short))

	fc.sequence_encoding2[len(fc.sequence_encoding2)-1] = fc.sequence_encoding2[1]
	fc.sequence_encoding2[0] = int(fc.length + add_length)

	/* finally, increase length property of the fold compound */
	fc.length = fc.length + add_length
}

/* vivek:
* Sets the `encoding3` and `encoding5` variables of a vrna_seq_t obj. Setting of the `encoding` field is explained in the docs of `vrna_seq_encode`.
*
* The original repo has additional functionality for a circular sequence which has been omitted from this function.
 */
func set_sequence(obj *vrna_seq_t, sequence string, md *vrna_md_t) {
	obj.sequence = strings.ToUpper(sequence)
	obj.length = uint32(len(sequence))

	obj.encoding = vrna_seq_encode(obj.sequence, md)
	// vivek: code commented out below is sets the `encoding3` and `encoding5` variables
	// but since they are not used for mfe calculations, I've commented them out.
	// obj.encoding3 = make([]int, obj.length+1)
	// obj.encoding5 = make([]int, obj.length+1)

	// if md.circ == 1 {
	// 	panic("encountered a circular sequence.")
	// } else {
	// 	// vivek: what is special about encoding5[1] and encoding3[obj.length]?
	// 	// what is encoding5[0]?
	// 	// where is encoding5 used?
	// 	obj.encoding5[1] = 0
	// 	obj.encoding3[obj.length] = 0
	// }

	// var i uint32
	// assign `encoding5` based on encoding
	// if
	// iterate through `encoding` and assign
	// for i = 1; i < obj.length; i++ {
	// 	if obj.encoding[i] == 0 {
	// 		obj.encoding5[i+1] = obj.encoding5[i]
	// 	} else {
	// 		obj.encoding5[i+1] = obj.encoding[i]
	// 	}
	// }

	// for i := obj.length; i > 1; i-- {
	// 	if obj.encoding[i] == 0 {
	// 		obj.encoding3[i-1] = obj.encoding3[i]
	// 	} else {
	// 		obj.encoding3[i-1] = obj.encoding[i]
	// 	}
	// }
}

/* vivek:
* returns the encoded sequence of `sequence` (see doc of `vrna_seq_encode_simple` for more
* details of the encoding process)
 */
func vrna_seq_encode(sequence string, md *vrna_md_t) []int {
	var i, l uint32
	var S []int

	if (sequence != "") && md != nil {
		S = vrna_seq_encode_simple(sequence, md)

		l = uint32(len(sequence))

		for i = 1; i <= l; i++ {
			S[i] = md.alias[S[i]]
		}

		S[l+1] = S[1]
		S[0] = S[l]
	}

	return S
}

// vivek:
// encodes `sequence` into `S[1:len(sequence)-1]`
// `S[0]` is the length of `sequence`
// `S[len(sequence) + 1]` or `S[-1]` == `S[1]` which makes the sequence circular
// not sure why the orig repo makes the sequence cirular
func vrna_seq_encode_simple(sequence string, md *vrna_md_t) []int {
	var i, l uint32
	var S []int

	if (sequence != "") && md != nil {
		l = uint32(len(sequence))
		S = make([]int, l+2)

		for i = 1; i <= l; i++ { /* make numerical encoding of sequence */
			S[i] = vrna_nucleotide_encode(([]rune(sequence))[i-1], md)
		}

		S[l+1] = S[1]
		S[0] = int(l)
	}

	// panic(fmt.Sprintf("%v", S))
	return S
}

// vivek: I don't know if this variable name is a joke or not...
// var Law_and_Order []rune = []rune("_ACGUTXKI")
var Law_and_Order string = "_ACGUTXKI"

/* return numerical representation of nucleotide used e.g. in vrna_md_t.pair[][] */ //
func vrna_nucleotide_encode(c rune, md *vrna_md_t) int {

	var code int = -1

	c = unicode.ToUpper(c)

	if md != nil {
		// vivek: md.energy_set defaults to 0
		if md.energy_set > 0 {
			code = int(c-'A') + 1
		} else {
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
		}
	}

	return code
}

func vrna_sequence_prepare(fc *vrna_fold_compound_t) {
	var strand_number, i uint32

	fc.strand_order = nil
	fc.strand_start = nil
	fc.strand_end = nil
	fc.strand_number = make([]uint32, fc.length+2)

	/* 1. store initial strand order */
	fc.strand_order = make([]uint32, fc.strands+1)
	for strand_number = 0; strand_number < fc.strands; strand_number++ {
		fc.strand_order[strand_number] = strand_number
	}

	/* 2. mark start and end positions of sequences */
	fc.strand_start = make([]uint32, fc.strands+1)
	fc.strand_end = make([]uint32, fc.strands+1)

	fc.strand_start[0] = 1
	fc.strand_end[0] = fc.strand_start[0] + fc.nucleotides[0].length - 1

	for strand_number = 1; strand_number < fc.strands; strand_number++ {
		fc.strand_start[strand_number] = fc.strand_end[strand_number-1] + 1
		fc.strand_end[strand_number] = fc.strand_start[strand_number] + fc.nucleotides[strand_number].length - 1
		for i = fc.strand_start[strand_number]; i <= fc.strand_end[strand_number]; i++ {
			fc.strand_number[i] = strand_number
		}
	}

	/* this sets pos. n + 1 as well */
	fc.strand_number[fc.length+1] = fc.strands - 1
}

func vrna_idx_row_wise(length uint32) []int {
	var i uint32
	var idx []int = make([]int, length+1)

	for i = 1; i <= length; i++ {
		idx[i] = int((((length + 1 - i) * (length - i)) / 2) + length + 1)
	}
	return idx
}

func vrna_idx_col_wise(length uint32) []int {
	var i uint32
	var idx []int = make([]int, length+1)

	for i = 1; i <= length; i++ {
		idx[i] = int((i * (i - 1)) / 2)
	}
	return idx
}

func energy_of_extLoop_pt(fc *vrna_fold_compound_t, i int, pt []int) int {
	// var sn []uint32
	var energy, mm5, mm3, bonus, p, q, q_prev, length, start int
	// var s, s1 []int
	// var S, S5, S3 [][]int
	// var a2s [][]uint32
	var P *vrna_param_t
	var md *vrna_md_t
	var sc *vrna_sc_t
	// var scs []*vrna_sc_t

	/* initialize vars */
	length = int(fc.length)
	P = fc.params
	md = P.model_details
	// dangle_model = md.dangles
	sc = fc.sc

	energy = 0
	bonus = 0
	p = 1
	start = 1
	q_prev = -1

	/* seek to opening base of first stem */
	for p <= length && (pt[p] == 0) {
		p++
	}

	// vivek: do we ever reach here? sc is not init-ed anywhere. Does C init it automatically?
	/* add soft constraints for first unpaired nucleotides */
	if sc != nil {
		if sc.energy_up != nil {
			bonus += sc.energy_up[start][p-start]
		}
		/* how do we handle generalized soft constraints here ? */
	}

	for p < length {
		var tt uint32
		/* p must have a pairing partner */
		q = int(pt[p])

		/* get type of base pair (p,q) */
		tt = vrna_get_ptype_md(fc.sequence_encoding2[p], fc.sequence_encoding2[q], md)
		if (fc.strand_number[p-1] == fc.strand_number[p]) && (p > 1) {
			mm5 = fc.sequence_encoding[p-1]
		} else {
			mm5 = -1
		}

		if (fc.strand_number[q] == fc.strand_number[q+1]) && (q < length) {
			mm3 = fc.sequence_encoding[q+1]
		} else {
			mm3 = -1
		}
		energy += vrna_E_ext_stem(tt, mm5, mm3, P)

		/* seek to the next stem */
		p = q + 1
		q_prev = q
		for p <= length && pt[p] == 0 {
			p++
		}

		/* add soft constraints for unpaired region */
		if sc != nil && (q_prev+1 <= length) {
			if sc.energy_up != nil {
				bonus += sc.energy_up[q_prev+1][p-q_prev-1]
			}
			/* how do we handle generalized soft constraints here ? */
		}

		if p == i {
			break /* cut was in loop */
		}
	}

	return energy + bonus
}

func vrna_E_ext_stem(type_1 uint32, n5d int, n3d int, p *vrna_param_t) int {
	var energy int = 0

	if n5d >= 0 && n3d >= 0 {
		energy += p.mismatchExt[type_1][n5d][n3d]
	} else if n5d >= 0 {
		energy += p.dangle5[type_1][n5d]
	} else if n3d >= 0 {
		energy += p.dangle3[type_1][n3d]
	}

	if type_1 > 2 {
		energy += p.TerminalAU
	}

	return energy
}

func stack_energy(fc *vrna_fold_compound_t, i int, pt []int) int {
	/* recursively calculate energy of substructure enclosed by (i,j) */
	var sn []uint32
	// so, ss
	var ee, energy, j, p, q int
	var sequence string
	var s []int
	var P *vrna_param_t
	var md *vrna_md_t

	sn = fc.strand_number
	// so = fc.strand_order
	// ss = fc.strand_start
	s = fc.sequence_encoding2
	P = fc.params
	md = P.model_details
	energy = 0

	j = pt[i]

	sequence = fc.sequence
	panic(fmt.Sprintf("%v", s))
	if md.pair[s[i]][s[j]] == 0 {
		panic(fmt.Sprintf("bases %v and %v (%v%v) can't pair!",
			i, j,
			string(sequence[i-1]),
			string(sequence[j-1])))
	}

	p = i
	q = j

	for p < q {
		/* process all stacks and interior loops */
		p++
		q--

		for pt[p] == 0 {
			p++
		}

		for pt[q] == 0 {
			q--
		}

		if (pt[q] != int(p)) || (p > q) {
			break
		}

		ee = 0

		if md.pair[s[q]][s[p]] == 0 {
			// if (verbosity_level > VRNA_VERBOSITY_QUIET) {
			panic(fmt.Sprintf("bases %d and %d (%c%c) can't pair!",
				p, q,
				sequence[p-1],
				sequence[q-1]))
			// }
		}

		ee = vrna_eval_int_loop(fc, i, j, p, q)

		// if (verbosity_level > 0) {
		//   vrna_cstr_print_eval_int_loop(output_stream,
		//                                 i, j,
		//                                 string[i - 1], string[j - 1],
		//                                 p, q,
		//                                 string[p - 1], string[q - 1],
		//                                 (vc.type == VRNA_FC_TYPE_COMPARATIVE) ?
		//                                 (int)ee / (int)vc.n_seq :
		//                                 ee);
		// }

		energy += ee
		i = p
		j = q
	} /* end while */

	/* p,q don't pair must have found hairpin or multiloop */

	if p > q {
		/* hairpin */
		ee = vrna_eval_hp_loop(fc, i, j)
		energy += ee

		// if (verbosity_level > 0) {
		//   vrna_cstr_print_eval_hp_loop(output_stream,
		//                                i, j,
		//                                string[i - 1], string[j - 1],
		//                                (vc.type == VRNA_FC_TYPE_COMPARATIVE) ?
		//                                (int)ee / (int)vc.n_seq :
		//                                ee);
		// }

		return energy
	}

	/* (i,j) is exterior pair of multiloop */
	for p < j {
		/* add up the contributions of the substructures of the ML */
		energy += stack_energy(fc, p, pt)
		p = pt[p]
		/* search for next base pair in multiloop */
		p++
		for pt[p] == 0 {
			p++
		}
	}

	ee = 0
	var err error

	var ii int = cut_in_loop(i, pt, sn)
	if ii == 0 {
		ee, err = energy_of_ml_pt(fc, i, pt)
		if err != nil {
			panic(err)
		}
	} else {
		ee = energy_of_extLoop_pt(fc, ii, pt)
	}

	energy += ee
	// if (verbosity_level > 0) {
	//   vrna_cstr_print_eval_mb_loop(output_stream,
	//                                i, j,
	//                                string[i - 1], string[j - 1],
	//                                (vc.type == VRNA_FC_TYPE_COMPARATIVE) ?
	//                                (int)ee / (int)vc.n_seq :
	//                                ee);
	// }

	return energy
}

func vrna_eval_int_loop(fc *vrna_fold_compound_t, i, j, k, l int) int {
	e := INF

	e = eval_int_loop(fc, i, j, k, l)

	return e
}

func eval_int_loop(fc *vrna_fold_compound_t, i, j, k, l int) int {
	// var n_seq, s uint32
	var sn, ss []uint32
	// var a2s [][]uint32
	var e int
	var type_1, type2 uint32
	// with_ud
	var rtype [8]int
	var S, S2 []int
	// var SS, S5, S3 [][]int
	var P *vrna_param_t
	var md *vrna_md_t
	// var domains_up *vrna_ud_t
	// struct sc_int_dat sc_wrapper;

	// n_seq = 1
	P = fc.params
	md = P.model_details
	sn = fc.strand_number
	ss = fc.strand_start
	rtype = md.rtype
	S = fc.sequence_encoding
	S2 = fc.sequence_encoding2
	// SS = nil
	// S5 = nil
	// S3 = nil
	// a2s = nil
	// domains_up  = fc.domains_up;
	// with_ud     = ((domains_up) && (domains_up.energy_cb)) ? 1 : 0;
	e = INF

	// init_sc_int(fc, &sc_wrapper);

	{
		var energy, u1, u2 int
		// e5, e3,

		energy = 0

		type_1 = vrna_get_ptype_md(S2[i], S2[j], md)
		type2 = vrna_get_ptype_md(S2[l], S2[k], md)

		u1 = k - i - 1
		u2 = j - l - 1

		if (sn[i] == sn[k]) && (sn[l] == sn[j]) {
			/* regular interior loop */
			energy = E_IntLoop(u1, u2, type_1, type2, S[i+1], S[j-1], S[k-1], S[l+1], P)
		} else {
			/* interior loop like cofold structure */
			var Si, Sj int
			if sn[i+1] == sn[i] {
				Si = S[i+1]
			} else {
				Si = -1
			}

			if sn[j] == sn[j-1] {
				Sj = S[j-1]
			} else {
				Sj = -1
			}

			energy = E_IntLoop_Co(rtype[type_1], rtype[type2],
				i, j, k, l,
				int(ss[fc.strand_order[1]]), /* serves as cut point substitute */
				Si, Sj,
				S[k-1], S[l+1],
				md.dangles,
				P)
		}

		/* add soft constraints */
		// vivek: do something about sc
		// if (sc_wrapper.pair)
		//   energy += sc_wrapper.pair(i, j, k, l, &sc_wrapper);

		e = energy

		// vivek: do something about sc
		//   if with_ud {
		//     e5, e3 = 0, 0

		//     u1  = k - i - 1
		//     u2  = j - l - 1

		//     if u1 > 0 {
		//       e5 = domains_up.energy_cb(fc,
		//                                  i + 1, k - 1,
		//                                  VRNA_UNSTRUCTURED_DOMAIN_INT_LOOP,
		//                                  domains_up.data);
		//     }

		//     if (u2 > 0) {
		//       e3 = domains_up.energy_cb(fc,
		//                                  l + 1, j - 1,
		//                                  VRNA_UNSTRUCTURED_DOMAIN_INT_LOOP,
		//                                  domains_up.data);
		//     }

		//     e = Min(e, energy + e5);
		//     e = Min(e, energy + e3);
		//     e = Min(e, energy + e5 + e3);
		//   }
		// }

		// free_sc_int(&sc_wrapper);
	}

	return e
}

func vrna_get_ptype_md(i, j int, md *vrna_md_t) uint32 {
	var tt uint32 = uint32(md.pair[i][j])

	if tt == 0 {
		return 7
	} else {
		return tt
	}
}

func E_IntLoop(n1, n2 int, type_1, type_2 uint32, si1, sj1, sp1, sq1 int, P *vrna_param_t) int {
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
				energy += Min(MAX_NINIO, (nl-ns)*P.ninio[2])
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

			energy += Min(MAX_NINIO, (nl-ns)*P.ninio[2])

			energy += P.mismatchI[type_1][si1][sj1] + P.mismatchI[type_2][sq1][sp1]
		}
	}

	return energy
}

func ON_SAME_STRAND(I, J, C int) bool {
	if (I >= C) || (J < C) {
		return true
	} else {
		return false
	}
}

func E_IntLoop_Co(type_1, type_2, i, j, p, q, cutpoint, si1, sj1, sp1, sq1, dangles int,
	P *vrna_param_t) int {
	var e, energy, d3, d5, d5_2, d3_2, tmm, tmm_2 int
	var ci, cj, cp, cq bool

	energy = 0
	if type_1 > 2 {
		energy += P.TerminalAU
	}

	if type_2 > 2 {
		energy += P.TerminalAU
	}

	if dangles == 1 {
		return energy
	}

	ci = ON_SAME_STRAND(i, i+1, cutpoint)
	cj = ON_SAME_STRAND(j-1, j, cutpoint)
	cp = ON_SAME_STRAND(p-1, p, cutpoint)
	cq = ON_SAME_STRAND(q, q+1, cutpoint)

	if ci {
		d3 = P.dangle3[type_1][si1]
	}

	if cj {
		d5 = P.dangle5[type_1][sj1]
	}

	if cp {
		d5_2 = P.dangle5[type_2][sp1]
	}

	if cq {
		d3_2 = P.dangle3[type_2][sq1]
	}

	if cj && ci {
		tmm = P.mismatchExt[type_1][sj1][si1]
	} else {
		tmm = d5 + d3
	}

	if cp && cq {
		tmm_2 = P.mismatchExt[type_2][sp1][sq1]
	} else {
		tmm_2 = d5_2 + d3_2
	}

	if dangles == 2 {
		return energy + tmm + tmm_2
	}

	/* now we may have non-double dangles only */
	if p-i > 2 {
		if j-q > 2 {
			/* all degrees of freedom */
			e = Min(tmm, d5)
			e = Min(e, d3)
			energy += e
			e = Min(tmm_2, d5_2)
			e = Min(e, d3_2)
			energy += e
		} else if j-q == 2 {
			/* all degrees of freedom in 5' part between i and p */
			e = Min(tmm+d5_2, d3+d5_2)
			e = Min(e, d5+d5_2)
			e = Min(e, d3+tmm_2)
			e = Min(e, d3+d3_2)
			e = Min(e, tmm_2) /* no dangles on enclosing pair */
			e = Min(e, d5_2)  /* no dangles on enclosing pair */
			e = Min(e, d3_2)  /* no dangles on enclosing pair */
			energy += e
		} else {
			/* no unpaired base between q and j */
			energy += d3 + d5_2
		}
	} else if p-i == 2 {
		if j-q > 2 {
			/* all degrees of freedom in 3' part between q and j */
			e = Min(tmm+d3_2, d5+d3_2)
			e = Min(e, d5+d3_2)
			e = Min(e, d3+d3_2)
			e = Min(e, d5+tmm_2)
			e = Min(e, tmm_2)
			e = Min(e, d5_2)
			e = Min(e, d3_2)
			energy += e
		} else if j-q == 2 {
			/* one possible dangling base between either side */
			e = Min(tmm, tmm_2)
			e = Min(e, d3)
			e = Min(e, d5)
			e = Min(e, d5_2)
			e = Min(e, d3_2)
			e = Min(e, d3+d3_2)
			e = Min(e, d5+d5_2)
			energy += e
		} else {
			/* one unpaired base between i and p */
			energy += Min(d3, d5_2)
		}
	} else {
		/* no unpaired base between i and p */
		if j-q > 2 {
			/* all degrees of freedom in 3' part between q and j */
			energy += d5 + d3_2
		} else if j-q == 2 {
			/* one unpaired base between q and j */
			energy += Min(d5, d3_2)
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
func vrna_eval_hp_loop(fc *vrna_fold_compound_t, i, j int) int {
	// char * *Ss
	// var Ss []string
	// var a2s [][]uint32
	var S, S2 []int
	// var SS, S5, S3 [][]int
	var sn []uint32
	var u, e, noGUclosure int
	var type_1 uint32
	// n_seq, en, s
	var P *vrna_param_t
	var md *vrna_md_t
	// //  vrna_ud_t         *domains_up;
	//  struct sc_hp_dat  sc_wrapper;

	P = fc.params
	md = P.model_details
	noGUclosure = md.noGUclosure
	sn = fc.strand_number
	//  domains_up  = fc.domains_up
	e = INF

	if sn[j] != sn[i] {
		return eval_hp_loop_fake(fc, i, j)
	}

	//  vivek: handle sc stuff
	//  init_sc_hp(fc, &sc_wrapper);

	/* regular hairpin loop */

	S = fc.sequence_encoding
	S2 = fc.sequence_encoding2
	u = j - i - 1
	type_1 = vrna_get_ptype_md(S2[i], S2[j], md)

	if !(noGUclosure == 1 && ((type_1 == 3) || (type_1 == 4))) {
		e = E_Hairpin(u, type_1, S[i+1], S[j-1], fc.sequence[i-1:], P)
	}

	//  if e != INF {
	// 	 if (sc_wrapper.pair)
	// 		 e += sc_wrapper.pair(i, j, &sc_wrapper);

	// 	 /* consider possible ligand binding */
	// 	 if (domains_up && domains_up->energy_cb) {
	// 		 en = domains_up->energy_cb(fc,
	// 																i + 1, j - 1,
	// 																VRNA_UNSTRUCTURED_DOMAIN_HP_LOOP,
	// 																domains_up->data);
	// 		 if (en != INF)
	// 			 en += e;

	// 		 e = MIN2(e, en);
	// 	 }
	//  }

	//  free_sc_hp(&sc_wrapper);

	return e
}

func eval_hp_loop_fake(fc *vrna_fold_compound_t, i, j int) int {
	var S, S2 []int
	var sn []uint32
	// var u, e, ij, type_1, en, noGUclosure int
	var e, noGUclosure int
	var type_1 uint32
	// var idx []int
	var P *vrna_param_t
	// var sc *vrna_sc_t
	var md *vrna_md_t
	// vrna_ud_t     *domains_up;

	// idx = fc.jindx
	P = fc.params
	md = P.model_details
	noGUclosure = md.noGUclosure
	sn = fc.strand_number
	// domains_up  = fc.domains_up
	e = INF

	S = fc.sequence_encoding
	S2 = fc.sequence_encoding2
	// sc = fc.sc
	// u = j - i - 1
	// ij = idx[j] + i
	type_1 = vrna_get_ptype_md(S2[j], S2[i], md)

	if noGUclosure == 1 && (type_1 == 3 || type_1 == 4) {
		return e
	}

	/* hairpin-like exterior loop (for cofolding) */
	var si, sj int

	if sn[i+1] == sn[i] {
		si = S[i+1]
	} else {
		si = -1
	}

	if sn[j] == sn[j-1] {
		sj = S[j-1]
	} else {
		sj = -1
	}

	switch md.dangles {
	case 0:
		e = vrna_E_ext_stem(type_1, -1, -1, P)

	case 2:
		e = vrna_E_ext_stem(type_1, sj, si, P)

	default:
		e = vrna_E_ext_stem(type_1, -1, -1, P)
		e = Min(e, vrna_E_ext_stem(type_1, sj, -1, P))
		e = Min(e, vrna_E_ext_stem(type_1, -1, si, P))
		e = Min(e, vrna_E_ext_stem(type_1, sj, si, P))
	}

	/* add soft constraints */
	// vivek: handle sc stuff
	// if (sc) {
	// 	if sc->energy_up {
	// 		e += sc->energy_up[i + 1][u];
	// 	}

	// 	if (sc->energy_bp)
	// 		e += sc->energy_bp[ij];

	// 	if (sc->f)
	// 		e += sc->f(i, j, i, j, VRNA_DECOMP_PAIR_HP, sc->data);
	// }

	/* consider possible ligand binding */
	// if (domains_up && domains_up->energy_cb) {
	// 	en = domains_up->energy_cb(fc,
	// 															i + 1, j - 1,
	// 															VRNA_UNSTRUCTURED_DOMAIN_HP_LOOP,
	// 															domains_up->data);
	// 	if (en != INF)
	// 		en += e;

	// 	e = MIN2(e, en);
	// }

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
func E_Hairpin(size int, type_1 uint32, si1, sj1 int, sequence string, P *vrna_param_t) int {
	var energy int

	if size <= 30 {
		energy = P.hairpin[size]
	} else {
		energy = P.hairpin[30] + int(P.lxc*math.Log(float64(size)/30.0))
	}

	if size < 3 {
		return energy /* should only be the case when folding alignments */
	}

	if P.model_details.special_hp == 1 {
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
	}

	energy += P.mismatchH[type_1][si1][sj1]

	return energy
}

func cut_in_loop(i int, pt []int, sn []uint32) int {
	/* walk around the loop;  return 5' pos of first pair after cut if
	 * cut_point in loop else 0 */
	var p, j int

	p = pt[i]
	j = pt[i]

	for {
		i = pt[p]
		p = i + 1
		for pt[p] == 0 {
			p++
		}

		if !((p != j) && (sn[i] == sn[p])) {
			break
		}
	}
	if sn[i] == sn[p] {
		return 0
	} else {
		return p
	}
}

/**
 *** i is the 5'-base of the closing pair
 ***
 *** since each helix can coaxially stack with at most one of its
 *** neighbors we need an auxiliarry variable  cx_energy
 *** which contains the best energy given that the last two pairs stack.
 *** energy  holds the best energy given the previous two pairs do not
 *** stack (i.e. the two current helices may stack)
 *** We don't allow the last helix to stack with the first, thus we have to
 *** walk around the Loop twice with two starting points and take the minimum
 ***/
func energy_of_ml_pt(vc *vrna_fold_compound_t, i int, pt []int) (int, error) {
	var sn []uint32
	var energy, cx_energy, tmp, tmp2, best_energy, dangle_model, logML int
	var n, tt uint32
	// circular, ss, n, n_seq
	var idx []int
	var rtype [8]int
	best_energy = INF

	var i1, j, p, q, q_prev, q_prev2, u, x, count, mm5, mm3, ld5, new_cx,
		dang5, dang3, dang int
	var type_1 uint32
	// uu
	var e_stem, e_stem5, e_stem3, e_stem53 int
	var mlintern [NBPAIRS + 1]int
	var s, s1 []int
	// var S, S5, S3 [][]int
	// var a2s [][]uint32
	var P *vrna_param_t
	var md *vrna_md_t
	var sc *vrna_sc_t
	// var scs []*vrna_sc_t

	/* helper variables for dangles == 1|5 case */
	var E_mm5_available int  /* energy of 5' part where 5' mismatch of current stem is available */
	var E_mm5_occupied int   /* energy of 5' part where 5' mismatch of current stem is unavailable */
	var E2_mm5_available int /* energy of 5' part where 5' mismatch of current stem is available with possible 3' dangle for enclosing pair (i,j) */
	var E2_mm5_occupied int  /* energy of 5' part where 5' mismatch of current stem is unavailable with possible 3' dangle for enclosing pair (i,j) */

	n = vc.length
	sn = vc.strand_number
	P = vc.params
	md = P.model_details
	idx = vc.jindx

	// circular = md.circ
	dangle_model = md.dangles
	logML = md.logML
	rtype = md.rtype
	s = vc.sequence_encoding2
	sc = vc.sc
	s1 = vc.sequence_encoding
	// S = nil
	// S5 = nil
	// S3 = nil
	// a2s = nil
	// n_seq = 1
	// scs = nil

	bonus := 0

	if i >= pt[i] {
		return INF, errors.New("energy_of_ml_pt: i is not 5' base of a closing pair!")
	}

	if i == 0 {
		j = int(n) + 1
	} else {
		j = int(pt[i])
	}

	if i != 0 {
		/* (i,j) is closing pair of multibranch loop, add soft constraints */
		if sc != nil {
			if sc.energy_bp != nil {
				bonus += sc.energy_bp[idx[j]+i]
			}
		}
	}

	/* init the variables */
	energy = 0
	u = 0 /* the total number of unpaired nucleotides */
	p = i + 1
	q_prev = i - 1
	q_prev2 = i

	for x = 0; x <= NBPAIRS; x++ {
		mlintern[x] = P.MLintern[x]
	}

	/* seek to opening base of first stem */
	for p <= j && pt[p] == 0 {
		p++
	}

	/* add bonus energies for first stretch of unpaired nucleotides */

	u += p - i - 1
	if sc != nil {
		if sc.energy_up != nil {
			bonus += sc.energy_up[i+1][u]
		}
	}

	switch dangle_model {
	case 0:
		for p < j {
			/* p must have a pairing partner */
			q = int(pt[p])
			/* get type of base pair (p,q) */
			tt = vrna_get_ptype_md(s[p], s[q], md)

			energy += E_MLstem(tt, -1, -1, P)

			/* seek to the next stem */
			p = q + 1
			q_prev = q
			q_prev2 = q
			for p < j && pt[p] == 0 {
				p++
			}
			u += p - q - 1 /* add unpaired nucleotides */

			if sc != nil {
				if sc.energy_up != nil {
					bonus += sc.energy_up[q+1][p-q-1]
				}
			}
		}

		/* now lets get the energy of the enclosing stem */
		if i > 0 {
			/* actual closing pair */
			tt = vrna_get_ptype_md(s[j], s[i], md)

			energy += E_MLstem(tt, -1, -1, P)
		} else {
			/* virtual closing pair */
			energy += E_MLstem(0, -1, -1, P)
		}

	case 2:
		for p < j {
			/* p must have a pairing partner */
			q = int(pt[p])
			/* get type of base pair (p,q) */
			tt = vrna_get_ptype_md(s[p], s[q], md)

			if sn[p-1] == sn[p] {
				mm5 = s1[p-1]
			} else {
				mm5 = -1
			}

			if sn[q] == sn[q+1] {
				mm3 = s1[q+1]
			} else {
				mm3 = -1
			}

			energy += E_MLstem(tt, mm5, mm3, P)

			/* seek to the next stem */
			p = q + 1
			q_prev = q
			q_prev2 = q
			for p < j && pt[p] == 0 {
				p++
			}
			u += p - q - 1 /* add unpaired nucleotides */

			if sc != nil {
				if sc.energy_up != nil {
					bonus += sc.energy_up[q+1][p-q-1]
				}
			}
		}
		if i > 0 {
			/* actual closing pair */
			tt = vrna_get_ptype_md(s[j], s[i], md)

			if sn[j-1] == sn[j] {
				mm5 = s1[j-1]
			} else {
				mm5 = -1
			}

			if sn[i] == uint32(s1[i+1]) {
				mm3 = s1[i+1]
			} else {
				mm3 = -1
			}

			energy += E_MLstem(tt, mm5, mm3, P)
		} else {
			/* virtual closing pair */
			energy += E_MLstem(0, -1, -1, P)
		}

	case 3: /* we treat helix stacking different */
		for count = 0; count < 2; count++ {
			/* do it twice */
			ld5 = 0 /* 5' dangle energy on prev pair (type) */
			if i == 0 {
				j = pt[0] + 1
				type_1 = 0 /* no pair */
			} else {
				j = pt[i]
				type_1 = vrna_get_ptype_md(s[j], s[i], md)

				/* prime the ld5 variable */
				if sn[j-1] == sn[j] {
					ld5 = P.dangle5[type_1][s1[j-1]]
					p = pt[j-2]

					if p == 1 && sn[j-2] == sn[j-1] {
						if P.dangle3[md.pair[s[p]][s[j-2]]][s1[j-1]] < ld5 {
							ld5 = 0
						}
					}
				}
			}

			i1 = i
			p = i + 1
			u = 0
			energy = 0
			cx_energy = INF

			for {
				/* walk around the multi-loop */
				new_cx = INF

				/* hop over unpaired positions */
				for p <= pt[0] && pt[p] == 0 {
					p++
				}

				/* memorize number of unpaired positions */
				u += p - i1 - 1

				if sc != nil {
					if sc.energy_up != nil {
						bonus += sc.energy_up[i1+1][p-i1-1]
					}
				}

				/* get position of pairing partner */
				if p == pt[0]+1 {
					q = 0
					tt = 0 /* virtual root pair */
				} else {
					q = pt[p]
					/* get type of base pair P->q */
					tt = vrna_get_ptype_md(s[p], s[q], md)
				}

				energy += mlintern[tt]
				cx_energy += mlintern[tt]

				dang5 = 0
				dang3 = 0
				if (sn[p-1] == sn[p]) && p > 1 {
					dang5 = P.dangle5[tt][s1[p-1]] /* 5'dangle of pq pair */
				}

				if (sn[i1] == sn[i1+1]) && i1 < s[0] {
					dang3 = P.dangle3[type_1][s1[i1+1]] /* 3'dangle of previous pair */
				}

				switch p - i1 - 1 {
				case 0:
					/* adjacent helices */
					if i1 != 0 {
						if sn[i1] == sn[p] {
							new_cx = energy + P.stack[rtype[type_1]][rtype[tt]]
							/* subtract 5'dangle and TerminalAU penalty */
							new_cx += -ld5 - mlintern[tt] - mlintern[type_1] + 2*mlintern[1]
						}

						ld5 = 0
						energy = Min(energy, cx_energy)
					}

				case 1: /* 1 unpaired base between helices */
					dang = Min(dang3, dang5)
					energy = energy + dang
					ld5 = dang - dang3
					/* may be problem here: Suppose
					* cx_energy>energy, cx_energy+dang5<energy
					* and the following helices are also stacked (i.e.
					* we'll subtract the dang5 again */
					if cx_energy+dang5 < energy {
						energy = cx_energy + dang5
						ld5 = dang5
					}

					new_cx = INF /* no coax stacking with mismatch for now */
				default: /* many unpaired base between helices */
					energy += dang5 + dang3
					energy = Min(energy, cx_energy+dang5)
					new_cx = INF /* no coax stacking possible */
					ld5 = dang5
				}
				type_1 = tt
				cx_energy = new_cx
				i1 = q
				p = q + 1

				if q == i {
					break
				}
			}
			best_energy = Min(energy, best_energy) /* don't use cx_energy here */
			/* skip a helix and start again */
			for pt[p] == 0 {
				p++
			}
			if i == pt[p] {
				break
			}

			i = pt[p]
		} /* end doing it twice */

		energy = best_energy

	default:
		E_mm5_available = INF
		E2_mm5_available = INF
		E_mm5_occupied = 0
		E2_mm5_occupied = 0
		for p < j {
			/* p must have a pairing partner */
			q = int(pt[p])
			/* get type of base pair (p,q) */
			tt = vrna_get_ptype_md(s[p], s[q], md)

			if q_prev+2 < p {
				E_mm5_available = Min(E_mm5_available, E_mm5_occupied)
				E_mm5_occupied = E_mm5_available
			}

			if q_prev2+2 < p {
				E2_mm5_available = Min(E2_mm5_available, E2_mm5_occupied)
				E2_mm5_occupied = E2_mm5_available
			}

			if (sn[p-1] == sn[p]) && pt[p-1] == 0 {
				mm5 = s1[p-1]
			} else {
				mm5 = -1
			}

			if (sn[q] == sn[q+1]) && pt[q+1] == 0 {
				mm3 = s1[q+1]
			} else {
				mm3 = -1
			}

			e_stem = E_MLstem(tt, -1, -1, P)
			e_stem5 = E_MLstem(tt, mm5, -1, P)
			e_stem3 = E_MLstem(tt, -1, mm3, P)
			e_stem53 = E_MLstem(tt, mm5, mm3, P)

			tmp = E_mm5_occupied + e_stem3
			tmp = Min(tmp, E_mm5_available+e_stem53)
			tmp = Min(tmp, E_mm5_available+e_stem3)
			tmp2 = E_mm5_occupied + e_stem
			tmp2 = Min(tmp2, E_mm5_available+e_stem5)
			tmp2 = Min(tmp2, E_mm5_available+e_stem)

			E_mm5_occupied = tmp
			E_mm5_available = tmp2

			tmp = E2_mm5_occupied + e_stem3
			tmp = Min(tmp, E2_mm5_available+e_stem53)
			tmp = Min(tmp, E2_mm5_available+e_stem3)
			tmp2 = E2_mm5_occupied + e_stem
			tmp2 = Min(tmp2, E2_mm5_available+e_stem5)
			tmp2 = Min(tmp2, E2_mm5_available+e_stem)

			E2_mm5_occupied = tmp
			E2_mm5_available = tmp2

			/* seek to the next stem */
			p = q + 1
			q_prev = q
			q_prev2 = q
			for p < j && pt[p] == 0 {
				p++
			}
			u += p - q - 1 /* add unpaired nucleotides */

			if sc != nil {
				if sc.energy_up != nil {
					bonus += sc.energy_up[q+1][p-q-1]
				}
			}
		}

		if i > 0 {
			/* actual closing pair */
			type_1 = vrna_get_ptype_md(s[j], s[i], md)

			if (sn[j-1] == sn[j]) && pt[j-1] == 0 {
				mm5 = s1[j-1]
			} else {
				mm5 = -1
			}

			if (sn[i] == sn[i+1]) && pt[i+1] == 0 {
				mm3 = s1[i+1]
			} else {
				mm3 = -1
			}

			if q_prev+2 < p {
				E_mm5_available = Min(E_mm5_available, E_mm5_occupied)
				E_mm5_occupied = E_mm5_available
			}

			if q_prev2+2 < p {
				E2_mm5_available = Min(E2_mm5_available, E2_mm5_occupied)
				E2_mm5_occupied = E2_mm5_available
			}

			e_stem = E_MLstem(type_1, -1, -1, P)
			e_stem5 = E_MLstem(type_1, mm5, -1, P)
			e_stem3 = E_MLstem(type_1, -1, mm3, P)
			e_stem53 = E_MLstem(type_1, mm5, mm3, P)
		} else {
			/* virtual closing pair */
			e_stem = E_MLstem(0, -1, -1, P)
			e_stem5 = e_stem
			e_stem3 = e_stem
			e_stem53 = e_stem
		}

		/* now lets see how we get the minimum including the enclosing stem */
		energy = E_mm5_occupied + e_stem
		energy = Min(energy, E_mm5_available+e_stem5)
		energy = Min(energy, E_mm5_available+e_stem)
		energy = Min(energy, E2_mm5_occupied+e_stem3)
		energy = Min(energy, E2_mm5_occupied+e_stem)
		energy = Min(energy, E2_mm5_available+e_stem53)
		energy = Min(energy, E2_mm5_available+e_stem3)
		energy = Min(energy, E2_mm5_available+e_stem5)
		energy = Min(energy, E2_mm5_available+e_stem)
	} /* end switch dangle_model */

	energy += P.MLclosing

	/*
	* logarithmic ML loop energy if logML
	* does this work for comparative predictions as well?
	 */
	if logML == 1 && (u > 6) {
		energy += 6*P.MLbase + int(P.lxc*math.Log(float64(u)/6.0))
	} else {
		energy += u * P.MLbase
	}

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
func E_MLstem(type_1 uint32, si1, sj1 int, P *vrna_param_t) int {
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
