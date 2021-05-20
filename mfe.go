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
 *  @warning  Reading/Writing from/to attributes that are not within the scope of the current type usually result
 *  in undefined behavior!
 *
 *  @see  #vrna_fold_compound_t.type, vrna_fold_compound(), vrna_fold_compound_comparative(), vrna_fold_compound_free(),
 *        #VRNA_FC_TYPE_SINGLE, #VRNA_FC_TYPE_COMPARATIVE
 */
type vrna_fold_compound_t struct {
	length uint32 /**<  @brief  The length of the sequence (or sequence alignment). vivek: In ViennaRNA, length is stored as unsigned int (a primitive C type which ranges from 0 to 4,294,967,295). uint32 is the equivalent type in Go. */
	// cutpoint int    /*  @brief  The position of the (cofold) cutpoint within the provided sequence. If there is no cutpoint, this field will be set to -1 */

	/*
	* @brief  The strand number a particular nucleotide is associated with
	* vivek:
	* size of this array is `len(sequence) + 2`
	* [0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0]
	 */
	strand_number []uint32
	/*
	* @brief  The strand order, i.e. permutation of current concatenated sequence
	* vivek:
	* size is `strands` + 1 (The list is zero-indexed)
	* [0 0]
	 */
	strand_order []uint32

	/*
	* @brief  The start position of a particular strand within the current concatenated sequence
	* vivek:
	* size is `strands` + 1 (The list is zero-indexed)
	* [1 0]
	 */
	strand_start []uint32

	/*
	* @brief  The end (last) position of a particular strand within the current concatenated sequence
	* vivek:
	* size is `strands` + 1 (The list is zero-indexed)
	* [30 0]
	 */
	strand_end []uint32

	/* vivek: The number of strands in this folding compound. Since we only calculate MFE for one compound, this is generally one.*/
	strands    uint32
	seq_struct *vrna_seq_t
	// nucleotides []vrna_seq_t /* */
	// alignment   *vrna_msa_t

	// hc *vrna_hc_t /**<  @brief  The hard constraints data structure used for structure prediction */

	// matrices *vrna_mx_mfe_t /**<  @brief  The MFE DP matrices */
	// exp_matrices *vrna_mx_pf_t  /**<  @brief  The PF DP matrices  */

	params *vrna_param_t /**<  @brief  The precomputed free energy contributions for each type of loop */
	// exp_params *vrna_exp_param_t    /**<  @brief  The precomputed free energy contributions as Boltzmann factors  */

	// iindx []int /**<  @brief  DP matrix accessor  */
	// jindx []int /**<  @brief  DP matrix accessor  */

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
	sequence          string /* The input sequence string */
	sequence_encoding []int  /**<  @brief  Numerical encoding of the input sequence
	 *    @see    vrna_sequence_encode()
	 */
	sequence_encoding2 []int
	// pair_type              string
	/**<  @brief  Pair type array
	 *
	 *    Contains the numerical encoding of the pair type for each pair (i,j) used
	 *    in MFE, Partition function and Evaluation computations.
	 *    @note This array is always indexed via jindx, in contrast to previously
	 *    different indexing between mfe and pf variants!
	 *    @warning   Only available if @verbatim type==VRNA_FC_TYPE_SINGLE @endverbatim
	 *    @see    vrna_idx_col_wise(), vrna_ptypes()
	 */

	// sc *vrna_sc_t
	/**<  @brief  The soft constraints for usage in structure prediction and evaluation
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
	// maxD1         uint32 /**<  @brief  Maximum allowed base pair distance to first reference */
	// maxD2         uint32 /**<  @brief  Maximum allowed base pair distance to second reference */
	// reference_pt1 []int  /**<  @brief  A pairtable of the first reference structure */
	// reference_pt2 []int  /**<  @brief  A pairtable of the second reference structure */

	// referenceBPs1 []uint32 /**<  @brief  Matrix containing number of basepairs of reference structure1 in interval [i,j] */
	// referenceBPs2 []uint32 /**<  @brief  Matrix containing number of basepairs of reference structure2 in interval [i,j] */
	// bpdist        []uint32 /**<  @brief  Matrix containing base pair distance of reference structure 1 and 2 on interval [i,j] */

	// mm1 []uint32 /**<  @brief  Maximum matching matrix, reference struct 1 disallowed */
	// mm2 []uint32 /**<  @brief  Maximum matching matrix, reference struct 2 disallowed */

	/**
	 *  @}
	 */

	/**
	 *  @name Additional data fields for local folding
	 *
	 *  These data fields are typically populated with meaningful data only if used in the context of local folding
	 *  @{
	 */
	window_size int /**<  @brief  window size for local folding sliding window approach */
	// ptype_local []string /**<  @brief  Pair type array (for local folding) */

	/**
	 *  @}
	 */
}

/* vivek:
* Returns a vrna_fold_compound_t with default options set for its variables.
* See the doc of `vrna_fold_compound_t` for explaination of the variables of the
* struct.
 */
func default_vrna_fold_compound_t() *vrna_fold_compound_t {
	return &vrna_fold_compound_t{
		length:  0,
		strands: 0,
		// cutpoint:      -1,
		strand_number: nil,
		strand_order:  nil,
		strand_start:  nil,
		strand_end:    nil,
		seq_struct:    nil,
		// nucleotides:  nil,
		// alignment:     nil,

		// hc:       nil,
		// matrices: nil,
		// exp_matrices: nil,
		params: nil,
		// exp_params:   nil,
		// iindx: nil,
		// jindx: nil,

		// stat_cb:      nil,
		// auxdata:      nil,
		// free_auxdata: nil,

		// domains_struc: nil,
		// domains_up:    nil,
		// aux_grammar:   nil,

		sequence:           "",
		sequence_encoding:  make([]int, 0),
		sequence_encoding2: nil,
		// ptype:              "",
		// sc: nil,

		// axD1:          0,
		// axD2:          0,
		// reference_pt1: nil,
		// reference_pt2: nil,
		// referenceBPs1: nil,
		// referenceBPs2: nil,
		// bpdist:        nil,
		// mm1:           nil,
		// mm2:           nil,

		window_size: -1,
		// ptype_local: nil,
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
// type vrna_sc_bp_storage_t struct {
// 	interval_start, interval_end uint32
// 	e                            int
// }

/**
 *  @brief  The soft constraints data structure
 *
 *  @ingroup soft_constraints
 */
// type vrna_sc_t struct {
// 	n uint32

// 	state rune

// 	energy_up     [][]int     /**<  @brief Energy contribution for stretches of unpaired nucleotides */
// 	exp_energy_up [][]float64 /**<  @brief Boltzmann Factors of the energy contributions for unpaired sequence stretches */

// 	up_storage []int /**<  @brief  Storage container for energy contributions per unpaired nucleotide */
// 	// bp_storage [][]vrna_sc_bp_storage_t /**<  @brief  Storage container for energy contributions per base pair */

// 	energy_bp     []int     /**<  @brief Energy contribution for base pairs */
// 	exp_energy_bp []float64 /**<  @brief Boltzmann Factors of the energy contribution for base pairs */

// 	energy_bp_local     [][]int     /**<  @brief Energy contribution for base pairs (sliding window approach) */
// 	exp_energy_bp_local [][]float64 /**<  @brief Boltzmann Factors of the energy contribution for base pairs (sliding window approach) */

// 	energy_stack     []int     /**<  @brief Pseudo Energy contribution per base pair involved in a stack */
// 	exp_energy_stack []float64 /**<  @brief Boltzmann weighted pseudo energy contribution per nucleotide involved in a stack */

// 	// /* generic soft contraints below */
// 	// vrna_callback_sc_energy     *f;     /**<  @brief  A function pointer used for pseudo
// 	//                                      *            energy contribution in MFE calculations
// 	//                                      *    @see    vrna_sc_add_f()
// 	//                                      */

// 	// vrna_callback_sc_backtrack  *bt;    /**<  @brief  A function pointer used to obtain backtraced
// 	//                                      *            base pairs in loop regions that were altered
// 	//                                      *            by soft constrained pseudo energy contributions
// 	//                                      *    @see    vrna_sc_add_bt()
// 	//                                      */

// 	// vrna_callback_sc_exp_energy *exp_f; /**<  @brief  A function pointer used for pseudo energy
// 	//                                      *            contribution boltzmann factors in PF
// 	//                                      *            calculations
// 	//                                      *    @see    vrna_sc_add_exp_f()
// 	//                                      */

// 	// void                        *data;  /**<  @brief  A pointer to the data object provided for
// 	//                                      *            for pseudo energy contribution functions of the
// 	//                                      *            generic soft constraints feature
// 	//                                      */
// 	// vrna_callback_free_auxdata  *free_data;
// }

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
	// log.Printf("%v", seq)
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
	// log.Printf("%v", fc)
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
	MAXALPHA int = 20 /** @brief Maximal length of alphabet */
	MAXLOOP  int = 30 /** The maximum loop length */

)

type vrna_md_t struct {
	temperature float64 /**<  @brief  The temperature used to scale the thermodynamic parameters */
	betaScale   float64 /**<  @brief  A scaling factor for the thermodynamic temperature of the Boltzmann factors */
	pf_smooth   int     /**<  @brief  A flat specifying whether energies in Boltzmann factors need to be smoothed */
	// dangles     int
	/**<  @brief  Specifies the dangle model used in any energy evaluation (0,1,2 or 3)
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
	// special_hp     int  /**<  @brief  Include special hairpin contributions for tri, tetra and hexaloops */
	// noLP           int  /**<  @brief  Only consider canonical structures, i.e. no 'lonely' base pairs */
	// noGU int /**<  @brief  Do not allow GU pairs */
	// noGUclosure    int  /**<  @brief  Do not allow loops to be closed by GU pair */
	// logML          int  /**<  @brief  Use logarithmic scaling for multiloops */
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
	window_size int                             /**<  @brief  Size of the sliding window for locally optimal structure prediction */
	oldAliEn    int                             /**<  @brief  Use old alifold energy model */
	ribo        int                             /**<  @brief  Use ribosum scoring table in alifold energy model */
	cv_fact     float64                         /**<  @brief  Co-variance scaling factor for consensus structure prediction */
	nc_fact     float64                         /**<  @brief  Scaling factor to weight co-variance contributions of non-canonical pairs */
	sfact       float64                         /**<  @brief  Scaling factor for partition function scaling */
	rtype       [8]int                          /**<  @brief  Reverse base pair type array */
	alias       [MAXALPHA + 1]int               /**<  @brief  alias of an integer nucleotide representation */
	pair        [MAXALPHA + 1][MAXALPHA + 1]int /**<  @brief  Integer representation of a base pair */
}

var (
	// Default temperature for structure prediction and free energy evaluation in $^\circ C$
	VRNA_MODEL_DEFAULT_TEMPERATURE float64 = 37.0
	VRNA_MODEL_DEFAULT_PF_SMOOTH   int     = 1
	// Default dangling end model
	// VRNA_MODEL_DEFAULT_DANGLES int = 2
	// Default model behavior for lookup of special tri-, tetra-, and hexa-loops
	// VRNA_MODEL_DEFAULT_SPECIAL_HP int = 1
	// Default model behavior for so-called 'lonely pairs'
	VRNA_MODEL_DEFAULT_NO_LP int = 1
	// Default model behavior for G-U base pairs
	// VRNA_MODEL_DEFAULT_NO_GU int = 0
	// Default model behavior for G-U base pairs closing a loop
	// VRNA_MODEL_DEFAULT_NO_GU_CLOSURE int = 0
	// Default model behavior on how to evaluate the energy contribution of multi-branch loops
	VRNA_MODEL_DEFAULT_LOG_ML int = 0
	// Default model behavior to treat a molecule as a circular RNA (DNA)
	VRNA_MODEL_DEFAULT_CIRC int = 0
	// Default model behavior regarding the treatment of G-Quadruplexes
	// VRNA_MODEL_DEFAULT_GQUAD int = 0
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
		temperature: VRNA_MODEL_DEFAULT_TEMPERATURE,
		betaScale:   1.,
		pf_smooth:   VRNA_MODEL_DEFAULT_PF_SMOOTH,
		// dangles:        VRNA_MODEL_DEFAULT_DANGLES,
		// special_hp:     VRNA_MODEL_DEFAULT_SPECIAL_HP,
		// noLP:           VRNA_MODEL_DEFAULT_NO_LP,
		// noGU: VRNA_MODEL_DEFAULT_NO_GU,
		// noGUclosure:    VRNA_MODEL_DEFAULT_NO_GU_CLOSURE,
		// logML:          VRNA_MODEL_DEFAULT_LOG_ML,
		circ: VRNA_MODEL_DEFAULT_CIRC,
		// gquad:          VRNA_MODEL_DEFAULT_GQUAD,
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
		alias: [MAXALPHA + 1]int{
			//   _  A  C  G  U  X  K  I
			/**/ 0, 1, 2, 3, 4, 3, 2, 0,
		},
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
				// if model_details.dangles > 0 {
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

var SHRT_MAX uint32 = 32767

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

	if uint32(len_seq) > SHRT_MAX {
		err := fmt.Sprintf("vrna_pair_table_from_string: Structure too long to be converted to pair table (n=%v, max=%v)", len_seq, SHRT_MAX)
		return nil, errors.New(err)
	}

	pair_table = make([]int, len_seq+2)
	pair_table[0] = len_seq

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
		if structure_char_slice[i] == open_bracket {
			// pop index of open bracket onto stack
			open_bracket_idx_stack[stack_idx] = i + 1
			stack_idx++
		} else if structure_char_slice[i] == close_bracket {
			stack_idx--

			if stack_idx < 0 {
				return nil, fmt.Errorf("%v\nunbalanced brackets '%v%v' found while extracting base pairs", structure,
					open_bracket, close_bracket)
			}

			open_bracket_idx := open_bracket_idx_stack[stack_idx]
			// current index of one-indexed sequence
			curr_idx := i + 1
			pair_table[curr_idx] = open_bracket_idx
			pair_table[open_bracket_idx] = curr_idx
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
	// var sn []uint32
	var i, length, energy int

	length = int(fc.length)

	energy = energy_of_extLoop_pair_table(fc, 0, pair_table)

	vrna_cstr_print_eval_ext_loop(energy)
	// vivek: assume input is only A, T, C, G, or U
	for i = 1; i <= length; i++ {
		if pair_table[i] == 0 {
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

func sanitize_bp_span(fc *vrna_fold_compound_t) {
	var md *vrna_md_t = fc.params.model_details

	/* non-local fold mode */
	md.window_size = int(fc.length)

	if md.max_bp_span <= 0 || md.max_bp_span > md.window_size {
		md.max_bp_span = md.window_size
	}
}

func set_fold_compound(fc *vrna_fold_compound_t) {
	// VRNA_OPTION_MFE | VRNA_OPTION_EVAL_ONLY
	var sequence string = fc.sequence
	// var sequences []string
	// var length, s int
	// var md_p *vrna_md_t

	// md_p = fc.params.model_details
	fc.sequence = ""
	fc.length = 0

	/* split input sequences at default delimiter '&' */
	// sequences = vrna_strsplit(sequence, NULL);

	vrna_sequence_add(fc, sequence)

	// if fc.strands > 1 {
	// 	fc.cutpoint = int(fc.nucleotides[0].length) + 1

	// 	if md_p.min_loop_size == TURN {
	// 		md_p.min_loop_size = 0 /* is it safe to set this here? */
	// 	}
	// }

	// if !(options & VRNA_OPTION_EVAL_ONLY) {
	// 	fc.ptype = (aux & WITH_PTYPE) ? vrna_ptypes(fc.sequence_encoding2, md_p) : NULL;
	// 	/* backward compatibility ptypes */
	// 	fc.ptype_pf_compat =
	// 			(aux & WITH_PTYPE_COMPAT) ? get_ptypes(fc.sequence_encoding2, md_p, 1) : NULL;
	// }

	vrna_sequence_prepare(fc)

	// fc.iindx = vrna_idx_row_wise(fc.length)
	// fc.jindx = vrna_idx_col_wise(fc.length)

}

/*
* Sets the `seq_struct`, `sequence_encoding` and `sequence_encoding2` fields of
* `fc`.
* `sequence_encoding` is the set to the same output as `vrna_seq_encode`
* `sequence_encoding2` is the same as `sequence_encoding`, except `sequence_encoding[0] = len(sequence)`
 */
func vrna_sequence_add(fc *vrna_fold_compound_t, sequence string) {
	// VRNA_SEQUENCE_RNA
	var add_length uint32 = uint32(len(sequence))

	/* add the sequence to the nucleotides container */
	// fc.nucleotides = (vrna_seq_t *)vrna_realloc(vc.nucleotides,
	// 																							sizeof(vrna_seq_t) *
	// 																							(vc.strands + 1));
	// fc.nucleotides = make([]vrna_seq_t, fc.strands+1)

	fc.seq_struct = &vrna_seq_t{
		sequence: strings.ToUpper(sequence),
		length:   uint32(len(sequence)),
		encoding: vrna_seq_encode(sequence),
	}
	// fc.seq_struct = set_sequence(sequence,
	// 	fc.params.model_details)

	/* increase strands counter */
	fc.strands++

	/* add new sequence to initial order of all strands */
	// fc.sequence = (char *)vrna_realloc(vc.sequence,
	// 																		sizeof(char) *
	// 																		(vc.length + add_length + 1));
	// vivek: should we add an ampersand between sequences?
	fc.sequence = sequence

	/* add encoding for new strand */
	// fc.sequence_encoding = make([]int, fc.length+add_length+2)

	// fc.sequence_encoding = fc.sequence_encoding + '\0' + fc.nucleotides[fc.strands - 1].encoding
	// log.Printf("struct encodinggg: %v", fc.seq_struct.encoding)
	fc.sequence_encoding = make([]int, len(sequence)+2)
	fc.sequence_encoding = append([]int{0}, fc.seq_struct.encoding[1:]...)
	// log.Printf("sequence encoding: %v", fc.sequence_encoding)
	fc.sequence_encoding[len(sequence)+1] = fc.sequence_encoding[1]
	fc.sequence_encoding[0] = fc.sequence_encoding[len(sequence)]
	// log.Printf("circular encoding: %v", fc.sequence_encoding)
	// fc.sequence_encoding[1:] = fc.seq_struct.encoding[1:]

	// fc.sequence_encoding = append(fc.sequence_encoding, fc.seq_struct.encoding[1:]...)

	/* restore circular encoding */
	// fc.sequence_encoding = append(fc.sequence_encoding, fc.sequence_encoding[1])
	// fmt.Println(len(fc.sequence_encoding))
	// fc.sequence_encoding[0] = fc.sequence_encoding[len(fc.sequence_encoding)-1]

	/* add encoding2 (simple encoding) for new strand */
	// fc.sequence_encoding2 = make([]int, fc.length+add_length+2)

	fc.sequence_encoding2 = make([]int, len(sequence)+2)
	enc := vrna_seq_encode_simple(fc.seq_struct.sequence)
	fc.sequence_encoding2 = append([]int{0}, enc[1:]...)
	fc.sequence_encoding2[len(sequence)+1] = fc.sequence_encoding2[1]
	fc.sequence_encoding2[0] = len(sequence)
	// log.Printf("circulaencoding2: %v", fc.sequence_encoding2)
	// fc.sequence_encoding2[1:] = enc[1:]
	// panic(fmt.Sprintf("%v", fc.sequence_encoding2))
	// fc.sequence_encoding2[fc.length+1:] = enc[1 : 1+add_length]
	// memcpy(vc.sequence_encoding2 + vc.length + 1,
	// 				enc + 1,
	// 				add_length * sizeof(short))

	// fc.sequence_encoding2[len(fc.sequence_encoding2)-1] = fc.sequence_encoding2[1]
	// fc.sequence_encoding2[0] = int(fc.length + add_length)

	/* finally, increase length property of the fold compound */
	fc.length = add_length
}

/* vivek:
* Setting of the `encoding` field is explained in the docs of `vrna_seq_encode`.
*
* The original repo has additional functionality for a circular sequence which has been omitted from this function.
 */
// func set_sequence(sequence string, md *vrna_md_t) *vrna_seq_t {
// 	return
// }

/* vivek:
* manipulates the encoded sequence of `sequence` (see doc of `vrna_seq_encode_simple` for more
* details of the encoding process) by setting `S[0]` to `S[len(sequence)]`
* If a sequence is axxx...xxxb where a and b are the first and last nucleotides
* of the sequence, and AXXX...XXXB is the encoded `sequence`, `vrna_seq_encode`
* returns BAXXX...XXXBA.
* thinking out loud:
* * could be done to make the sequence circular
* * may be needed in functions that require the previous and next values to
* compute a value
 */
func vrna_seq_encode(sequence string) []int {
	var l uint32
	var S []int

	if sequence != "" {
		// S[1:len(sequence)+1] contains the encoding of sequence
		// S[0] is the lenght of sequence
		// S[len(sequence)+1] = S[1]
		S = vrna_seq_encode_simple(sequence)

		l = uint32(len(sequence))
		// log.Printf("here")
		// log.Printf("%v", S)
		// for i = 1; i <= l; i++ {
		// 	/* _  A  C  G  U  X  K  I */
		// 	S[i] = md.alias[S[i]]
		// }
		// log.Printf("%v", S)

		S[l+1] = S[1]
		S[0] = S[l]
		// log.Printf("%v", S)
	}

	return S
}

// vivek:
// encodes `sequence` based on `vrna_nucleotide_encode` into `S[1:len(sequence)-1]`
// `S[0]` is the length of `sequence`
// `S[len(sequence) + 1]` (or `S[-1]`) == `S[1]` which makes the sequence circular
// not sure why the orig repo makes the sequence cirular

/* Contains the numerical encoding of the pair type for each pair (i,j) used
 *    @note This array is always indexed via jindx
 */
func vrna_seq_encode_simple(sequence string) []int {
	var i, l uint32
	var S []int

	if sequence != "" {
		l = uint32(len(sequence))
		S = make([]int, l+2)

		for i = 1; i <= l; i++ { /* make numerical encoding of sequence */
			S[i] = vrna_nucleotide_encode(([]rune(sequence))[i-1])
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

/* vivek:
* Maps:
*	A -> 1
* C -> 2
* G -> 3
* T / U -> 4
 */
func vrna_nucleotide_encode(c rune) int {
	/* return numerical representation of nucleotide used e.g. in vrna_md_t.pair[][] */ //
	var code int = -1

	c = unicode.ToUpper(c)

	// if md != nil {
	// vivek: md.energy_set defaults to 0
	// if md.energy_set > 0 {
	// 	code = int(c-'A') + 1
	// } else {
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
	// }
	// }

	return code
}

/* vivek:
* Sets the `strand_order`, `strand_start`, `strand_end` and `strand_number` fields
* of `fc`.
 */
// func vrna_sequence_prepare(fc *vrna_fold_compound_t) {
// 	var strand_number, i uint32

// 	fc.strand_order = nil
// 	fc.strand_start = nil
// 	fc.strand_end = nil
// 	fc.strand_number = make([]uint32, fc.length+2)

// 	/* 1. store initial strand order */
// 	fc.strand_order = make([]uint32, fc.strands+1)
// 	for strand_number = 0; strand_number < fc.strands; strand_number++ {
// 		fc.strand_order[strand_number] = strand_number
// 	}
// 	// log.Printf("strand_order: %v", fc.strand_order)

// 	/* 2. mark start and end positions of sequences */
// 	fc.strand_start = make([]uint32, fc.strands+1)
// 	fc.strand_end = make([]uint32, fc.strands+1)

// 	fc.strand_start[0] = 1
// 	fc.strand_end[0] = fc.strand_start[0] + fc.seq_struct.length - 1

// 	for strand_number = 1; strand_number < fc.strands; strand_number++ {
// 		fc.strand_start[strand_number] = fc.strand_end[strand_number-1] + 1
// 		fc.strand_end[strand_number] = fc.strand_start[strand_number] + fc.seq_struct.length - 1
// 		for i = fc.strand_start[strand_number]; i <= fc.strand_end[strand_number]; i++ {
// 			fc.strand_number[i] = strand_number
// 		}
// 	}

// 	/* this sets pos. n + 1 as well */
// 	fc.strand_number[fc.length+1] = fc.strands - 1
// }

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
	for pair_table_idx <= length && (pair_table[pair_table_idx] == 0) {
		pair_table_idx++
	}

	for pair_table_idx < length {
		/* pair_table_idx  must have a pairing partner */
		pair_close_idx := int(pair_table[pair_table_idx])

		/* get type of base pair (p,q) */
		pair_type := vrna_get_pair_type_md(fc.sequence_encoding2[pair_table_idx],
			fc.sequence_encoding2[pair_close_idx],
			fc.params.model_details)

		// vivek: cryptic variable names. 5 and 3 have to do with the 5' and 3' ends
		// of a DNA/RNA strand
		var mm5, mm3 int
		if pair_table_idx > 1 {
			mm5 = fc.sequence_encoding[pair_table_idx-1]
		} else {
			mm5 = -1
		}

		if pair_close_idx < length {
			mm3 = fc.sequence_encoding[pair_close_idx+1]
		} else {
			mm3 = -1
		}

		energy += vrna_E_ext_stem(pair_type, mm5, mm3, fc.params)

		/* seek to the next stem */
		pair_table_idx = pair_close_idx + 1
		for pair_table_idx <= length && pair_table[pair_table_idx] == 0 {
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
func vrna_E_ext_stem(pair_type uint32, n5d int, n3d int, p *vrna_param_t) int {
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

	if fc.params.model_details.pair[fc.sequence_encoding2[pair_open_idx]][fc.sequence_encoding2[pair_close_idx]] == 0 {
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
		for pair_table[n5p_iterator] == 0 {
			n5p_iterator++
		}

		// seek to closing bracket from 3' end
		n3p_iterator--
		for pair_table[n3p_iterator] == 0 {
			n3p_iterator--
		}

		if (pair_table[n3p_iterator] != n5p_iterator) || (n5p_iterator > n3p_iterator) {
			break
		}

		// vivek: should this be pair[n5p_iterator][n3p_iterator] or is the current order correct?
		if fc.params.model_details.pair[fc.sequence_encoding2[n3p_iterator]][fc.sequence_encoding2[n5p_iterator]] == 0 {
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
		for pair_table[n5p_iterator] == 0 {
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
	// var rtype [8]int = fc.params.model_details.rtype
	var u1, u2 int
	energy := 0

	u1 = k - i - 1
	u2 = j - l - 1

	ij_pair_type := vrna_get_pair_type_md(fc.sequence_encoding2[i], fc.sequence_encoding2[j], fc.params.model_details)
	kl_pair_type := vrna_get_pair_type_md(fc.sequence_encoding2[l], fc.sequence_encoding2[k], fc.params.model_details)

	energy = E_IntLoop(u1, u2, ij_pair_type, kl_pair_type, fc.sequence_encoding[i+1], fc.sequence_encoding[j-1], fc.sequence_encoding[k-1], fc.sequence_encoding[l+1], fc.params)

	return energy
}

/* vivek:
*	Couldn't find any documentation about what `md.pair` is.
* Created an issue to get more info: https://github.com/ViennaRNA/ViennaRNA/issues/124
* retruns type of base pair (p,q)
 */
func vrna_get_pair_type_md(i, j int, md *vrna_md_t) uint32 {
	var tt uint32 = uint32(md.pair[i][j])

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
	var u, e int
	var type_1 uint32
	var P *vrna_param_t
	var md *vrna_md_t

	P = fc.params
	md = P.model_details
	e = INF

	u = j - i - 1
	type_1 = vrna_get_pair_type_md(fc.sequence_encoding2[i], fc.sequence_encoding2[j], md)

	if !((type_1 == 3) || (type_1 == 4)) {
		e = E_Hairpin(u, type_1, fc.sequence_encoding[i+1], fc.sequence_encoding[j-1], fc.sequence[i-1:], P)
	}

	return e
}

func eval_hp_loop_fake(fc *vrna_fold_compound_t, i, j int) int {
	var e int
	var type_1 uint32
	var P *vrna_param_t
	var md *vrna_md_t

	P = fc.params
	md = P.model_details

	e = INF

	type_1 = vrna_get_pair_type_md(fc.sequence_encoding2[j], fc.sequence_encoding2[i], md)

	/* hairpin-like exterior loop (for cofolding) */
	var si, sj int
	si = fc.sequence_encoding[i+1]
	sj = fc.sequence_encoding[j-1]


	e = vrna_E_ext_stem(type_1, sj, si, P)

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

func cut_in_loop(pair_open_idx int, pair_table []int, sn []uint32) int {
	/* walk around the loop;  return 5' pos of first pair after cut if
	 * cut_point in loop else 0 */
	var n5p_iterator, pair_close_idx int

	n5p_iterator = pair_table[pair_open_idx]
	pair_close_idx = pair_table[pair_open_idx]

	for {
		pair_open_idx = pair_table[n5p_iterator]
		n5p_iterator = pair_open_idx + 1
		for pair_table[n5p_iterator] == 0 {
			n5p_iterator++
		}

		if n5p_iterator == pair_close_idx {
			break
		}
	}
	if sn[pair_open_idx] == sn[n5p_iterator] {
		return 0
	} else {
		return n5p_iterator
	}
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

	for n5p_iterator <= pair_close_idx && pair_table[n5p_iterator] == 0 {
		n5p_iterator++
	}

	/* add bonus energies for first stretch of unpaired nucleotides */
	num_unpaired_nucs += n5p_iterator - pair_open_idx - 1

	for n5p_iterator < pair_close_idx {
		/* p must have a pairing partner */
		n5p_iterator_close_idx = pair_table[n5p_iterator]
		/* get type of base pair (p,q) */
		pair_type := vrna_get_pair_type_md(vc.sequence_encoding2[n5p_iterator], vc.sequence_encoding2[n5p_iterator_close_idx], vc.params.model_details)

		mm5 = vc.sequence_encoding[n5p_iterator-1]
		mm3 = vc.sequence_encoding[n5p_iterator_close_idx+1]

		energy += E_MLstem(pair_type, mm5, mm3, vc.params)

		/* seek to the next stem */
		n5p_iterator = n5p_iterator_close_idx + 1

		for n5p_iterator < pair_close_idx && pair_table[n5p_iterator] == 0 {
			n5p_iterator++
		}
		num_unpaired_nucs += n5p_iterator - n5p_iterator_close_idx - 1 /* add unpaired nucleotides */
	}

	if pair_open_idx > 0 {
		/* actual closing pair */
		pair_type := vrna_get_pair_type_md(vc.sequence_encoding2[pair_close_idx], vc.sequence_encoding2[pair_open_idx], vc.params.model_details)
		mm5 = vc.sequence_encoding[pair_close_idx-1]
		mm3 = vc.sequence_encoding[pair_open_idx+1]

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
