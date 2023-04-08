package fold

import (
	"fmt"
	"math"
	"strings"

	"github.com/TimothyStiles/poly/checks"
)

// MultibranchEnergies holds the a, b, c, d in a linear multi-branch energy
// change function.
// Inferred from:
// Supplemental Material: Annu.Rev.Biophs.Biomol.Struct.33:415-40
// doi: 10.1146/annurev.biophys.32.110601.141800
// The Thermodynamics of DNA Structural Motifs, SantaLucia and Hicks, 2004
type MultibranchEnergies struct {
	A, B, C, D float64
}

// Energy holds two energies, enthaply and entropy
// SantaLucia & Hicks (2004), Annu. Rev. Biophys. Biomol. Struct 33: 415-440
type Energy struct {
	// enthaply
	H float64
	// entropy
	S float64
}

// BPEnergy is the energy of matching base pairs
type BPEnergy map[string]Energy

// LoopEnergy is a map[int]Energy where the int is the length of the loop
type LoopEnergy map[int]Energy

// Energies holds the needed energy maps, BpEnergy and LoopEnergy, to compute
// the folding, it also holds the complement map for the kind of sequence, RNA
// or DNA
type Energies struct {
	BulgeLoops    LoopEnergy
	Complement    map[byte]byte
	DE            BPEnergy
	HairpinLoops  LoopEnergy
	Multibranch   MultibranchEnergies
	InternalLoops LoopEnergy
	INTERNAL_MM   BPEnergy
	NN            BPEnergy
	TERMINAL_MM   BPEnergy
	TriTetraLoops BPEnergy
}

// Subsequence represent an interval of bases in the sequence that can contain
// a inward structure.
type Subsequence struct {
	Start, End int
}

// A single structure with a free energy, description, and inward children
type NucleicAcidStructure struct {
	// Desc is the description of the NucleicAcidStructure
	Desc string
	// Inner is list of inner structure represented as intervals in the
	// sequence
	Inner []Subsequence
	// E is the energy of the NucleicAcidStructure
	E float64
}

func (s NucleicAcidStructure) Equal(other NucleicAcidStructure) bool {
	if len(s.Inner) != len(other.Inner) {
		return false
	}
	for i, val := range s.Inner {
		if val != other.Inner[i] {
			return false
		}
	}
	return s.E == other.E
}

func (s NucleicAcidStructure) Valid() bool {
	return s.E != math.Inf(1) && s.E != math.Inf(-1)
}

func (s NucleicAcidStructure) String() string {
	i, j := "", ""
	if len(s.Inner) > 0 {
		i, j = fmt.Sprint(s.Inner[0].Start), fmt.Sprint(s.Inner[0].End)
	}
	return fmt.Sprintf("%4s %4s % 6.2f  %-15s", i, j, s.E, s.Desc)
}

// MultiString returns all the fields as strings, this is useful to output the
// result in a tabular fashion using the same format string.
func (s NucleicAcidStructure) MultiString() (string, string, string, string) {
	i, j := "", ""
	if len(s.Inner) > 0 {
		i, j = fmt.Sprint(s.Inner[0].Start), fmt.Sprint(s.Inner[0].End)
	}
	return i, j, fmt.Sprintf("%6.2f", s.E), s.Desc
}

// STRUCT_DEFAULT is the default (zero value) nucleic acid structure, it used
// mostly to initialize the caches, see FoldingContext
var STRUCT_DEFAULT = NucleicAcidStructure{
	Desc: "",
	E:    math.Inf(-1),
}

// STRUCT_NULL represent an invalid nucleic acid structure
var STRUCT_NULL = NucleicAcidStructure{
	Desc: "",
	E:    math.Inf(1),
}

// FoldContext holds the energy caches, energy maps, sequence, and temperature
// needed in order to compute the folding energy and structures.
type FoldContext struct {
	Energies Energies
	Seq      string
	V        [][]NucleicAcidStructure
	W        [][]NucleicAcidStructure
	T        float64
}

// NewFoldingContext returns a FoldContext ready to use, in case of error
// the returned FoldingContext is empty.
func NewFoldingContext(seq string, temp float64) (FoldContext, error) {
	seq = strings.ToUpper(seq)

	// figure out whether it's DNA or RNA, choose energy map
	var emap Energies
	switch {
	case checks.IsDNA(seq):
		emap = DNAEnergies
	case checks.IsRNA(seq):
		emap = RNAEnergies
	default:
		return FoldContext{}, fmt.Errorf("the sequence %s is not RNA or DNA", seq)
	}

	var (
		n       = len(seq)
		v_cache = make([][]NucleicAcidStructure, n)
		w_cache = make([][]NucleicAcidStructure, n)
		row     = make([]NucleicAcidStructure, n)
	)
	for i := 0; i < n; i++ {
		row[i] = STRUCT_DEFAULT
	}
	for j := 0; j < n; j++ {
		v_cache[j] = make([]NucleicAcidStructure, n)
		copy(v_cache[j], row)

		w_cache[j] = make([]NucleicAcidStructure, n)
		copy(w_cache[j], row)
	}
	ret := FoldContext{
		Energies: emap,
		Seq:      seq,
		V:        v_cache,
		W:        w_cache,
		T:        temp + 273.15, // kelvin
	}

	// fill the cache
	_, err := W(0, n-1, ret)
	if err != nil {
		return FoldContext{}, fmt.Errorf("error filling the caches for the FoldingContext: %w", err)
	}
	return ret, nil

}
