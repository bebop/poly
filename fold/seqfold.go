package fold

import (
	"fmt"
	"math"
	"strings"

	"github.com/TimothyStiles/poly/checks"
)

/*
MultibranchEnergies holds the a, b, c, d in a linear multi-branch energy

change function.

A - number of helices in the loop
B - number of unpaired nucleotides in the loop
C - coxial stacking in the loop
D - terminal mismatch contributions
E - base composition of the unpaired nucleotides (probably neglible?)
Inferred from:
Supplemental Material: Annu.Rev.Biophs.Biomol.Struct.33:415-40
doi: 10.1146/annurev.biophys.32.110601.141800
The Thermodynamics of DNA Structural Motifs, SantaLucia and Hicks, 2004
*/
type MultibranchEnergies struct {
	A, B, C, D float64
}

// Energy holds two energies, enthaply and entropy
// SantaLucia & Hicks (2004), Annu. Rev. Biophys. Biomol. Struct 33: 415-440
type Energy struct {
	// enthalpy
	EnthalpyH float64
	// entropy
	EntropyS float64
}

// MatchingBasepairEnergy is the energy of matching base pairs
type MatchingBasepairEnergy map[string]Energy

// LoopEnergy is a map[int]Energy where the int is the length of the loop
type LoopEnergy map[int]Energy

// Energies holds the needed energy maps, BpEnergy and LoopEnergy, to compute
// the folding, it also holds the complement map for the kind of sequence, RNA
// or DNA
type Energies struct {
	BulgeLoops         LoopEnergy
	Complement         map[byte]byte
	DanglingEnds       MatchingBasepairEnergy
	HairpinLoops       LoopEnergy
	Multibranch        MultibranchEnergies
	InternalLoops      LoopEnergy
	InternalMismatches MatchingBasepairEnergy
	NearestNeighbors   MatchingBasepairEnergy
	TerminalMismatches MatchingBasepairEnergy
	TriTetraLoops      MatchingBasepairEnergy
}

// Subsequence represent an interval of bases in the sequence that can contain
// a inward structure.
type Subsequence struct {
	Start, End int
}

// A single structure with a free energy, description, and inward children
type NucleicAcidStructure struct {
	// Description is the description of the NucleicAcidStructure
	Description string
	// Inner is list of inner structure represented as intervals in the
	// sequence
	Inner []Subsequence
	// Energy is the energy of the NucleicAcidStructure
	Energy float64
}

func (structure NucleicAcidStructure) Equal(other NucleicAcidStructure) bool {
	if len(structure.Inner) != len(other.Inner) {
		return false
	}
	for i, val := range structure.Inner {
		if val != other.Inner[i] {
			return false
		}
	}
	return structure.Energy == other.Energy
}

func (structure NucleicAcidStructure) Valid() bool {
	return structure.Energy != math.Inf(1) && structure.Energy != math.Inf(-1)
}

func (structure NucleicAcidStructure) String() string {
	i, j := "", ""
	if len(structure.Inner) > 0 {
		i, j = fmt.Sprint(structure.Inner[0].Start), fmt.Sprint(structure.Inner[0].End)
	}
	return fmt.Sprintf("%4s %4s % 6.2f  %-15s", i, j, structure.Energy, structure.Description)
}

// MultiString returns all the fields as strings, this is useful to output the
// result in a tabular fashion using the same format string.
func (structure NucleicAcidStructure) MultiString() (string, string, string, string) {
	i, j := "", ""
	if len(structure.Inner) > 0 {
		i, j = fmt.Sprint(structure.Inner[0].Start), fmt.Sprint(structure.Inner[0].End)
	}
	return i, j, fmt.Sprintf("%6.2f", structure.Energy), structure.Description
}

// DefaultStructure is the default (zero value) nucleic acid structure, it used
// mostly to initialize the caches, see FoldingContext
var DefaultStructure = NucleicAcidStructure{
	Description: "",
	Energy:      math.Inf(-1),
}

// InvalidStructure represent an invalid nucleic acid structure
var InvalidStructure = NucleicAcidStructure{
	Description: "",
	Energy:      math.Inf(1),
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
	var energyMap Energies
	switch {
	case checks.IsDNA(seq):
		energyMap = DNAEnergies
	case checks.IsRNA(seq):
		energyMap = RNAEnergies
	default:
		return FoldContext{}, fmt.Errorf("the sequence %s is not RNA or DNA", seq)
	}

	var (
		sequenceLength = len(seq)
		vCache         = make([][]NucleicAcidStructure, sequenceLength)
		wCache         = make([][]NucleicAcidStructure, sequenceLength)
		row            = make([]NucleicAcidStructure, sequenceLength)
	)
	for nucleicAcidIndex := 0; nucleicAcidIndex < sequenceLength; nucleicAcidIndex++ {
		row[nucleicAcidIndex] = DefaultStructure
	}
	for j := 0; j < sequenceLength; j++ {
		vCache[j] = make([]NucleicAcidStructure, sequenceLength)
		copy(vCache[j], row)

		wCache[j] = make([]NucleicAcidStructure, sequenceLength)
		copy(wCache[j], row)
	}
	ret := FoldContext{
		Energies: energyMap,
		Seq:      seq,
		V:        vCache,
		W:        wCache,
		T:        temp + 273.15, // kelvin
	}

	// fill the cache
	_, err := W(0, sequenceLength-1, ret)
	if err != nil {
		return FoldContext{}, fmt.Errorf("error filling the caches for the FoldingContext: %w", err)
	}
	return ret, nil

}
