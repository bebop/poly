package fold

import (
	"fmt"
	"math"
	"strings"

	"github.com/TimothyStiles/poly/checks"
)

/*
multibranchEnergies holds the a, b, c, d in a linear multi-branch energy

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
type multibranchEnergies struct {
	helicesCount, unpairedCount, coaxialStackCount, terminalMismatchCount float64
}

// Energy holds two energies, enthaply and entropy
// SantaLucia & Hicks (2004), Annu. Rev. Biophys. Biomol. Struct 33: 415-440
type Energy struct {
	// enthalpy
	EnthalpyH float64
	// entropy
	EntropyS float64
}

// matchingBasepairEnergy is the energy of matching base pairs
type matchingBasepairEnergy map[string]Energy

// loopEnergy is a map[int]Energy where the int is the length of the loop
type loopEnergy map[int]Energy

// complementFunc is function to translate a base in its complement
type complementFunc func(rune) rune

// energies holds the needed energy maps, BpEnergy and loopEnergy, to compute
// the folding, it also holds the complement map for the kind of sequence, rna
// or DNA
type energies struct {
	bulgeLoops         loopEnergy
	complement         complementFunc
	danglingEnds       matchingBasepairEnergy
	hairpinLoops       loopEnergy
	multibranch        multibranchEnergies
	internalLoops      loopEnergy
	internalMismatches matchingBasepairEnergy
	nearestNeighbors   matchingBasepairEnergy
	terminalMismatches matchingBasepairEnergy
	triTetraLoops      matchingBasepairEnergy
}

// Subsequence represent an interval of bases in the sequence that can contain
// a inward structure.
type Subsequence struct {
	Start, End int
}

// NucleicAcidStructure is single structure with a free energy, description, and inward children
type NucleicAcidStructure struct {
	// Description is the description of the NucleicAcidStructure
	Description string
	// Inner is list of inner structure represented as intervals in the
	// sequence
	Inner []Subsequence
	// Energy is the energy of the NucleicAcidStructure
	Energy float64
}

// Equal returns true if two NucleicAcidStructures are equal
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

// Valid returns true if the NucleicAcidStructure is valid
func (structure NucleicAcidStructure) Valid() bool {
	return structure.Energy != math.Inf(1) && structure.Energy != math.Inf(-1)
}

// defaultStructure is the default (zero value) nucleic acid structure, it used
// mostly to initialize the caches, see FoldingContext
var defaultStructure = NucleicAcidStructure{
	Description: "",
	Energy:      math.Inf(-1),
}

// invalidStructure represent an invalid nucleic acid structure
var invalidStructure = NucleicAcidStructure{
	Description: "",
	Energy:      math.Inf(1),
}

// context holds the energy caches, energy maps, sequence, and temperature
// needed in order to compute the folding energy and structures.
type context struct {
	energies                   energies
	Seq                        string
	pairedMinimumFreeEnergyV   [][]NucleicAcidStructure
	unpairedMinimumFreeEnergyW [][]NucleicAcidStructure
	T                          float64
}

// newFoldingContext returns a context ready to use, in case of error
// the returned FoldingContext is empty.
func newFoldingContext(seq string, temp float64) (context, error) {
	seq = strings.ToUpper(seq)

	// figure out whether it's DNA or rna, choose energy map
	var energyMap energies
	switch {
	case checks.IsDNA(seq):
		energyMap = dnaEnergies
	case checks.IsRNA(seq):
		energyMap = rnaEnergies
	default:
		return context{}, fmt.Errorf("the sequence %s is not RNA or DNA", seq)
	}

	var (
		sequenceLength = len(seq)
		vCache         = make([][]NucleicAcidStructure, sequenceLength)
		wCache         = make([][]NucleicAcidStructure, sequenceLength)
		row            = make([]NucleicAcidStructure, sequenceLength)
	)
	for nucleicAcidIndex := 0; nucleicAcidIndex < sequenceLength; nucleicAcidIndex++ {
		row[nucleicAcidIndex] = defaultStructure
	}
	for j := 0; j < sequenceLength; j++ {
		vCache[j] = make([]NucleicAcidStructure, sequenceLength)
		copy(vCache[j], row)

		wCache[j] = make([]NucleicAcidStructure, sequenceLength)
		copy(wCache[j], row)
	}
	ret := context{
		energies:                   energyMap,
		Seq:                        seq,
		pairedMinimumFreeEnergyV:   vCache,
		unpairedMinimumFreeEnergyW: wCache,
		T:                          temp + 273.15, // kelvin
	}

	// fill the cache
	_, err := unpairedMinimumFreeEnergyW(0, sequenceLength-1, ret)
	if err != nil {
		return context{}, fmt.Errorf("error filling the caches for the FoldingContext: %w", err)
	}
	return ret, nil

}
