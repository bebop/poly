package secondary_structure

/******************************************************************************
July, 01, 2021

This file defines the structs needed to contain information about a RNA's
secondary structure.

The main RNA secondary structures are `MultiLoop`s, `Hairpin`s, and
`SingleStrandedRegion`s. Hairpins and Multi-loops both can optionally have a
`Stem`.

A `Stem` consists of a list of `StemStructure`s. A `StemStructure` is delimited
by a closing and enclosed base pair with the requirement that there are no base
pairs between the closing and enclosed base pair. A `StemStrcture` is classified
into a `StemStructureType` based on the number of unpaired nucleotides between
the closing and enclosed base pairs.

******************************************************************************/

// SecondaryStructure is composed of a list of `MultiLoop`s, `Hairpin`s,
// and `SingleStrandedRegion`s. Note that since Go doesn't support inheritance,
// we use `interface{}` as the type for the structures list, but the only types
// that are allowed/used are `MultiLoop`, `Hairpin` and `SingleStrandedRegion`.
type SecondaryStructure struct {
	Structures          []interface{}
	Length              int
	ExteriorLoopsEnergy int // the free energy (in dcal / mol) of the exterior loops present in the secondary structure
	Energy              int // the free energy (in dcal / mol) of the entire secondary structure
}

// MultiLoop contains all the information needed to denote a multi-loop in a
// RNA's secondary structure. It consists of a `Stem` and a list of
// substructures. Note that a `Multiloop` will always contain atleast one
// substructure.
type MultiLoop struct {
	Stem                                                  Stem
	SubstructuresFivePrimeIdx, SubstructuresThreePrimeIdx int
	Substructures                                         SecondaryStructure
	Energy                                                int // free energy (in dcal / mol) only from the multi-loop (doesn't include free energy from the substrctures or stem)
}

// Hairpin contains all the information needed to denote a hairpin loop in a
// RNA's secondary structure. It consists of a `Stem` and a single stranded
// region that forms the loop of the structure.
// Sometimes a `Hairpin` may only consist of a Stem without a single stranded
// region. In such cases, the `SingleStrandedFivePrimeIdx` and
// `SingleStrandedThreePrimeIdx` of the hairpin are set to `-1`.
type Hairpin struct {
	Stem                                                    Stem
	SingleStrandedFivePrimeIdx, SingleStrandedThreePrimeIdx int
	Energy                                                  int // free energy in dcal / mol only from the hairpin loop (doesn't include free energy from the stem)
}

// SingleStrandedRegion contains all the information needed to denote a
// single stranded region in a RNA's secondary structure.
type SingleStrandedRegion struct {
	FivePrimeIdx, ThreePrimeIdx int
}

// Stem contains all the information needed to denote the stems of a `Hairpin`
// or `Multiloop`. It is not a "top-level" structure of a RNA and only exists
// as part of a `Hairpin` or `Multiloop`. The closing pairs denote where the
// stem starts and enclosed pairs where the stem ends. The actual stem
// consists of a list of `StemStructure`s.
//
// Note that a `Stem` may not contain any stem structures. This occurs in cases
// where there is only one base pair that delimits a `Hairpin` or `MultiLoop`.
// (For example, "eee(hhh)eee" would be a Hairpin (Stem at indexs 3 and 7) and
// "ee(mmm((hh))mm)ee" would be a Multiloop (Stem at indexs 2 and 14) whose
// stems don't contain any structures)
// In such cases, the stem will only have its `ClosingFivePrimeIdx` and
// `ClosingThreePrimeIdx` set. The `EnclosedFivePrimeIdx` and
// `EnclosedThreePrimeIdx` will be set to -1, and the list of `StemStructres`
// will be empty.
type Stem struct {
	ClosingFivePrimeIdx, EnclosedFivePrimeIdx   int
	EnclosedThreePrimeIdx, ClosingThreePrimeIdx int
	Structures                                  []StemStructure
	Energy                                      int // free energy (in dcal / mol) all the substructures present in the stem
}

// StemStructure contains all the information needed to denote the substructures
// present in a `Stem`. A `StemStructure` is delimited by a closing and
// enclosed base pair with the requirement that there are no base pairs between
// the closing and enclosed base pair.
//
// A `StemStrcture` is classified into a `StemStructureType` based on the
// number of unpaired nucleotides between the closing and enclosed base pairs.
// (See (*StemStructure).setStructureType for more details)
type StemStructure struct {
	ClosingFivePrimeIdx, EnclosedFivePrimeIdx   int
	EnclosedThreePrimeIdx, ClosingThreePrimeIdx int
	Type                                        StemStructureType
	Energy                                      int // free energy in dcal / mol
}

// StemStructureType denotes the type of a `StemStructure`.
type StemStructureType int

const (
	// StackingPair is the type of a `StemStructure` where there are no unpaired
	// nucleotides between the closing and enclosed base pairs of the
	// `StemStructure`.
	StackingPair StemStructureType = iota
	// Bulge is the type of a `StemStructure` where there is more than one
	// unpaired nucleotide on one 'side' of the `StemStructure` and no unpaired
	// nucleotides on the other 'side'.
	Bulge
	// Interior1x1Loop is the type of a `StemStructure` where there is one
	// unpaired nucleotide on both 'sides' of the `StemStructure`.
	Interior1x1Loop
	// Interior2x1Loop is the type of a `StemStructure` where there are two
	// unpaired nucleotides on one 'side' and one unpaired nucleotides on the
	// other 'side' of the `StemStructure`.
	Interior2x1Loop
	// Interior1xnLoop is the type of a `StemStructure` where there is one
	// unpaired nucleotides on one 'side' and more than two unpaired nucleotides
	// on the other 'side' of the `StemStructure`.
	Interior1xnLoop
	// Interior2x2Loop is the type of a `StemStructure` where there are two
	// unpaired nucleotides on both 'sides' of the `StemStructure`.
	Interior2x2Loop
	// Interior2x3Loop is the type of a `StemStructure` where there are two
	// unpaired nucleotides on one 'side' and three unpaired nucleotides on the
	// other 'side' of the `StemStructure`.
	Interior2x3Loop
	// GenericInteriorLoop is the type of a `StemStructure` which is not denoted
	// by `StackingPair`, `Bulge`, `Interior1x1Loop`, `Interior2x1Loop`,
	// `Interior1xnLoop`, `Interior2x2Loop`, or `Interior2x3Loop`.
	// Thus, the `StemStructure` can be one of:
	// * two unpaired nucleotides on one 'side' and more than three on the other
	//   (2x4, 2x5, ..., 2xn interior loops)
	// * three unpaired nucleotides on one 'side' and three or more on the other
	//   (3x3, 3x4, ..., 3xn interior loops)
	GenericInteriorLoop
)

// setStructureType sets the `Type` field of a `StemStructure`.
func (structure *StemStructure) setStructureType() {
	nbUnpairedFivePrime := structure.EnclosedFivePrimeIdx - structure.ClosingFivePrimeIdx - 1
	nbUnpairedThreePrime := structure.ClosingThreePrimeIdx - structure.EnclosedThreePrimeIdx - 1

	var nbUnpairedLarger, nbUnpairedSmaller int

	if nbUnpairedFivePrime > nbUnpairedThreePrime {
		nbUnpairedLarger = nbUnpairedFivePrime
		nbUnpairedSmaller = nbUnpairedThreePrime
	} else {
		nbUnpairedLarger = nbUnpairedThreePrime
		nbUnpairedSmaller = nbUnpairedFivePrime
	}

	switch nbUnpairedSmaller {
	case 0:
		if nbUnpairedLarger == 0 {
			structure.Type = StackingPair
		} else {
			structure.Type = Bulge
		}
	case 1:
		switch nbUnpairedLarger {
		case 1:
			structure.Type = Interior1x1Loop
		case 2:
			structure.Type = Interior2x1Loop
		default:
			structure.Type = Interior1xnLoop
		}
	case 2:
		switch nbUnpairedLarger {
		case 2:
			structure.Type = Interior2x2Loop
		case 3:
			structure.Type = Interior2x3Loop
		default:
			structure.Type = GenericInteriorLoop
		}

	default:
		structure.Type = GenericInteriorLoop
	}
}

// NewStemStructure is a wrapper to create a `StemStructure` and call the
// functions (`(*StemStructure).setStructureType`) required to initalise the
// struct
func NewStemStructure(closingFivePrimeIdx, closingThreePrimeIdx,
	enclosedFivePrimeIdx, enclosedThreePrimeIdx int) StemStructure {

	stemStructure := StemStructure{
		ClosingFivePrimeIdx:   closingFivePrimeIdx,
		EnclosedFivePrimeIdx:  enclosedFivePrimeIdx,
		EnclosedThreePrimeIdx: enclosedThreePrimeIdx,
		ClosingThreePrimeIdx:  closingThreePrimeIdx,
	}

	stemStructure.setStructureType()

	return stemStructure
}
