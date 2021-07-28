/*
Package secondary_structure provides the structs needed to contain
information about a RNA's secondary structure

Overview of the structs

The struct that contains information of a RNA's secondary structure is
`SecondaryStructure`. The field `Structures` contains a list of the main
RNA secondary structures (`*MultiLoop`, `*Hairpin`, and `*SingleStrandedRegion`).
`Hairpin`s and `MultiLoop`s both can optionally have a `Stem`.

A `Stem` consists of a list of `StemStructure`s. A `StemStructure` consists
of a closing and enclosed base pair with the requirement that there are
no base pairs between the closing and enclosed base pair.

See the declaration of the structs for detailed information on their
definition.


Explanation of the energy fields of the structs

The energy fields of the structs are only used in the `SecondaryStructure`
returned from the func `MinimumFreeEnergy` in the subpackage `mfe` in
`poly`. The func `MinimumFreeEnergy` reads energy paramater values from files
specified in the `RNAfold parameter file v2.0` file format in the
`energy_params` subpackage in `poly`. This file format specifies energy
values in the unit deca-cal / mol (with the `int` type).

Due to Golang's inaccuracies with adding and subtracting `float64` values,
a design decision was made to specify the energy fields of the structs
defined in this file in the unit deca-cal / mol (with the `int` type).

Thus, to convert to the 'standard' energy unit of kcal/mol, the energy
values have to be converted to type `float32` or `float64` and divided by `100`.


*/
package secondary_structure

// SecondaryStructure is composed of a list of `MultiLoop`s, `Hairpin`s,
// and `SingleStrandedRegion`s. Note that since Go doesn't support inheritance,
// we use `interface{}` as the type for the structures list, but the only types
// that are allowed/used are `*MultiLoop`, `*Hairpin` and `*SingleStrandedRegion`.
//
// The free energy of the entire secondary structure (including the energy
// of the exterior loops of the secondary structure) is stored in the `Energy`
// field of the struct. Since the energy of the exterior loops of the
// secondary structure is not part of any of the `SecondaryStructure`'s
// `Structures` (and belongs to the entire SecondaryStructure struct), it is
// stored in the field `ExteriorLoopsEnergy`.
type SecondaryStructure struct {
	Structures          []interface{}
	Length              int
	ExteriorLoopsEnergy int // the free energy (in dcal / mol) of the exterior loops present in the secondary structure
	Energy              int // the free energy (in dcal / mol) of the entire secondary structure
}

// MultiLoop contains all the information needed to denote a multi-loop in a
// RNA's secondary structure. It consists of a `Stem` and a list of
// substructures. A `Multiloop` will always contain at least one substructure.
// The substructures that can be present in a Multiloop are `*Hairpin`s,
// `*SingleStrandedRegion`s and `*MultiLoop`s.
//
// The total free energy of a Multiloop can be calculated by summing the energy
// of the Multiloop's stem (present in the `Energy` field of the `Stem`),
// the energy of the substructures (the `Energy` field of `Substructures`),
// and the energy of the Multiloop (the `Energy` field of this struct).
type MultiLoop struct {
	Stem                Stem
	Substructures       []interface{}
	Energy              int // free energy (in dcal / mol) only from the multi-loop (doesn't include free energy from the substrctures or stem)
	SubstructuresEnergy int // free energy (in dcal / mol) only from the substructures of the multi-loop (doesn't include free energy from the stem or multi-loop)
}

// Hairpin contains all the information needed to denote a hairpin loop in a
// RNA's secondary structure. It consists of a `Stem` and a single stranded
// region that forms the loop of the structure.
//
// Sometimes a `Hairpin` may only consist of a Stem without a single stranded
// region. In such cases, the `SingleStrandedFivePrimeIdx` and
// `SingleStrandedThreePrimeIdx` of the hairpin are set to `-1`.
//
// The total free energy of a Hairpin can be calculated by summing the energy
// of the Hairpin's stem (present in the `Energy` field of the `Stem`),
// and the energy of the single stranded region of a Hairpin's loop (the
// `Energy` field of this struct).
type Hairpin struct {
	Stem                                                    Stem
	SingleStrandedFivePrimeIdx, SingleStrandedThreePrimeIdx int
	Energy                                                  int // free energy in dcal / mol only from the hairpin loop (doesn't include free energy from the stem)
}

// SingleStrandedRegion contains all the information needed to denote a
// single stranded region in a RNA's secondary structure.
// At the minimum, a `SingleStrandedRegion` consists of of a single
// unpaired nucleotide. In such a case, `FivePrimeIdx` == `ThreePrimeIdx`.
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
// For example,
// dot-bracket structure:
// . . . ( . . . ) . .
// annotated structure:
// e e e ( h h h ) e e
// index:
// 0 1 2 3 4 5 6 7 8 9
// would be a Hairpin (with a Stem with the closing base pair at indexs 3 and
// 7) whose stem don't contain any structures.
//
// For example,
// dot-bracket structure:
// . . ( . . . ( ( . . )  )  .  .  )  .  .
// annotated structure:
// . . ( m m m ( ( h h )  )  m  m  )  e  e
// index:
// 0 1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16
// would be a Multiloop (with a Stem with the closing base pair at indexs 2 and
// 14) whose stem don't contain any structures.
// In such cases, the stem will only have its `ClosingFivePrimeIdx` and
// `ClosingThreePrimeIdx` set. The `EnclosedFivePrimeIdx` and
// `EnclosedThreePrimeIdx` will be set to -1, and the list of `StemStructures`
// will be empty.
type Stem struct {
	ClosingFivePrimeIdx, EnclosedFivePrimeIdx   int
	EnclosedThreePrimeIdx, ClosingThreePrimeIdx int
	Structures                                  []StemStructure
	Energy                                      int // free energy (in dcal / mol) all the substructures present in the stem
}

// StemStructure contains all the information needed to denote the substructures
// present in a `Stem`. A `StemStructure` always consists of a closing and
// enclosed base pair with the requirement that there are no base pairs between
// the closing and enclosed base pair.
//
// A `StemStructure` is classified into a `StemStructureType` based on the
// number of unpaired nucleotides between the closing and enclosed base pairs.
// (See (*StemStructure).setStructureType for more details)
type StemStructure struct {
	ClosingFivePrimeIdx, EnclosedFivePrimeIdx   int
	EnclosedThreePrimeIdx, ClosingThreePrimeIdx int
	NBUnpairedFivePrime, NBUnpairedThreePrime   int // the number of unpaired nucleotides on the five and three prime ends
	Type                                        StemStructureType
	Energy                                      int // free energy of the `StemStructure` in dcal / mol
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

// setStructureType sets the `Type` field of a `StemStructure` based on the
// number of unpaired nucleotides between the closing and enclosed base pairs.
func (structure *StemStructure) setStructureType() {
	nbUnpairedFivePrime := structure.EnclosedFivePrimeIdx - structure.ClosingFivePrimeIdx - 1
	structure.NBUnpairedFivePrime = nbUnpairedFivePrime
	nbUnpairedThreePrime := structure.ClosingThreePrimeIdx - structure.EnclosedThreePrimeIdx - 1
	structure.NBUnpairedThreePrime = nbUnpairedThreePrime

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
// functions (`(*StemStructure).setStructureType`) required to initialize the
// struct.
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
