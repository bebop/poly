package seconday_structure

// A RNA's SecondaryStructure is composed of a list of MultiLoops, HairpinLoops,
// and SingleStrandedRegions. Note that since Go doesn't support inheritance,
// we use `interface{}` as the type for the list, but the only types that are
// allowed/used are `MultiLoop`, `HairpinLoop` and `SingleStrandedRegion`
type SecondaryStructure struct {
	Structures []interface{}
	Length     int
}

type MultiLoop struct {
	StemFivePrimeIdx, StemThreePrimeIdx                   int
	Stem                                                  Stem
	SubstructuresFivePrimeIdx, SubstructuresThreePrimeIdx int
	Substructures                                         SecondaryStructure
}

type HairpinLoop struct {
	StemFivePrimeIdx, StemThreePrimeIdx                     int
	Stem                                                    Stem
	SingleStrandedFivePrimeIdx, SingleStrandedThreePrimeIdx int
}

type SingleStrandedRegion struct {
	FivePrimeIdx, ThreePrimeIdx int
}

type Stem struct {
	ClosingFivePrimeIdx, EnclosedFivePrimeIdx   int
	EnclosedThreePrimeIdx, ClosingThreePrimeIdx int
	structures                                  []StemStructure
}

type StemStructureType int

const (
	StackingPair StemStructureType = iota
	Bulge
	Interior1x1Loop
	Interior2x1Loop
	Interior1xnLoop
	Interior2x2Loop
	Interior2x3Loop
	GenericInteriorLoop
)

type StemStructure struct {
	ClosingFivePrimeIdx, EnclosedFivePrimeIdx   int
	EnclosedThreePrimeIdx, ClosingThreePrimeIdx int
	StructureType                               StemStructureType
}

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

	if nbUnpairedLarger == 0 {
		// stacking pair
		structure.StructureType = StackingPair
		return
	}

	if nbUnpairedSmaller == 0 {
		// bulge
		structure.StructureType = Bulge
		return
	} else {
		// interior loop
		if nbUnpairedSmaller == 1 {
			if nbUnpairedLarger == 1 {
				// 1x1 loop
				structure.StructureType = Interior1x1Loop
				return
			}

			if nbUnpairedLarger == 2 {
				// 2x1 loop
				structure.StructureType = Interior2x1Loop
				return
			} else {
				// 1xn loop
				structure.StructureType = Interior1xnLoop
				return
			}
		} else if nbUnpairedSmaller == 2 {
			if nbUnpairedLarger == 2 {
				// 2x2 loop
				structure.StructureType = Interior2x2Loop
				return
			} else if nbUnpairedLarger == 3 {
				// 2x3 loop
				structure.StructureType = Interior2x3Loop
				return
			}
		}

		{
			/* generic interior loop (no else here!)*/
			structure.StructureType = GenericInteriorLoop
			return
		}
	}
}
