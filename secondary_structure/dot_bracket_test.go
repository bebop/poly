package secondary_structure

import (
	"fmt"
)

func ExampleSecondaryStructureFromDotBracket() {
	dotBracket := "..((((...))))...((........)).."

	// dotBracket with index
	// . . ( ( ( ( . . . )  )  )  )  .  .  .  (  (  .  .  .  .  .  .  .  .  )  )  .  .
	// 0 1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22 23 24 25 26 27 28 29

	actualSecondaryStructure := SecondaryStructure{
		Length:              30,
		ExteriorLoopsEnergy: 0,
		Energy:              0,
		Structures: []interface{}{
			&SingleStrandedRegion{
				FivePrimeIdx:  0,
				ThreePrimeIdx: 1,
			},
			&Hairpin{
				Stem: Stem{
					ClosingFivePrimeIdx:   2,
					EnclosedFivePrimeIdx:  5,
					EnclosedThreePrimeIdx: 9,
					ClosingThreePrimeIdx:  12,
					Structures: []StemStructure{
						StemStructure{
							ClosingFivePrimeIdx:   2,
							EnclosedFivePrimeIdx:  3,
							EnclosedThreePrimeIdx: 11,
							ClosingThreePrimeIdx:  12,
							Type:                  StackingPair, Energy: 0,
						},
						StemStructure{
							ClosingFivePrimeIdx:   3,
							EnclosedFivePrimeIdx:  4,
							EnclosedThreePrimeIdx: 10,
							ClosingThreePrimeIdx:  11,
							Type:                  StackingPair, Energy: 0,
						},
						StemStructure{
							ClosingFivePrimeIdx:   4,
							EnclosedFivePrimeIdx:  5,
							EnclosedThreePrimeIdx: 9,
							ClosingThreePrimeIdx:  10,
							Type:                  StackingPair, Energy: 0,
						},
					},
					Energy: 0,
				},
				SingleStrandedFivePrimeIdx:  6,
				SingleStrandedThreePrimeIdx: 8,
				Energy:                      0,
			},
			&SingleStrandedRegion{
				FivePrimeIdx:  13,
				ThreePrimeIdx: 15,
			},
			&Hairpin{
				Stem: Stem{
					ClosingFivePrimeIdx:   16,
					EnclosedFivePrimeIdx:  17,
					EnclosedThreePrimeIdx: 26,
					ClosingThreePrimeIdx:  27,
					Structures: []StemStructure{
						StemStructure{
							ClosingFivePrimeIdx:   16,
							EnclosedFivePrimeIdx:  17,
							EnclosedThreePrimeIdx: 26,
							ClosingThreePrimeIdx:  27,
							Type:                  StackingPair, Energy: 0,
						},
					},
					Energy: 0,
				},
				SingleStrandedFivePrimeIdx:  18,
				SingleStrandedThreePrimeIdx: 25,
				Energy:                      0,
			},
			&SingleStrandedRegion{
				FivePrimeIdx:  28,
				ThreePrimeIdx: 29,
			},
		},
	}

	annotatedStructure, secondaryStructure, err := SecondaryStructureFromDotBracket(dotBracket)
	if err != nil {
		panic(err)
	}

	fmt.Println(annotatedStructure)
	// check if predicted and actual structures are same
	fmt.Println(isSecondaryStructuresEqual(actualSecondaryStructure, *secondaryStructure))

	// Output:
	// ee((((hhh))))eee((hhhhhhhh))ee
	// true
}

func ExampleSecondaryStructureFromDotBracket_multiLoop() {
	dotBracket := "..(..((((...))))...((........))..)."

	// dotBracket with index
	// . . ( . . ( ( ( ( .  .  .  )  )  )  )  .  .  .  (  (  .  .  .  .  .  .  .  .  )  )  .  .  )  .
	// 0 1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22 23 24 25 26 27 28 29 30 31 32 33 34

	actualSecondaryStructure := SecondaryStructure{
		Length:              35,
		ExteriorLoopsEnergy: 0,
		Energy:              0,
		Structures: []interface{}{
			&SingleStrandedRegion{
				FivePrimeIdx:  0,
				ThreePrimeIdx: 1,
			},
			&MultiLoop{
				Energy: 0,
				Stem: Stem{
					ClosingFivePrimeIdx:   2,
					EnclosedFivePrimeIdx:  -1,
					EnclosedThreePrimeIdx: -1,
					ClosingThreePrimeIdx:  33,
					Structures:            []StemStructure{},
					Energy:                0,
				},
				Substructures: []interface{}{
					&SingleStrandedRegion{
						FivePrimeIdx:  3,
						ThreePrimeIdx: 4,
					},
					&Hairpin{
						Stem: Stem{
							ClosingFivePrimeIdx:   5,
							EnclosedFivePrimeIdx:  8,
							EnclosedThreePrimeIdx: 12,
							ClosingThreePrimeIdx:  15,
							Structures: []StemStructure{
								StemStructure{
									ClosingFivePrimeIdx:   5,
									EnclosedFivePrimeIdx:  6,
									EnclosedThreePrimeIdx: 14,
									ClosingThreePrimeIdx:  15,
									Type:                  StackingPair, Energy: 0,
								},
								StemStructure{
									ClosingFivePrimeIdx:   6,
									EnclosedFivePrimeIdx:  7,
									EnclosedThreePrimeIdx: 13,
									ClosingThreePrimeIdx:  14,
									Type:                  StackingPair, Energy: 0,
								},
								StemStructure{
									ClosingFivePrimeIdx:   7,
									EnclosedFivePrimeIdx:  8,
									EnclosedThreePrimeIdx: 12,
									ClosingThreePrimeIdx:  13,
									Type:                  StackingPair, Energy: 0,
								},
							},
							Energy: 0,
						},
						SingleStrandedFivePrimeIdx:  9,
						SingleStrandedThreePrimeIdx: 11,
						Energy:                      0,
					},
					&SingleStrandedRegion{
						FivePrimeIdx:  16,
						ThreePrimeIdx: 18,
					},
					&Hairpin{
						Stem: Stem{
							ClosingFivePrimeIdx:   19,
							EnclosedFivePrimeIdx:  20,
							EnclosedThreePrimeIdx: 29,
							ClosingThreePrimeIdx:  30,
							Structures: []StemStructure{
								StemStructure{
									ClosingFivePrimeIdx:   19,
									EnclosedFivePrimeIdx:  20,
									EnclosedThreePrimeIdx: 29,
									ClosingThreePrimeIdx:  30,
									Type:                  StackingPair, Energy: 0,
								},
							},
							Energy: 0,
						},
						SingleStrandedFivePrimeIdx:  21,
						SingleStrandedThreePrimeIdx: 28,
						Energy:                      0,
					},
					&SingleStrandedRegion{
						FivePrimeIdx:  31,
						ThreePrimeIdx: 32,
					},
				},
			},
			&SingleStrandedRegion{
				FivePrimeIdx:  34,
				ThreePrimeIdx: 34,
			},
		},
	}

	annotatedStructure, secondaryStructure, err := SecondaryStructureFromDotBracket(dotBracket)
	if err != nil {
		panic(err)
	}

	fmt.Println(annotatedStructure)
	// check if predicted and actual structures are same
	fmt.Println(isSecondaryStructuresEqual(actualSecondaryStructure, *secondaryStructure))

	// Output:
	// ee(mm((((hhh))))mmm((hhhhhhhh))mm)e
	// true
}

func ExampleSecondaryStructureFromDotBracket_hairpinWithoutSingleStrandedRegion() {
	dotBracket := "..((())).."

	// dotBracket with index
	// . . ( ( ( ) ) ) . .
	// 0 1 2 3 4 5 6 7 8 9

	actualSecondaryStructure := SecondaryStructure{
		Length:              10,
		ExteriorLoopsEnergy: 0,
		Energy:              0,
		Structures: []interface{}{
			&SingleStrandedRegion{
				FivePrimeIdx:  0,
				ThreePrimeIdx: 1,
			},
			&Hairpin{
				Stem: Stem{
					ClosingFivePrimeIdx:   2,
					EnclosedFivePrimeIdx:  4,
					EnclosedThreePrimeIdx: 5,
					ClosingThreePrimeIdx:  7,
					Structures: []StemStructure{
						StemStructure{
							ClosingFivePrimeIdx:   2,
							EnclosedFivePrimeIdx:  3,
							EnclosedThreePrimeIdx: 6,
							ClosingThreePrimeIdx:  7,
							Type:                  StackingPair, Energy: 0,
						},
						StemStructure{
							ClosingFivePrimeIdx:   3,
							EnclosedFivePrimeIdx:  4,
							EnclosedThreePrimeIdx: 5,
							ClosingThreePrimeIdx:  6,
							Type:                  StackingPair, Energy: 0,
						},
					},
					Energy: 0,
				},
				SingleStrandedFivePrimeIdx:  -1,
				SingleStrandedThreePrimeIdx: -1,
				Energy:                      0,
			},
			&SingleStrandedRegion{
				FivePrimeIdx:  8,
				ThreePrimeIdx: 9,
			},
		},
	}

	annotatedStructure, secondaryStructure, err := SecondaryStructureFromDotBracket(dotBracket)
	if err != nil {
		panic(err)
	}

	fmt.Println(annotatedStructure)
	// check if predicted and actual structures are same
	fmt.Println(isSecondaryStructuresEqual(actualSecondaryStructure, *secondaryStructure))

	// Output:
	// ee((()))ee
	// true
}

func ExampleSecondaryStructureFromDotBracket_interiorBulge() {
	dotBracket := "(((.).....).)"

	// dotBracket with index
	// ( ( ( . ) . . . . .  )  .  )
	// 0 1 2 3 4 5 6 7 8 9 10 11 12

	actualSecondaryStructure := SecondaryStructure{
		Length:              13,
		ExteriorLoopsEnergy: 0,
		Energy:              0,
		Structures: []interface{}{
			&Hairpin{
				Stem: Stem{
					ClosingFivePrimeIdx:   0,
					EnclosedFivePrimeIdx:  2,
					EnclosedThreePrimeIdx: 4,
					ClosingThreePrimeIdx:  12,
					Structures: []StemStructure{
						StemStructure{
							ClosingFivePrimeIdx:   0,
							EnclosedFivePrimeIdx:  1,
							EnclosedThreePrimeIdx: 10,
							ClosingThreePrimeIdx:  12,
							Type:                  Bulge,
							Energy:                0,
						},
						StemStructure{
							ClosingFivePrimeIdx:   1,
							EnclosedFivePrimeIdx:  2,
							EnclosedThreePrimeIdx: 4,
							ClosingThreePrimeIdx:  10,
							Type:                  Bulge,
							Energy:                0,
						},
					},
					Energy: 0,
				},
				SingleStrandedFivePrimeIdx:  3,
				SingleStrandedThreePrimeIdx: 3,
				Energy:                      0,
			},
		},
	}

	annotatedStructure, secondaryStructure, err := SecondaryStructureFromDotBracket(dotBracket)
	if err != nil {
		panic(err)
	}

	fmt.Println(annotatedStructure)
	// check if predicted and actual structures are same
	fmt.Println(isSecondaryStructuresEqual(actualSecondaryStructure, *secondaryStructure))

	// Output:
	// (((h)iiiii)i)
	// true
}

func ExampleSecondaryStructureFromDotBracket_interior1xnLoops() {
	dotBracket := "(.(.(...(..).)..).)"

	// dotBracket with index
	// ( . ( . ( . . . ( .  .  )  .  )  .  .  )  .  )
	// 0 1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18

	actualSecondaryStructure := SecondaryStructure{
		Length:              19,
		ExteriorLoopsEnergy: 0,
		Energy:              0,
		Structures: []interface{}{
			&Hairpin{
				Stem: Stem{
					ClosingFivePrimeIdx:   0,
					EnclosedFivePrimeIdx:  8,
					EnclosedThreePrimeIdx: 11,
					ClosingThreePrimeIdx:  18,
					Structures: []StemStructure{
						StemStructure{
							ClosingFivePrimeIdx:   0,
							EnclosedFivePrimeIdx:  2,
							EnclosedThreePrimeIdx: 16,
							ClosingThreePrimeIdx:  18,
							Type:                  Interior1x1Loop,
							Energy:                0,
						},
						StemStructure{
							ClosingFivePrimeIdx:   2,
							EnclosedFivePrimeIdx:  4,
							EnclosedThreePrimeIdx: 13,
							ClosingThreePrimeIdx:  16,
							Type:                  Interior2x1Loop,
							Energy:                0,
						},
						StemStructure{
							ClosingFivePrimeIdx:   4,
							EnclosedFivePrimeIdx:  8,
							EnclosedThreePrimeIdx: 11,
							ClosingThreePrimeIdx:  13,
							Type:                  Interior1xnLoop,
							Energy:                0,
						},
					},
					Energy: 0,
				},
				SingleStrandedFivePrimeIdx:  9,
				SingleStrandedThreePrimeIdx: 10,
				Energy:                      0,
			},
		},
	}

	annotatedStructure, secondaryStructure, err := SecondaryStructureFromDotBracket(dotBracket)
	if err != nil {
		panic(err)
	}

	fmt.Println(annotatedStructure)
	// check if predicted and actual structures are same
	fmt.Println(isSecondaryStructuresEqual(actualSecondaryStructure, *secondaryStructure))

	// Output:
	// (i(i(iii(hh)i)ii)i)
	// true
}

func ExampleSecondaryStructureFromDotBracket_interior2xnLoops() {
	dotBracket := "(..(...(..(..)....)..)..)"

	// dotBracket with index
	// ( . . ( . . . ( . .  (  .  .  )  .  .  .  .  )  .  .  )  .  .  )
	// 0 1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22 23 24

	actualSecondaryStructure := SecondaryStructure{
		Length:              25,
		ExteriorLoopsEnergy: 0,
		Energy:              0,
		Structures: []interface{}{
			&Hairpin{
				Stem: Stem{
					ClosingFivePrimeIdx:   0,
					EnclosedFivePrimeIdx:  10,
					EnclosedThreePrimeIdx: 13,
					ClosingThreePrimeIdx:  24,
					Structures: []StemStructure{
						StemStructure{
							ClosingFivePrimeIdx:   0,
							EnclosedFivePrimeIdx:  3,
							EnclosedThreePrimeIdx: 21,
							ClosingThreePrimeIdx:  24,
							Type:                  Interior2x2Loop,
							Energy:                0,
						},
						StemStructure{
							ClosingFivePrimeIdx:   3,
							EnclosedFivePrimeIdx:  7,
							EnclosedThreePrimeIdx: 18,
							ClosingThreePrimeIdx:  21,
							Type:                  Interior2x3Loop,
							Energy:                0,
						},
						StemStructure{
							ClosingFivePrimeIdx:   7,
							EnclosedFivePrimeIdx:  10,
							EnclosedThreePrimeIdx: 13,
							ClosingThreePrimeIdx:  18,
							Type:                  GenericInteriorLoop,
							Energy:                0,
						},
					},
					Energy: 0,
				},
				SingleStrandedFivePrimeIdx:  11,
				SingleStrandedThreePrimeIdx: 12,
				Energy:                      0,
			},
		},
	}

	annotatedStructure, secondaryStructure, err := SecondaryStructureFromDotBracket(dotBracket)
	if err != nil {
		panic(err)
	}

	fmt.Println(annotatedStructure)
	// check if predicted and actual structures are same
	fmt.Println(isSecondaryStructuresEqual(actualSecondaryStructure, *secondaryStructure))

	// Output:
	// (ii(iii(ii(hh)iiii)ii)ii)
	// true
}

func ExampleSecondaryStructureFromDotBracket_genericInteriorLoops() {
	dotBracket := "(..(...(....().....)...)....)"

	// dotBracket with index
	// ( . . ( . . . ( . .  .  .  (  )  .  .  .  .  .  )  .  .  .  )  .  .  .  .  )
	// 0 1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22 23 24 25 26 27 28

	actualSecondaryStructure := SecondaryStructure{
		Length:              29,
		ExteriorLoopsEnergy: 0,
		Energy:              0,
		Structures: []interface{}{
			&Hairpin{
				Stem: Stem{
					ClosingFivePrimeIdx:   0,
					EnclosedFivePrimeIdx:  12,
					EnclosedThreePrimeIdx: 13,
					ClosingThreePrimeIdx:  28,
					Structures: []StemStructure{
						StemStructure{
							ClosingFivePrimeIdx:   0,
							EnclosedFivePrimeIdx:  3,
							EnclosedThreePrimeIdx: 23,
							ClosingThreePrimeIdx:  28,
							Type:                  GenericInteriorLoop,
							Energy:                0,
						},
						StemStructure{
							ClosingFivePrimeIdx:   3,
							EnclosedFivePrimeIdx:  7,
							EnclosedThreePrimeIdx: 19,
							ClosingThreePrimeIdx:  23,
							Type:                  GenericInteriorLoop,
							Energy:                0,
						},
						StemStructure{
							ClosingFivePrimeIdx:   7,
							EnclosedFivePrimeIdx:  12,
							EnclosedThreePrimeIdx: 13,
							ClosingThreePrimeIdx:  19,
							Type:                  GenericInteriorLoop,
							Energy:                0,
						},
					},
					Energy: 0,
				},
				SingleStrandedFivePrimeIdx:  -1,
				SingleStrandedThreePrimeIdx: -1,
				Energy:                      0,
			},
		},
	}

	annotatedStructure, secondaryStructure, err := SecondaryStructureFromDotBracket(dotBracket)
	if err != nil {
		panic(err)
	}

	fmt.Println(annotatedStructure)
	// check if predicted and actual structures are same
	fmt.Println(isSecondaryStructuresEqual(actualSecondaryStructure, *secondaryStructure))

	// Output:
	// (ii(iii(iiii()iiiii)iii)iiii)
	// true
}

func ExampleSecondaryStructureFromDotBracket_stemWithoutEnclosedPair() {
	dotBracket := ".(.)."

	// dotBracket with index
	// . ( . ) .
	// 0 1 2 3 4

	actualSecondaryStructure := SecondaryStructure{
		Length:              5,
		ExteriorLoopsEnergy: 0,
		Energy:              0,
		Structures: []interface{}{
			SingleStrandedRegion{
				FivePrimeIdx:  0,
				ThreePrimeIdx: 0,
			},
			Hairpin{
				Stem: Stem{
					ClosingFivePrimeIdx:   1,
					EnclosedFivePrimeIdx:  -1,
					EnclosedThreePrimeIdx: -1,
					ClosingThreePrimeIdx:  3,
					Structures:            []StemStructure{},
					Energy:                0,
				},
				SingleStrandedFivePrimeIdx:  2,
				SingleStrandedThreePrimeIdx: 2,
				Energy:                      0,
			},
			SingleStrandedRegion{
				FivePrimeIdx:  4,
				ThreePrimeIdx: 4,
			},
		},
	}

	annotatedStructure, secondaryStructure, err := SecondaryStructureFromDotBracket(dotBracket)
	if err != nil {
		panic(err)
	}

	fmt.Println(annotatedStructure)
	// check if predicted and actual structures are same
	fmt.Println(isSecondaryStructuresEqual(actualSecondaryStructure, *secondaryStructure))

	// Output:
	// e(h)e
	// true
}

func isSecondaryStructuresEqual(structA, structB SecondaryStructure) bool {
	// convert predicted and actual structures to their string representation
	// so that they can be compared
	if structA.Length != structB.Length || structA.ExteriorLoopsEnergy != structB.ExteriorLoopsEnergy {
		return false
	}

	return isStructuresEqual(structA.Structures, structB.Structures)
}

func isStructuresEqual(structAStructures, structBStructures []interface{}) bool {
	if len(structAStructures) != len(structBStructures) {
		return false
	}
	for i := 0; i < len(structAStructures); i++ {
		structAStruct, structBStruct := structAStructures[i], structBStructures[i]

		switch structAStruct := structAStruct.(type) {
		case *SingleStrandedRegion:
			structAStructStr := fmt.Sprintf("%+v", *structAStruct)
			structBStructStr := fmt.Sprintf("%+v", *structBStruct.(*SingleStrandedRegion))
			if structAStructStr != structBStructStr {
				return false
			}
		case *Hairpin:
			structAStructStr := fmt.Sprintf("%+v", *structAStruct)
			structBStructStr := fmt.Sprintf("%+v", *structBStruct.(*Hairpin))
			if structAStructStr != structBStructStr {
				return false
			}
		case *MultiLoop:
			structAStructStemStr := fmt.Sprintf("%+v", structAStruct.Stem)
			structBStructStemStr := fmt.Sprintf("%+v", structBStruct.(*MultiLoop).Stem)
			if structAStructStemStr != structBStructStemStr {
				return false
			}
			if !isStructuresEqual(structAStruct.Substructures, structBStruct.(*MultiLoop).Substructures) {
				return false
			}
		}
	}
	return true
}
