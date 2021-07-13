package secondary_structure

import (
	"fmt"
)

// dotBracket: . . ( ( ( ( . . . )  )  )  )  .  .  .  (  (  .  .  .  .  .  .  .  .  )  )  .  .
// 			index: 0 1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22 23 24 25 26 27 28 29
func ExampleFromDotBracket() {
	dotBracket := "..((((...))))...((........)).."

	actualSecondaryStructure := SecondaryStructure{
		Length:              30,
		ExteriorLoopsEnergy: 0,
		Energy:              0,
		Structures: []interface{}{
			SingleStrandedRegion{
				FivePrimeIdx:  0,
				ThreePrimeIdx: 1,
			},
			Hairpin{
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
			SingleStrandedRegion{
				FivePrimeIdx:  13,
				ThreePrimeIdx: 15,
			},
			Hairpin{
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
			SingleStrandedRegion{
				FivePrimeIdx:  28,
				ThreePrimeIdx: 29,
			},
		},
	}

	annotatedStructure, secondaryStructure, err := FromDotBracket(dotBracket)
	if err != nil {
		panic(err)
	}

	predictedSecondaryStructureStr := fmt.Sprintf("%+v", *secondaryStructure)
	actualSecondaryStructureStr := fmt.Sprintf("%+v", actualSecondaryStructure)

	fmt.Println(annotatedStructure)
	fmt.Println(predictedSecondaryStructureStr == actualSecondaryStructureStr)

	// Output:
	// ee((((hhh))))eee((hhhhhhhh))ee
	// true
}

// dotBracket: . . ( . . ( ( ( ( .  .  .  )  )  )  )  .  .  .  (  (  .  .  .  .  .  .  .  .  )  )  .  .  )  .
// 			index: 0 1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22 23 24 25 26 27 28 29 30 31 32 33 34
func ExampleFromDotBracket_MultiLoop() {
	dotBracket := "..(..((((...))))...((........))..)."

	actualSecondaryStructure := SecondaryStructure{
		Length:              35,
		ExteriorLoopsEnergy: 0,
		Energy:              0,
		Structures: []interface{}{
			SingleStrandedRegion{
				FivePrimeIdx:  0,
				ThreePrimeIdx: 1,
			},
			MultiLoop{
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
					SingleStrandedRegion{
						FivePrimeIdx:  3,
						ThreePrimeIdx: 4,
					},
					Hairpin{
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
					SingleStrandedRegion{
						FivePrimeIdx:  16,
						ThreePrimeIdx: 18,
					},
					Hairpin{
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
					SingleStrandedRegion{
						FivePrimeIdx:  31,
						ThreePrimeIdx: 32,
					},
				},
			},
			SingleStrandedRegion{
				FivePrimeIdx:  34,
				ThreePrimeIdx: 34,
			},
		},
	}

	annotatedStructure, secondaryStructure, err := FromDotBracket(dotBracket)
	if err != nil {
		panic(err)
	}

	predictedSecondaryStructureStr := fmt.Sprintf("%+v", *secondaryStructure)
	actualSecondaryStructureStr := fmt.Sprintf("%+v", actualSecondaryStructure)

	fmt.Println(annotatedStructure)
	fmt.Println(predictedSecondaryStructureStr == actualSecondaryStructureStr)

	// Output:
	// ee(mm((((hhh))))mmm((hhhhhhhh))mm)e
	// true
}

// dotBracket: . . ( ( ( ) ) ) . .
// 			index: 0 1 2 3 4 5 6 7 8 9
func ExampleFromDotBracket_HairpinWithoutSingleStrandedRegion() {
	dotBracket := "..((())).."

	actualSecondaryStructure := SecondaryStructure{
		Length:              10,
		ExteriorLoopsEnergy: 0,
		Energy:              0,
		Structures: []interface{}{
			SingleStrandedRegion{
				FivePrimeIdx:  0,
				ThreePrimeIdx: 1,
			},
			Hairpin{
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
			SingleStrandedRegion{
				FivePrimeIdx:  8,
				ThreePrimeIdx: 9,
			},
		},
	}

	annotatedStructure, secondaryStructure, err := FromDotBracket(dotBracket)
	if err != nil {
		panic(err)
	}

	predictedSecondaryStructureStr := fmt.Sprintf("%+v", *secondaryStructure)
	actualSecondaryStructureStr := fmt.Sprintf("%+v", actualSecondaryStructure)

	fmt.Println(annotatedStructure)
	fmt.Println(predictedSecondaryStructureStr == actualSecondaryStructureStr)

	// Output:
	// ee((()))ee
	// true
}

func ExampleFromDotBracket_MultiLoop2() {
	dotBracket := "(((((((((...((((((.........))))))........((((((.......))))))..)))))))))"
	annotatedStructure, _, err := FromDotBracket(dotBracket)
	if err != nil {
		panic(err)
	}
	fmt.Println(annotatedStructure)

	// Output:
	// (((((((((mmm((((((hhhhhhhhh))))))mmmmmmmm((((((hhhhhhh))))))mm)))))))))
}

func ExampleFromDotBracket_InteriorUnpairedNucleotides() {
	dotBracket := "((((.((.((......))))((((...))....)).))))"
	annotatedStructure, _, err := FromDotBracket(dotBracket)
	if err != nil {
		panic(err)
	}
	fmt.Println(annotatedStructure)

	// Output:
	// ((((m((i((hhhhhh))))((((hhh))iiii))m))))
}
