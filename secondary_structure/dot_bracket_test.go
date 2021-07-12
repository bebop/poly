package secondary_structure

import "fmt"

func ExampleGetSecondaryStructure() {
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

	fmt.Println(annotatedStructure)
	predictedSecondaryStructureStr := fmt.Sprintf("%+v", *secondaryStructure)
	actualSecondaryStructureStr := fmt.Sprintf("%+v", actualSecondaryStructure)
	fmt.Println(predictedSecondaryStructureStr == actualSecondaryStructureStr)
	// Output:
	// ee((((hhh))))eee((hhhhhhhh))ee
	// true
}

func ExampleGetSecondaryStructureMultiLoop() {
	dotBracket := "..(..((((...))))...((........))..)."

	actualSecondaryStructure := SecondaryStructure{
		Length:              30,
		ExteriorLoopsEnergy: 0,
		Energy:              0,
		Structures: []interface{}{
			SingleStrandedRegion{
				FivePrimeIdx:  0,
				ThreePrimeIdx: 1,
			},
			MultiLoop{
				SubstructuresFivePrimeIdx:  0,
				SubstructuresThreePrimeIdx: -2,
				Energy:                     0,
				Stem: Stem{
					ClosingFivePrimeIdx: 2,
					EnclosedFivePrimeIdx:-1,
					EnclosedThreePrimeIdx:-1,
					ClosingThreePrimeIdx:33,
					Structures: []StemStructure{},
					Energy:0,
				},
				Substructures:
					SecondaryStructure{
						Length:              30,
						ExteriorLoopsEnergy: 0,
						Energy:              0,
						Structures: []interface{}{
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
		},
	}

	annotatedStructure, secondaryStructure, err := FromDotBracket(dotBracket)
	if err != nil {
		panic(err)
	}

	fmt.Println(annotatedStructure)
	predictedSecondaryStructureStr := fmt.Sprintf("%+v", *secondaryStructure)
	// actualSecondaryStructureStr := fmt.Sprintf("%v+", actualSecondaryStructure)
	// fmt.Println(predictedSecondaryStructureStr == actualSecondaryStructureStr)
	fmt.Println(predictedSecondaryStructureStr)
	// Output:
	// ee(mm((((hhh))))mmm((hhhhhhhh))mm)e
	// true
}

func ExampleGetSecondaryStructure2() {
	dotBracket := "(((((((((...((((((.........))))))........((((((.......))))))..)))))))))"
	annotatedStructure, secondaryStructure, err := FromDotBracket(dotBracket)
	if err != nil {
		panic(err)
	}
	fmt.Println(annotatedStructure)
	fmt.Println(dotBracketFromSecondaryStructure(*secondaryStructure) == dotBracket)
	// Output:
	// (((((((((mmm((((((hhhhhhhhh))))))mmmmmmmm((((((hhhhhhh))))))mm)))))))))
	// true
}

func ExampleGetSecondaryStructure3() {
	dotBracket := "((((.((((......))))((((...))....)).))))"
	annotatedStructure, secondaryStructure, err := FromDotBracket(dotBracket)
	if err != nil {
		panic(err)
	}
	fmt.Println(annotatedStructure)
	fmt.Println(dotBracketFromSecondaryStructure(*secondaryStructure) == dotBracket)
	// Output:
	// ((((m((((hhhhhh))))((((hhh))iiii))m))))
	// true
}

func ExampleGetSecondaryStructure4() {
	dotBracket := "......((((((.(((..(((((.(((....(((((......)))))..))).))))).)))....))))))..................."
	annotatedStructure, secondaryStructure, err := FromDotBracket(dotBracket)
	if err != nil {
		panic(err)
	}
	fmt.Println(annotatedStructure)
	fmt.Println(dotBracketFromSecondaryStructure(*secondaryStructure) == dotBracket)
	// Output:
	// eeeeee((((((i(((ii(((((i(((iiii(((((hhhhhh)))))ii)))i)))))i)))iiii))))))eeeeeeeeeeeeeeeeeee
	// true
}
