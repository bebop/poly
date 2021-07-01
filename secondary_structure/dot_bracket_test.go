package secondary_structure

import "fmt"

func ExampleGetSecondaryStructure() {
	dotBracket := "..((((...))))...((........)).."
	annotatedStructure, secondaryStructure, err := FromDotBracket(dotBracket)
	if err != nil {
		panic(err)
	}
	fmt.Println(annotatedStructure)
	fmt.Println(dotBracketFromSecondaryStructure(*secondaryStructure) == dotBracket)
	// Output:
	// ee((((hhh))))eee((hhhhhhhh))ee
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
