package secondary_structure

import "fmt"

func ExampleGetSecondaryStructure() {
	dotBracket := "..((((...))))...((........)).."
	annotatedStructure, secondaryStructure := FromDotBracket(dotBracket)
	fmt.Println(annotatedStructure)
	fmt.Println(DotBracket(secondaryStructure, 0) == dotBracket)
	// Output:
	// ee((((hhh))))eee((hhhhhhhh))ee
	// true
}

func ExampleGetSecondaryStructure2() {
	dotBracket := "(((((((((...((((((.........))))))........((((((.......))))))..)))))))))"
	annotatedStructure, secondaryStructure := FromDotBracket(dotBracket)
	fmt.Println(annotatedStructure)
	fmt.Println(DotBracket(secondaryStructure, 0) == dotBracket)
	// Output:
	// (((((((((mmm((((((hhhhhhhhh))))))mmmmmmmm((((((hhhhhhh))))))mm)))))))))
	// true
}

func ExampleGetSecondaryStructure3() {
	dotBracket := "((((.((((......))))((((...))....)).))))"
	annotatedStructure, secondaryStructure := FromDotBracket(dotBracket)
	fmt.Println(annotatedStructure)
	fmt.Println(DotBracket(secondaryStructure, 0) == dotBracket)
	// Output:
	// ((((m((((hhhhhh))))((((hhh))iiii))m))))
	// true
}

func ExampleGetSecondaryStructure4() {
	dotBracket := "......((((((.(((..(((((.(((....(((((......)))))..))).))))).)))....))))))..................."
	annotatedStructure, secondaryStructure := FromDotBracket(dotBracket)
	fmt.Println(annotatedStructure)
	fmt.Println(DotBracket(secondaryStructure, 0) == dotBracket)
	// Output:
	// eeeeee((((((i(((ii(((((i(((iiii(((((hhhhhh)))))ii)))i)))))i)))iiii))))))eeeeeeeeeeeeeeeeeee
	// true
}
