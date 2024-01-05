package fold

import (
	"fmt"
	"math"
	"testing"
)

func TestResult_MinimumFreeEnergy_LengthZero(t *testing.T) {
	result := Result{} // Create a Result instance with empty structs

	expectedEnergy := math.Inf(1)
	actualEnergy := result.MinimumFreeEnergy()

	if actualEnergy != expectedEnergy {
		t.Errorf("expected energy to be %f, but got %f", expectedEnergy, actualEnergy)
	}
}

func TestResult_DotBracket_LengthZero(t *testing.T) {
	result := Result{} // Create a Result instance with empty structs

	expectedDotBracket := ""
	actualDotBracket := result.DotBracket()

	if actualDotBracket != expectedDotBracket {
		t.Errorf("expected dot bracket to be %s, but got %s", expectedDotBracket, actualDotBracket)
	}
}

func TestNewFoldingContext_InvalidSequence(t *testing.T) {
	seq := "XYZ"
	temp := 37.0

	_, err := newFoldingContext(seq, temp)
	if err == nil {
		t.Errorf("expected error, but got nil")
	}
	expectedError := fmt.Errorf("the sequence %s is not RNA or DNA", seq)
	if err.Error() != expectedError.Error() {
		t.Errorf("expected error message to be %q, but got %q", expectedError.Error(), err.Error())
	}
}
