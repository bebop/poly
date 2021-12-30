package checks_test

import (
	"testing"

	"github.com/TimothyStiles/poly/checks"
)

// This also needs an example test.
func TestIsPalindromic(t *testing.T) {
	ecori := checks.IsPalindromic("GAATTC")
	if ecori != true {
		t.Errorf("IsPalindromic failed to call EcoRI a palindrome")
	}
	bsai := checks.IsPalindromic("GGTCTC")
	if bsai != false {
		t.Errorf("IsPalindromic failed call BsaI NOT a palindrome")
	}
}

func TestGcContent(t *testing.T) {
	content := checks.GcContent("GGTATC")
	if content != 0.5 {
		t.Errorf("GcContent did not properly calculate GC content")
	}
}
