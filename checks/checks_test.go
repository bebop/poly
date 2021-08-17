package checks

import "testing"

// This also needs an example test.
func TestIsPalindromic(t *testing.T) {
	ecori := IsPalindromic("GAATTC")
	if ecori != true {
		t.Errorf("IsPalindromic failed to call EcoRI a palindrome")
	}
	bsai := IsPalindromic("GGTCTC")
	if bsai != false {
		t.Errorf("IsPalindromic failed call BsaI NOT a palindrome")
	}
}

func TestGcContent(t *testing.T) {
	content := GcContent("GGTATC")
	if content != 0.5 {
		t.Errorf("GcContent did not properly calculate GC content")
	}
}
