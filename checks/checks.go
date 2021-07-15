package checks

import "github.com/Open-Science-Global/poly/transform"

// IsPalindromic accepts a sequence of even length and returns if it is
// palindromic. More here - https://en.wikipedia.org/wiki/Palindromic_sequence
func IsPalindromic(sequence string) bool {
	return sequence == transform.ReverseComplement(sequence)
}

func gcContent(sequence string) float64 {

	if len(sequence) < 1 {
	  return 0.0
	}
  
	const c = 67
	const g = 71
	count := 0.0
  
	for i, character := range(sequence) {
	  if character == g || character == c {
		count += 1.0 
	  }
	}
  
	return count / float64(len(sequence))
  }
