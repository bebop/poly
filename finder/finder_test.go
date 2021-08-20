package finder

import (
	"fmt"
	"testing"

	"github.com/TimothyStiles/poly/io/genbank"
)

func ExampleForbiddenSequencesFinder() {
	sequence := "AAAAAATCGGTCGTAAAAAATT"
	var functions []func(string) []Match
	functions = append(functions, ForbiddenSequence([]string{"AAAAAA"}))

	matches := Find(sequence, functions)
	fmt.Println(matches)
	// Output: [{0 6 Forbidden sequence | AAAAAA} {14 20 Forbidden sequence | AAAAAA}]
}

func ExampleMatchSequences() {
	sequence := "AAAAAATCGGTCGTAAGGTCTCAAAATTGAGACC"
	var functions []func(string) []Match
	functions = append(functions, MatchSequences(map[string]string{"GGTCTC": "BsaI restriction binding site"}))

	matches := Find(sequence, functions)
	fmt.Println(matches)
	// Output: [{16 22 BsaI restriction binding site | GGTCTC} {28 34 BsaI restriction binding site | GAGACC}]
}

func ExampleRemoveRepeat() {
	sequence := "AAAAAATCGGTCGTAAGGTCTCAAAATTGAGACC"
	var functions []func(string) []Match
	functions = append(functions, RemoveRepeat(5))

	matches := Find(sequence, functions)
	fmt.Println(matches)
	// Output: [{1 6 Repeated sequence | AAAAA} {22 27 Repeated sequence | AAAAT}]

}

func ExampleGlobalRemoveRepeat() {
	sequence := "ATGAGTATTCAACATTTCCGTGTCGCCCTTATTCCCTTTTTTGCGGCATTTTGCCTTCCTGTTTTTGCTCACCCAGAAACGCTGGTGAAAGTAAAAGATGCTGAAGATCAGTTGGGTGCACGAGTGGGTTACATCGAACTGGATCTCAACAGCGGTAAGATCCTTGAGAGTTTTCGCCCCGAAGAACGTTTTCCAATGATGAGCACTTTTAAAGTTCTGCTATGTGGCGCGGTATTATCCCGTATTGACGCCGGGCAAGAGCAACTCGGTCGCCGCATACACTATTCTCAGAATGACTTGGTTGAGTACTCACCAGTCACAGAAAAGCATCTTACGGATGGCATGACAGTAAGAGAATTATGCAGTGCTGCCATAACCATGAGTGATAACACTGCGGCCAACTTACTTCTGACAACGATCGGAGGACCGAAGGAGCTAACCGCTTTTTTGCACAACATGGGGGATCATGTAACTCGCCTTGATCGTTGGGAACCGGAGCTGAATGAAGCCATACCAAACGACGAGCGTGACACCACGATGCCTGTAGCAATGGCAACAACGTTGCGCAAACTATTAACTGGCGAACTACTTACTCTAGCTTCCCGGCAACAATTAATAGACTGGATGGAGGCGGATAAAGTTGCAGGACCACTTCTGCGCTCGGCCCTTCCGGCTGGCTGGTTTATTGCTGATAAATCTGGAGCCGGTGAGCGTGGGTCTCGCGGTATCATTGCAGCACTGGGGCCAGATGGTAAGCCCTCCCGTATCGTAGTTATCTACACGACGGGGAGTCAGGCAACTATGGATGAACGAAATAGACAGATCGCTGAGATAGGTGCCTCACTGATTAAGCATTGGTAA"
	var functions []func(string) []Match
	functions = append(functions, GlobalRemoveRepeat(33, "ATGAGTATTCAACATTTCCGTGTCGCCCTTATT"))

	matches := Find(sequence, functions)
	fmt.Println(matches)
	// Output: [{0 33 Global repeated sequence | ATGAGTATTCAACATTTCCGTGTCGCCCTTATT}]
}

func TestAddMatchesToSequence(t *testing.T) {
	sequence := genbank.Read("../data/benchling.gb")
	var functions []func(string) []Match
	functions = append(functions, ForbiddenSequence([]string{"AAAAAA"}))

	matches := Find(sequence.Sequence, functions)
	sequenceWithMatchesAdded := AddMatchesToSequence(matches, sequence)

	if len(sequence.Features) == len(sequenceWithMatchesAdded.Features) {
		t.Errorf("Expected sequence to have different quantity of feature, but have same size.")
	}
	if len(sequenceWithMatchesAdded.Features) != 21 {
		t.Errorf("Expected sequence to have 21 feature, got: %d", len(sequence.Features))
	}
}
