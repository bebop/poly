package poly

import (
	"fmt"
	"sort"
	"strconv"
	"testing"
)

func TestGetSequenceMethods(t *testing.T) {

	gbk := ReadGbk("data/t4_intron.gb")

	// Check to see if GetSequence method works on Annotated struct
	if gbk.GetSequence() != gbk.Sequence {
		t.Errorf(" Sequence GetSequence method has failed'. Got this:\n%s instead of \n%s", gbk.GetSequence(), gbk.Sequence)
	}

	// Check to see if GetSequence method works on Features struct
	feature := gbk.Features[1].GetSequence()
	seq := "atgagattacaacgccagagcatcaaagattcagaagttagaggtaaatggtattttaatatcatcggtaaagattctgaacttgttgaaaaagctgaacatcttttacgtgatatgggatgggaagatgaatgcgatggatgtcctctttatgaagacggagaaagcgcaggattttggatttaccattctgacgtcgagcagtttaaagctgattggaaaattgtgaaaaagtctgtttga"
	if feature != seq {
		t.Errorf("Feature GetSequence method has failed. Got this:\n%s instead of \n%s", feature, seq)
	}

	// Check to see if GetSequence method works on Sequence struct
	if gbk.GetSequence() != gbk.Sequence {
		t.Errorf("Sequence GetSequence method has failed.. Got this:\n%s instead of \n%s", gbk.GetSequence(), gbk.Sequence)
	}

}

func TestLocationParser(t *testing.T) {
	gbk := ReadGbk("data/t4_intron.gb")

	// Read 1..243
	feature := gbk.Features[1].GetSequence()
	seq := "atgagattacaacgccagagcatcaaagattcagaagttagaggtaaatggtattttaatatcatcggtaaagattctgaacttgttgaaaaagctgaacatcttttacgtgatatgggatgggaagatgaatgcgatggatgtcctctttatgaagacggagaaagcgcaggattttggatttaccattctgacgtcgagcagtttaaagctgattggaaaattgtgaaaaagtctgtttga"
	if feature != seq {
		t.Errorf("Feature sequence parser has changed on test '1..243'. Got this:\n%s instead of \n%s", feature, seq)
	}

	// Read join(893..1441,2459..2770)
	featureJoin := gbk.Features[6].GetSequence()
	seqJoin := "atgaaacaataccaagatttaattaaagacatttttgaaaatggttatgaaaccgatgatcgtacaggcacaggaacaattgctctgttcggatctaaattacgctgggatttaactaaaggttttcctgcggtaacaactaagaagctcgcctggaaagcttgcattgctgagctaatatggtttttatcaggaagcacaaatgtcaatgatttacgattaattcaacacgattcgttaatccaaggcaaaacagtctgggatgaaaattacgaaaatcaagcaaaagatttaggataccatagcggtgaacttggtccaatttatggaaaacagtggcgtgattttggtggtgtagaccaaattatagaagttattgatcgtattaaaaaactgccaaatgataggcgtcaaattgtttctgcatggaatccagctgaacttaaatatatggcattaccgccttgtcatatgttctatcagtttaatgtgcgtaatggctatttggatttgcagtggtatcaacgctcagtagatgttttcttgggtctaccgtttaatattgcgtcatatgctacgttagttcatattgtagctaagatgtgtaatcttattccaggggatttgatattttctggtggtaatactcatatctatatgaatcacgtagaacaatgtaaagaaattttgaggcgtgaacctaaagagctttgtgagctggtaataagtggtctaccttataaattccgatatctttctactaaagaacaattaaaatatgttcttaaacttaggcctaaagatttcgttcttaacaactatgtatcacaccctcctattaaaggaaagatggcggtgtaa"
	if featureJoin != seqJoin {
		t.Errorf("Feature sequence parser has changed on test 'join(893..1441,2459..2770)'. Got this:\n%s instead of \n%s", featureJoin, seqJoin)
	}

	// Read complement(2791..3054)
	featureComplement := gbk.Features[10].GetSequence()
	seqComplement := "ttattcactacccggcatagacggcccacgctggaataattcgtcatattgtttttccgttaaaacagtaatatcgtagtaacagtcagaagaagttttaactgtggaaattttattatcaaaatactcacgagtcattttatgagtatagtattttttaccataaatggtaataggctgttctggtcctggaacttctaactcgcttgggttaggaagtgtaaaaagaactacaccagaagtatctttaaatcgtaaaatcat"
	if featureComplement != seqComplement {
		t.Errorf("Feature sequence parser has changed on test 'complement(2791..3054)'. Got this:\n%s instead of \n%s", featureComplement, seqComplement)
	}

	// Read join(complement(315..330),complement(339..896))
	// Note: it is known that some software, like Snapgene, assumes that since both strands are in the reverse direction
	// that the first sequence should be appended to the reverse sequence, instead of the second sequence
	// getting appended to the first. Biopython appends the second sequence to the first, and that is logically
	// the most obvious thing to do, so we are implementing it that way.
	featureJoinComplement := gbk.Features[3].GetSequence()
	seqJoinComplement := "ataccaatttaatcattcatttatatactgattccgtaagggttgttacttcatctattttataccaatgcgtttcaaccatttcacgcttgcttatatcatcaagaaaacttgcgtctaattgaactgttgaattaacacgatgccttttaacgatgcgagaaacaactacttcatctgcataaggtaatgcagcatataacagagcaggcccgccaattacacttactttagaattctgatcaagcatagtttcgaatggtgcattagggcttgacacttgaatttcgccgccagaaatgtaagttatatattgctcccaagtaatatagaaatgtgctaaatcgccgtctttagttacaggataatcacgcgcaaggtcacacaccacaatatggctacgaccaggaagtaatgtaggcaatgactggaacgttttagcacccataatcataattgtgccttcagtacgagctttaaaattctggaggtcctttttaactcgtccccatggtaaaccatcacctaaaccgaatgctaattcattaaagccgtcgaccgttttagttggaga"
	if featureJoinComplement != seqJoinComplement {
		t.Errorf("Feature sequence parser has changed on test 'join(complement(315..330),complement(339..896))'. Got this:\n%s instead of \n%s", featureJoinComplement, seqJoinComplement)
	}

	// Read complement(join(893..1098,1101..2770))
	featureComplementJoin := gbk.Features[5].GetSequence()
	seqComplementJoin := "ttacaccgccatctttcctttaataggagggtgtgatacatagttgttaagaacgaaatctttaggcctaagtttaagaacatattttaattgttctttagtagaaagatatcggaatttataaggtagaccacttattaccagctcacaaagctctttaggttcacgcctcaaaatttctttacattgttctacgtgattcatatagatatgagtattaccaccagaaaatatcaaatcccctggaataagattacacatcttagctacaatatgaactaacgtagcatatgacgcaatattaaacggtagcattatgttcagataaggtcgttaatcttaccccggaattatatccagctgcatgtcaccatgcagagcagactatatctccaacttgttaaagcaagttgtctatcgtttcgagtcacttgaccctactccccaaagggatagtcgttaggcatttatgtagaaccaattccatttatcagattttacacgataagtaactaatccagacgaaattttaaaatgtctagctgcatctgctgcacaatcaaaaataaccccatcacatgaaatctttttaatattactaggctttttacctttcatcttttctgatattttagatttagttatgtctgaatgcttatgattaaagaatgaattattttcacctgaacgatttctgcatttactacaagtataagcagaagtttgtatgcgaacaccgcacttacaaaacttatgggtttctggattccaacgcccgtttttacttccgggtttactgtaaagagctttccgaccatcaggtccaagtttaagcatcttagctttaacagtttcagaacgtttcttaataatttcttcttttaatggatgcgtagaacatgtatcaccaaacgttgcatcagcaatattgtatccattaattttagaattaagctctttaatccaaaaattttctcgttcaataatcaaatctttctcatatggaatttcttccaaaatagaacattcaaacacattaccatgtttgttaaaagacctctgaagttttatagaagaatggcatcctttttctaaatctttaaaatgcctcttccatctcttttcaaaatctttagcacttcctacatatactttattgtttaaagtatttttaatctgataaattccgcttttcataaatacctctttaaatatagaagtatttattaaagggcaagtcctacaatttagcacgggattgtctactagagaggttccccgtttagatagattacaagtataagtcaccttatactcaggcctcaattaacccaagaaaacatctactgagcgttgataccactgcaaatccaaatagccattacgcacattaaactgatagaacatatgacaaggcggtaatgccatatatttaagttcagctggattccatgcagaaacaatttgacgcctatcatttggcagttttttaatacgatcaataacttctataatttggtctacaccaccaaaatcacgccactgttttccataaattggaccaagttcaccgctatggtatcctaaatcttttgcttgattttcgtaattttcatcccagactgttttgccttggattaacgaatcgtgttgaattaatcgtaaatcatacatttgtgcttcctgataaaaaccatattagctcagcaatgcaagctttccaggcgagcttcttagttgttaccgcaggaaaacctttagttaaatcccagcgtaatttagatccgaacagagcaattgttcctgtgcctgtacgatcatcggtttcataaccattttcaaaaatgtctttaattaaatcttggtattgtttcat"
	if featureComplementJoin != seqComplementJoin {
		t.Errorf("Feature sequence parser has changed on test 'complement(join(893..1098,1101..2770))'. Got this:\n%s instead of \n%s", featureComplementJoin, seqComplementJoin)
	}
}

func TestIUPAC(t *testing.T) {
	testSeq := "ATN"
	testVariants := []string{"ATG", "ATA", "ATT", "ATC"}
	testVariantsIUPAC := AllVariantsIUPAC(testSeq)

	sort.Strings(testVariants)
	sort.Strings(testVariantsIUPAC)

	for index := range testVariants {
		if testVariants[index] != testVariantsIUPAC[index] {
			t.Errorf("IUPAC variant has changed on test 'allIUPAC('ATN')'. Got this:\n%s instead of \n%s", testVariantsIUPAC, testVariants)
		}
	}
}

func ExampleRandomProteinSequence() {
	// RandomProteinSequence builds a Protein Sequence by only passing through arguments a length and a seed that will be use to generate a randomly the sequence. The length needs to be greater than two because every sequence already have a start and stop codon. Seed makes this test deterministic.
	randomProtein, _ := RandomProteinSequence(15, 2)
	fmt.Println(randomProtein)
	// Output: MHHPAFRMFNTMYG*
}

func TestRandomProteinSequence(t *testing.T) {
	const length = 10
	const seed = 2
	sequence, _ := RandomProteinSequence(length, seed)

	if sequence[0] != 'M' {
		t.Errorf("Random sequence doesn't have the correct initial aminoacid in sequence 'RandomSequence(10, 2)'. Got this: \n%s instead of \n%s", string(sequence[0]), "M")
	}

	if sequence[len(sequence)-1] != '*' {
		t.Errorf("Random sequence doesn't have correct last aminoacid in sequence 'RandomSequence(10, 2)'. Got this: \n%s instead of \n%s", string(sequence[len(sequence)-1]), "*")
	}

	if len(sequence) != length {
		t.Errorf("Random sequence doesn't have the sequence size equal parameter passed through n 'RandomSequence(10, 2)'. Got this: \n%s instead of \n%s", strconv.Itoa(len(sequence)), strconv.Itoa(length))
	}
}

// Write a new case of test when you have a n inferior than 3
func TestRandomProteinSequenceError(t *testing.T) {
	const length = 2
	const seed = 4
	sequence, _ := RandomProteinSequence(length, seed)

	if len(sequence) != 0 {
		t.Errorf("Random sequence must have sequence size equals 0 'RandomSequence(2, 4)'. Got this: \n%s instead of \n%s", strconv.Itoa(len(sequence)), strconv.Itoa(length))
	}
}
