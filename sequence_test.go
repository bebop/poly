package poly

import (
	"fmt"
	"testing"

	"github.com/google/go-cmp/cmp"
	"github.com/google/go-cmp/cmp/cmpopts"
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

func TestSequenceMethods(t *testing.T) {
	sequence := ReadGbk("data/puc19.gbk")
	originalSequenceString := sequence.Sequence
	deletedSubSequence := sequence.GetRange(0, 10)

	sequenceLength := len(sequence.Sequence)
	sequence.deleteRange(0, 10)
	alteredSequenceLength := len(sequence.Sequence)
	sequenceLengthDiff := sequenceLength - alteredSequenceLength
	expectedLength := 10

	if sequenceLengthDiff != expectedLength {
		t.Errorf("deleteRange Sequence method has changed on test. Got length: \n%q instead of \n%q", sequenceLengthDiff, expectedLength)
	}

	sequence.insertSequence(0, deletedSubSequence)

	if sequence.Sequence != originalSequenceString {
		t.Errorf("sequence methods have changed on test Got this: \n%s instead of \n\n%s", sequence.Sequence, originalSequenceString)
	}

}

func TestDeleteMeta(t *testing.T) {
	sequence := ReadGbk("data/puc19.gbk")

	// backupFeatures := sequence.Features
	for _, feature := range sequence.Features {
		sequence.RemoveFeature(feature)
	}

	if len(sequence.Features) != 0 {
		t.Errorf("sequence methods have changed on test Got this: \n%d instead of \n%d", len(sequence.Features), 0)

	}
}

func TestSingleInsert(t *testing.T) {
	sequence := ReadGbk("./data/puc19.gbk")
	sequenceCopy := ReadGbk("./data/puc19.gbk")

	// feature := sequence.Features[4]
	featureCopy := sequenceCopy.Features[4]
	sequence.Delete(featureCopy)

	for _, feature := range sequence.Features {
		if Equals(feature, featureCopy) {
			t.Errorf("feature was not deleted from Features.")
		}
	}

	// checks to see if subsequence was deleted
	if diff := cmp.Diff(sequence.GetSequence(), sequenceCopy.GetSequence(), cmpopts.IgnoreUnexported(Feature{})); diff == "" {
		t.Errorf(" mismatch (-want +got):\n%s", diff)
	}

	sequence.Insert(featureCopy)

	// checks to see if a feature was added
	if len(sequence.Features) != len(sequenceCopy.Features) {
		t.Errorf("sequence feature slices have different lengths. original: %d \n copy: %d", len(sequence.Features), len(sequenceCopy.Features))
	}

	// checks to see if subsequence was reinserted
	if diff := cmp.Diff(sequence.GetSequence(), sequenceCopy.GetSequence(), cmpopts.IgnoreUnexported(Feature{})); diff != "" {
		t.Errorf(" mismatch (-want +got):\n%s", diff)
	}

	sequence.sortFeatures()
	sequenceCopy.sortFeatures()

	for index := range sequence.Features {
		original := sequence.Features[index]
		copy := sequenceCopy.Features[index]
		if !Equals(original, copy) {
			fmt.Println(index)
		}
	}

	// checks to see
	if diff := cmp.Diff(sequence, sequenceCopy, cmpopts.IgnoreUnexported(Feature{})); diff != "" {
		t.Errorf(" mismatch (-want +got):\n%s", diff)
	}
}

func TestPlasmidEdgeInDel(t *testing.T) {

	sequence := ReadGbk("./data/puc19.gbk")
	sequenceCopy := ReadGbk("./data/puc19.gbk")

	featureCopy := sequenceCopy.Features[20]
	primerBindCopy := sequenceCopy.Features[1]
	sequence.Delete(featureCopy)

	for _, feature := range sequence.Features {
		if Equals(feature, featureCopy) {
			t.Errorf("feature was not deleted from Features.")
		}
	}

	// checks to see if subsequence was deleted
	if diff := cmp.Diff(sequence.GetSequence(), sequenceCopy.GetSequence(), cmpopts.IgnoreUnexported(Feature{})); diff == "" {
		t.Errorf(" mismatch (-want +got):\n%s", diff)
	}

	sequence.Insert(featureCopy)
	sequence.Annotate(primerBindCopy)

	// checks to see if a feature was added
	if len(sequence.Features) != len(sequenceCopy.Features) {
		t.Errorf("sequence feature slices have different lengths. original: %d \n copy: %d", len(sequence.Features), len(sequenceCopy.Features))
	}

	// checks to see if subsequence was reinserted
	if diff := cmp.Diff(sequence.GetSequence(), sequenceCopy.GetSequence(), cmpopts.IgnoreUnexported(Feature{})); diff != "" {
		t.Errorf(" mismatch (-want +got):\n%s", diff)
	}

	sequence.sortFeatures()
	sequenceCopy.sortFeatures()

	for index := range sequence.Features {
		original := sequence.Features[index]
		copy := sequenceCopy.Features[index]
		if !Equals(original, copy) {
			fmt.Println(index)
		}
	}

	// checks to see
	if diff := cmp.Diff(sequence, sequenceCopy, cmpopts.IgnoreUnexported(Feature{})); diff != "" {
		t.Errorf(" mismatch (-want +got):\n%s", diff)
	}

}
