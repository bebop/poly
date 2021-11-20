package poly_test

import (
	"github.com/TimothyStiles/poly/io/genbank"
	"github.com/TimothyStiles/poly/synthesis/codon"
	"github.com/TimothyStiles/poly/synthesis/fix"
)

func Example_getting_started() {
	// So we're going to design a puc19 derived plasmid that can be transformed
	// into and expressed in good ol' DH5A e.Coli

	/*** First we need to read in our data ***/
	dh5a := genbank.Read("/data/dh5a.gbk")  //dh5a e.coli
	puc19 := genbank.Read("data/puc19.gbk") // puc19 vector

	// hard coding usually is bad practice so don't do this at home kiddos
	gfpDnaSequence := "ATGGCTAGCAAAGGAGAAGAACTTTTCACTGGAGTTGTCCCAATTCTTGTTGAATTAGATGGTGATGTTAATGGGCACAAATTTTCTGTCAGTGGAGAGGGTGAAGGTGATGCTACATACGGAAAGCTTACCCTTAAATTTATTTGCACTACTGGAAAACTACCTGTTCCATGGCCAACACTTGTCACTACTTTCTCTTATGGTGTTCAATGCTTTTCCCGTTATCCGGATCATATGAAACGGCATGACTTTTTCAAGAGTGCCATGCCCGAAGGTTATGTACAGGAACGCACTATATCTTTCAAAGATGACGGGAACTACAAGACGCGTGCTGAAGTCAAGTTTGAAGGTGATACCCTTGTTAATCGTATCGAGTTAAAAGGTATTGATTTTAAAGAAGATGGAAACATTCTCGGACACAAACTCGAGTACAACTATAACTCACACAATGTATACATCACGGCAGACAAACAAAAGAATGGAATCAAAGCTAACTTCAAAATTCGCCACAACATTGAAGATGGATCCGTTCAACTAGCAGACCATTATCAACAAAATACTCCAATTGGCGATGGCCCTGTCCTTTTACCAGACAACCATTACCTGTCGACACAATCTGCCCTTTCGAAAGATCCCAACGAAAAGCGTGACCACATGGTCCTTCTTGAGTTTGTAACTGCTGCTGGGATTACACATGGCATGGATGAGCTCTACAAATAA"

	/*** Now we need to optimize our gfp sequence for expression ***/

	// first we need to define our coding table via genbanks weird number index definition
	codonTable := codon.GetCodonTable(11)

	// this pulls all of dh5a's coding regions and join them into a single string.
	codingRegions := codon.GetCodingRegions(dh5a)

	// we weight our optimization table using codon frequency within the coding regions of DH5A
	optimizationTable := codonTable.OptimizeTable(codingRegions)

	// we translate our DNA sequence into an amino acid sequence that our optimizer can work with
	gfpTranslation, _ := codon.Translate(gfpDnaSequence, codonTable)

	// we use our new optimization table to optimize our GFP sequence for expression in DH5A
	codonOptimizedSequence, _ := codon.Optimize(gfpTranslation, optimizationTable)

	// optimize our CDS for synthesis
	optimizedSequence, _, _ := fix.CdsSimple(codonOptimizedSequence, codonTable, []string{})

}
