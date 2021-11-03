package codon_test

import (
	"fmt"
	"io/ioutil"
	"os"

	"github.com/TimothyStiles/poly/io/genbank"
	"github.com/TimothyStiles/poly/synthesis/codon"
)

func ExampleTranslate() {

	gfpTranslation := "MASKGEELFTGVVPILVELDGDVNGHKFSVSGEGEGDATYGKLTLKFICTTGKLPVPWPTLVTTFSYGVQCFSRYPDHMKRHDFFKSAMPEGYVQERTISFKDDGNYKTRAEVKFEGDTLVNRIELKGIDFKEDGNILGHKLEYNYNSHNVYITADKQKNGIKANFKIRHNIEDGSVQLADHYQQNTPIGDGPVLLPDNHYLSTQSALSKDPNEKRDHMVLLEFVTAAGITHGMDELYK*"
	gfpDnaSequence := "ATGGCTAGCAAAGGAGAAGAACTTTTCACTGGAGTTGTCCCAATTCTTGTTGAATTAGATGGTGATGTTAATGGGCACAAATTTTCTGTCAGTGGAGAGGGTGAAGGTGATGCTACATACGGAAAGCTTACCCTTAAATTTATTTGCACTACTGGAAAACTACCTGTTCCATGGCCAACACTTGTCACTACTTTCTCTTATGGTGTTCAATGCTTTTCCCGTTATCCGGATCATATGAAACGGCATGACTTTTTCAAGAGTGCCATGCCCGAAGGTTATGTACAGGAACGCACTATATCTTTCAAAGATGACGGGAACTACAAGACGCGTGCTGAAGTCAAGTTTGAAGGTGATACCCTTGTTAATCGTATCGAGTTAAAAGGTATTGATTTTAAAGAAGATGGAAACATTCTCGGACACAAACTCGAGTACAACTATAACTCACACAATGTATACATCACGGCAGACAAACAAAAGAATGGAATCAAAGCTAACTTCAAAATTCGCCACAACATTGAAGATGGATCCGTTCAACTAGCAGACCATTATCAACAAAATACTCCAATTGGCGATGGCCCTGTCCTTTTACCAGACAACCATTACCTGTCGACACAATCTGCCCTTTCGAAAGATCCCAACGAAAAGCGTGACCACATGGTCCTTCTTGAGTTTGTAACTGCTGCTGGGATTACACATGGCATGGATGAGCTCTACAAATAA"
	testTranslation, _ := codon.Translate(gfpDnaSequence, codon.GetCodonTable(11)) // need to specify which codons map to which amino acids per NCBI table

	fmt.Println(gfpTranslation == testTranslation)
	// output: true
}

func ExampleOptimize() {

	gfpTranslation := "MASKGEELFTGVVPILVELDGDVNGHKFSVSGEGEGDATYGKLTLKFICTTGKLPVPWPTLVTTFSYGVQCFSRYPDHMKRHDFFKSAMPEGYVQERTISFKDDGNYKTRAEVKFEGDTLVNRIELKGIDFKEDGNILGHKLEYNYNSHNVYITADKQKNGIKANFKIRHNIEDGSVQLADHYQQNTPIGDGPVLLPDNHYLSTQSALSKDPNEKRDHMVLLEFVTAAGITHGMDELYK*"

	sequence := genbank.Read("../../data/puc19.gbk")
	codonTable := codon.GetCodonTable(11)
	codingRegions := codon.GetCodingRegions(sequence)

	optimizationTable := codonTable.OptimizeTable(codingRegions)

	optimizedSequence, _ := codon.Optimize(gfpTranslation, optimizationTable)
	optimizedSequenceTranslation, _ := codon.Translate(optimizedSequence, optimizationTable)

	fmt.Println(optimizedSequenceTranslation == gfpTranslation)
	// output: true
}

func ExampleGetCodingRegions() {

	gfpTranslation := "MASKGEELFTGVVPILVELDGDVNGHKFSVSGEGEGDATYGKLTLKFICTTGKLPVPWPTLVTTFSYGVQCFSRYPDHMKRHDFFKSAMPEGYVQERTISFKDDGNYKTRAEVKFEGDTLVNRIELKGIDFKEDGNILGHKLEYNYNSHNVYITADKQKNGIKANFKIRHNIEDGSVQLADHYQQNTPIGDGPVLLPDNHYLSTQSALSKDPNEKRDHMVLLEFVTAAGITHGMDELYK*"

	sequence := genbank.Read("../../data/puc19.gbk")
	codonTable := codon.GetCodonTable(11)

	// GetCodingRegions returns a single concatenated string of all coding regions.
	codingRegions := codon.GetCodingRegions(sequence)

	optimizationTable := codonTable.OptimizeTable(codingRegions)

	optimizedSequence, _ := codon.Optimize(gfpTranslation, optimizationTable)
	optimizedSequenceTranslation, _ := codon.Translate(optimizedSequence, optimizationTable)

	fmt.Println(optimizedSequenceTranslation == gfpTranslation)
	// output: true

}

func ExampleReadCodonJSON() {
	codontable := codon.ReadCodonJSON("../../data/bsub_codon_test.json")

	fmt.Println(codontable.AminoAcids[0].Codons[0].Weight)
	//output: 28327
}

func ExampleParseCodonJSON() {
	file, _ := ioutil.ReadFile("../../data/bsub_codon_test.json")
	codontable := codon.ParseCodonJSON(file)

	fmt.Println(codontable.AminoAcids[0].Codons[0].Weight)
	//output: 28327
}

func ExampleWriteCodonJSON() {
	codontable := codon.ReadCodonJSON("../../data/bsub_codon_test.json")
	codon.WriteCodonJSON(codontable, "../../data/codon_test.json")
	testCodonTable := codon.ReadCodonJSON("../../data/codon_test.json")

	// cleaning up test data
	os.Remove("../../data/codon_test.json")

	fmt.Println(testCodonTable.AminoAcids[0].Codons[0].Weight)
	//output: 28327
}

func ExampleCompromiseCodonTable() {
	sequence := genbank.Read("../../data/puc19.gbk")
	codonTable := codon.GetCodonTable(11)
	codingRegions := codon.GetCodingRegions(sequence)
	optimizationTable := codonTable.OptimizeTable(codingRegions)

	sequence2 := genbank.Read("../../data/phix174.gb")
	codonTable2 := codon.GetCodonTable(11)
	codingRegions2 := codon.GetCodingRegions(sequence2)
	optimizationTable2 := codonTable2.OptimizeTable(codingRegions2)

	finalTable, _ := codon.CompromiseCodonTable(optimizationTable, optimizationTable2, 0.1)
	for _, aa := range finalTable.AminoAcids {
		for _, codon := range aa.Codons {
			if codon.Triplet == "TAA" {
				fmt.Println(codon.Weight)
			}
		}
	}
	//output: 2727
}

func ExampleAddCodonTable() {
	sequence := genbank.Read("../../data/puc19.gbk")
	codonTable := codon.GetCodonTable(11)
	codingRegions := codon.GetCodingRegions(sequence)
	optimizationTable := codonTable.OptimizeTable(codingRegions)

	sequence2 := genbank.Read("../../data/phix174.gb")
	codonTable2 := codon.GetCodonTable(11)
	codingRegions2 := codon.GetCodingRegions(sequence2)
	optimizationTable2 := codonTable2.OptimizeTable(codingRegions2)

	finalTable := codon.AddCodonTable(optimizationTable, optimizationTable2)
	for _, aa := range finalTable.AminoAcids {
		for _, codon := range aa.Codons {
			if codon.Triplet == "GGC" {
				fmt.Println(codon.Weight)
			}
		}
	}
	//output: 90
}
