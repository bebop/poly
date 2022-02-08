package codon_test

import (
	"fmt"
	"io/ioutil"
	"os"
	"strings"

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

	sequence, _ := genbank.Read("../../data/puc19.gbk")
	codonTable := codon.GetCodonTable(11)

	// a string builder to build a single concatenated string of all coding regions
	var codingRegionsBuilder strings.Builder

	// iterate through the features of the genbank file and if the feature is a coding region, append the sequence to the string builder
	for _, feature := range sequence.Features {
		if feature.Type == "CDS" {
			sequence, _ := feature.GetSequence()
			codingRegionsBuilder.WriteString(sequence)
		}
	}

	// get the concatenated sequence string of the coding regions
	codingRegions := codingRegionsBuilder.String()

	// weight our codon optimization table using the regions we collected from the genbank file above
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
	sequence, _ := genbank.Read("../../data/puc19.gbk")
	codonTable := codon.GetCodonTable(11)

	// a string builder to build a single concatenated string of all coding regions
	var codingRegionsBuilder strings.Builder

	// iterate through the features of the genbank file and if the feature is a coding region, append the sequence to the string builder
	for _, feature := range sequence.Features {
		if feature.Type == "CDS" {
			sequence, _ := feature.GetSequence()
			codingRegionsBuilder.WriteString(sequence)
		}
	}

	// get the concatenated sequence string of the coding regions
	codingRegions := codingRegionsBuilder.String()

	// weight our codon optimization table using the regions we collected from the genbank file above
	optimizationTable := codonTable.OptimizeTable(codingRegions)

	sequence2, _ := genbank.Read("../../data/phix174.gb")
	codonTable2 := codon.GetCodonTable(11)
	// a string builder to build a single concatenated string of all coding regions
	var codingRegionsBuilder2 strings.Builder

	// iterate through the features of the genbank file and if the feature is a coding region, append the sequence to the string builder
	for _, feature := range sequence2.Features {
		if feature.Type == "CDS" {
			sequence, _ := feature.GetSequence()
			codingRegionsBuilder.WriteString(sequence)
		}
	}

	// get the concatenated sequence string of the coding regions
	codingRegions2 := codingRegionsBuilder2.String()

	// weight our codon optimization table using the regions we collected from the genbank file above
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
	sequence, _ := genbank.Read("../../data/puc19.gbk")
	codonTable := codon.GetCodonTable(11)

	// a string builder to build a single concatenated string of all coding regions
	var codingRegionsBuilder strings.Builder

	// iterate through the features of the genbank file and if the feature is a coding region, append the sequence to the string builder
	for _, feature := range sequence.Features {
		if feature.Type == "CDS" {
			sequence, _ := feature.GetSequence()
			codingRegionsBuilder.WriteString(sequence)
		}
	}

	// get the concatenated sequence string of the coding regions
	codingRegions := codingRegionsBuilder.String()

	// weight our codon optimization table using the regions we collected from the genbank file above
	optimizationTable := codonTable.OptimizeTable(codingRegions)

	sequence2, _ := genbank.Read("../../data/phix174.gb")
	codonTable2 := codon.GetCodonTable(11)

	// a string builder to build a single concatenated string of all coding regions
	var codingRegionsBuilder2 strings.Builder

	// iterate through the features of the genbank file and if the feature is a coding region, append the sequence to the string builder
	for _, feature := range sequence2.Features {
		if feature.Type == "CDS" {
			sequence, _ := feature.GetSequence()
			codingRegionsBuilder.WriteString(sequence)
		}
	}

	// get the concatenated sequence string of the coding regions
	codingRegions2 := codingRegionsBuilder2.String()

	// weight our codon optimization table using the regions we collected from the genbank file above
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
