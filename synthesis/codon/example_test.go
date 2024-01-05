package codon_test

import (
	"fmt"
	"os"

	"github.com/bebop/poly/io/genbank"
	"github.com/bebop/poly/synthesis/codon"
)

func ExampleTranslationTable_Translate() {
	gfpTranslation := "MASKGEELFTGVVPILVELDGDVNGHKFSVSGEGEGDATYGKLTLKFICTTGKLPVPWPTLVTTFSYGVQCFSRYPDHMKRHDFFKSAMPEGYVQERTISFKDDGNYKTRAEVKFEGDTLVNRIELKGIDFKEDGNILGHKLEYNYNSHNVYITADKQKNGIKANFKIRHNIEDGSVQLADHYQQNTPIGDGPVLLPDNHYLSTQSALSKDPNEKRDHMVLLEFVTAAGITHGMDELYK*"
	gfpDnaSequence := "ATGGCTAGCAAAGGAGAAGAACTTTTCACTGGAGTTGTCCCAATTCTTGTTGAATTAGATGGTGATGTTAATGGGCACAAATTTTCTGTCAGTGGAGAGGGTGAAGGTGATGCTACATACGGAAAGCTTACCCTTAAATTTATTTGCACTACTGGAAAACTACCTGTTCCATGGCCAACACTTGTCACTACTTTCTCTTATGGTGTTCAATGCTTTTCCCGTTATCCGGATCATATGAAACGGCATGACTTTTTCAAGAGTGCCATGCCCGAAGGTTATGTACAGGAACGCACTATATCTTTCAAAGATGACGGGAACTACAAGACGCGTGCTGAAGTCAAGTTTGAAGGTGATACCCTTGTTAATCGTATCGAGTTAAAAGGTATTGATTTTAAAGAAGATGGAAACATTCTCGGACACAAACTCGAGTACAACTATAACTCACACAATGTATACATCACGGCAGACAAACAAAAGAATGGAATCAAAGCTAACTTCAAAATTCGCCACAACATTGAAGATGGATCCGTTCAACTAGCAGACCATTATCAACAAAATACTCCAATTGGCGATGGCCCTGTCCTTTTACCAGACAACCATTACCTGTCGACACAATCTGCCCTTTCGAAAGATCCCAACGAAAAGCGTGACCACATGGTCCTTCTTGAGTTTGTAACTGCTGCTGGGATTACACATGGCATGGATGAGCTCTACAAATAA"
	table, err := codon.NewTranslationTable(11)
	if err != nil {
		fmt.Printf("error running example: %s\n", err)
		return
	}

	testTranslation, _ := table.Translate(gfpDnaSequence) // need to specify which codons map to which amino acids per NCBI table

	fmt.Println(gfpTranslation == testTranslation)
	// output: true
}

func ExampleTranslationTable_UpdateWeights() {
	gfpTranslation := "MASKGEELFTGVVPILVELDGDVNGHKFSVSGEGEGDATYGKLTLKFICTTGKLPVPWPTLVTTFSYGVQCFSRYPDHMKRHDFFKSAMPEGYVQERTISFKDDGNYKTRAEVKFEGDTLVNRIELKGIDFKEDGNILGHKLEYNYNSHNVYITADKQKNGIKANFKIRHNIEDGSVQLADHYQQNTPIGDGPVLLPDNHYLSTQSALSKDPNEKRDHMVLLEFVTAAGITHGMDELYK*"

	sequenceWithCustomWeights := "ATGGCGAGCAAGGGCGAAGAGCTTTTTACTGGAGTGGTACCCATCCTTGTGGAGCTGGATGGGGATGTTAATGGGCACAAGTTTTCTGTGTCCGGTGAGGGGGAGGGTGACGCGACCTATGGCAAACTAACGTTGAAGTTTATCTGCACCACCGGCAAGCTCCCTGTCCCTTGGCCGACGCTGGTAACCACTTTTTCATACGGAGTGCAATGCTTTTCACGATACCCAGACCACATGAAACGGCACGACTTCTTCAAGAGCGCGATGCCAGAAGGTTATGTGCAAGAGCGTACGATCTCATTCAAGGACGACGGGAATTATAAGACAAGAGCAGAGGTGAAATTTGAGGGGGACACGTTAGTAAATCGGATTGAATTAAAGGGAATCGACTTTAAGGAGGATGGGAACATACTTGGTCACAAACTGGAATATAATTACAATTCACACAATGTTTACATCACTGCCGACAAGCAAAAAAATGGGATTAAAGCAAATTTCAAAATTCGGCATAATATTGAGGATGGTAGTGTCCAGCTCGCGGATCACTATCAGCAAAACACACCTATCGGAGACGGACCCGTTTTACTACCGGATAATCATTACTTAAGCACCCAATCAGCGTTATCCAAAGATCCGAACGAAAAACGTGACCACATGGTTCTCTTGGAGTTCGTCACCGCAGCTGGAATAACTCATGGAATGGACGAACTATACAAATAA"

	table, err := codon.NewTranslationTable(11)
	if err != nil {
		fmt.Printf("error running example: %s\n", err)
		return
	}

	// this example is using custom weights for different codons for Arginine. Use this if you would rather use your own
	// codon weights, they can also be computed for you with `UpdateWeightsWithSequence`.

	err = table.UpdateWeights([]codon.AminoAcid{
		{
			Letter: "R",
			Codons: []codon.Codon{
				{
					Triplet: "CGU",
					Weight:  1,
				},
				{
					Triplet: "CGA",
					Weight:  2,
				},
				{
					Triplet: "CGG",
					Weight:  4,
				},
				{
					Triplet: "AGA",
					Weight:  6,
				},
				{
					Triplet: "AGG",
					Weight:  2,
				},
			},
		},
	})
	if err != nil {
		fmt.Println("Could not update weights in example")
	}

	optimizedSequence, err := table.Optimize(gfpTranslation, 1)
	if err != nil {
		fmt.Printf("error running example: %s\n", err)
		return
	}

	fmt.Println(optimizedSequence == sequenceWithCustomWeights)
	// output: true
}

func ExampleTranslationTable_Optimize() {
	gfpTranslation := "MASKGEELFTGVVPILVELDGDVNGHKFSVSGEGEGDATYGKLTLKFICTTGKLPVPWPTLVTTFSYGVQCFSRYPDHMKRHDFFKSAMPEGYVQERTISFKDDGNYKTRAEVKFEGDTLVNRIELKGIDFKEDGNILGHKLEYNYNSHNVYITADKQKNGIKANFKIRHNIEDGSVQLADHYQQNTPIGDGPVLLPDNHYLSTQSALSKDPNEKRDHMVLLEFVTAAGITHGMDELYK*"

	sequence, _ := genbank.Read("../../data/puc19.gbk")
	codonTable, err := codon.NewTranslationTable(11)
	if err != nil {
		fmt.Printf("error running example: %s\n", err)
		return
	}

	_ = codonTable.UpdateWeightsWithSequence(sequence)

	// Here, we double check if the number of genes is equal to the number of stop codons
	stopCodonCount := 0
	for _, aa := range codonTable.AminoAcids {
		if aa.Letter == "*" {
			for _, codon := range aa.Codons {
				stopCodonCount = stopCodonCount + codon.Weight
			}
		}
	}

	if stopCodonCount != codonTable.Stats.GeneCount {
		fmt.Println("Stop codons don't equal number of genes!")
	}

	optimizedSequence, _ := codonTable.Optimize(gfpTranslation)
	optimizedSequenceTranslation, _ := codonTable.Translate(optimizedSequence)

	fmt.Println(optimizedSequenceTranslation == gfpTranslation)
	// output: true
}

func ExampleReadCodonJSON() {
	codontable := codon.ReadCodonJSON("../../data/bsub_codon_test.json")

	fmt.Println(codontable.GetWeightedAminoAcids()[0].Codons[0].Weight)
	//output: 28327
}

func ExampleParseCodonJSON() {
	file, _ := os.ReadFile("../../data/bsub_codon_test.json")
	codontable := codon.ParseCodonJSON(file)

	fmt.Println(codontable.GetWeightedAminoAcids()[0].Codons[0].Weight)
	//output: 28327
}

func ExampleWriteCodonJSON() {
	codontable := codon.ReadCodonJSON("../../data/bsub_codon_test.json")
	codon.WriteCodonJSON(codontable, "../../data/codon_test.json")
	testCodonTable := codon.ReadCodonJSON("../../data/codon_test.json")

	// cleaning up test data
	os.Remove("../../data/codon_test.json")

	fmt.Println(testCodonTable.GetWeightedAminoAcids()[0].Codons[0].Weight)
	//output: 28327
}

func ExampleCompromiseCodonTable() {
	sequence, _ := genbank.Read("../../data/puc19.gbk")

	// weight our codon optimization table using the regions we collected from the genbank file above
	optimizationTable, err := codon.NewTranslationTable(11)
	if err != nil {
		fmt.Printf("error running example: %s\n", err)
		return
	}

	err = optimizationTable.UpdateWeightsWithSequence(sequence)
	if err != nil {
		panic(fmt.Errorf("got unexpected error in an example: %w", err))
	}

	sequence2, _ := genbank.Read("../../data/phix174.gb")
	optimizationTable2, err := codon.NewTranslationTable(11)
	if err != nil {
		fmt.Printf("error running example: %s\n", err)
		return
	}

	err = optimizationTable2.UpdateWeightsWithSequence(sequence2)
	if err != nil {
		panic(fmt.Errorf("got unexpected error in an example: %w", err))
	}

	finalTable, _ := codon.CompromiseCodonTable(optimizationTable, optimizationTable2, 0.1)
	for _, aa := range finalTable.GetWeightedAminoAcids() {
		for _, codon := range aa.Codons {
			if codon.Triplet == "TAA" {
				fmt.Println(codon.Weight)
			}
		}
	}
	//output: 3863
}

func ExampleAddCodonTable() {
	sequence, _ := genbank.Read("../../data/puc19.gbk")

	// weight our codon optimization table using the regions we collected from the genbank file above
	optimizationTable, err := codon.NewTranslationTable(11)
	if err != nil {
		fmt.Printf("error running example: %s\n", err)
		return
	}

	err = optimizationTable.UpdateWeightsWithSequence(sequence)
	if err != nil {
		panic(fmt.Errorf("got unexpected error in an example: %w", err))
	}

	sequence2, _ := genbank.Read("../../data/phix174.gb")
	optimizationTable2, err := codon.NewTranslationTable(11)
	if err != nil {
		fmt.Printf("error running example: %s\n", err)
		return
	}

	err = optimizationTable2.UpdateWeightsWithSequence(sequence2)
	if err != nil {
		panic(fmt.Errorf("got unexpected error in an example: %w", err))
	}

	finalTable, err := codon.AddCodonTable(optimizationTable, optimizationTable2)
	if err != nil {
		panic(fmt.Errorf("got error in adding codon table example: %w", err))
	}

	for _, aa := range finalTable.AminoAcids {
		for _, codon := range aa.Codons {
			if codon.Triplet == "GGC" {
				fmt.Println(codon.Weight)
			}
		}
	}
	//output: 51
}
