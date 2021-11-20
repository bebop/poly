package poly_test

import (
	"fmt"

	"github.com/TimothyStiles/poly/io/genbank"
	"github.com/TimothyStiles/poly/io/poly"
	"github.com/TimothyStiles/poly/primers/pcr"
	"github.com/TimothyStiles/poly/synthesis/codon"
)

func Example_getting_started() {
	// Here's how you'd optimize and insert in Poly to be expressed and synthesized

	/*** First we need to read in our data ***/
	bsub := genbank.Read("./data/bsub.gbk") //bsub

	// hard coding usually is bad practice so don't do this at home kiddos
	gfpDnaSequence := "ATGGCTAGCAAAGGAGAAGAACTTTTCACTGGAGTTGTCCCAATTCTTGTTGAATTAGATGGTGATGTTAATGGGCACAAATTTTCTGTCAGTGGAGAGGGTGAAGGTGATGCTACATACGGAAAGCTTACCCTTAAATTTATTTGCACTACTGGAAAACTACCTGTTCCATGGCCAACACTTGTCACTACTTTCTCTTATGGTGTTCAATGCTTTTCCCGTTATCCGGATCATATGAAACGGCATGACTTTTTCAAGAGTGCCATGCCCGAAGGTTATGTACAGGAACGCACTATATCTTTCAAAGATGACGGGAACTACAAGACGCGTGCTGAAGTCAAGTTTGAAGGTGATACCCTTGTTAATCGTATCGAGTTAAAAGGTATTGATTTTAAAGAAGATGGAAACATTCTCGGACACAAACTCGAGTACAACTATAACTCACACAATGTATACATCACGGCAGACAAACAAAAGAATGGAATCAAAGCTAACTTCAAAATTCGCCACAACATTGAAGATGGATCCGTTCAACTAGCAGACCATTATCAACAAAATACTCCAATTGGCGATGGCCCTGTCCTTTTACCAGACAACCATTACCTGTCGACACAATCTGCCCTTTCGAAAGATCCCAACGAAAAGCGTGACCACATGGTCCTTCTTGAGTTTGTAACTGCTGCTGGGATTACACATGGCATGGATGAGCTCTACAAATAA"

	/*** Now we need to optimize our gfp sequence for expression ***/

	// first we need to define our coding table via genbank's weird number index definition
	codonTable := codon.GetCodonTable(11)

	// this pulls all of bsub's coding regions and join them into a single string.
	codingRegions := codon.GetCodingRegions(bsub)

	// we weight our optimization table using codon frequency within the coding regions of bsub
	optimizationTable := codonTable.OptimizeTable(codingRegions)

	// we translate our DNA sequence into an amino acid sequence that our optimizer can work with
	gfpTranslation, _ := codon.Translate(gfpDnaSequence, codonTable)

	// we use our new optimization table to optimize our GFP sequence for expression in bsub
	codonOptimizedSequence, _ := codon.Optimize(gfpTranslation, optimizationTable)

	// since our optimizedSequence output is stochastic we'll translate it and compare to gfpTranslation
	// to see if everything worked correctly
	optimizedTranslation, _ := codon.Translate(codonOptimizedSequence, codonTable)
	fmt.Println(optimizedTranslation)
	// Output: MASKGEELFTGVVPILVELDGDVNGHKFSVSGEGEGDATYGKLTLKFICTTGKLPVPWPTLVTTFSYGVQCFSRYPDHMKRHDFFKSAMPEGYVQERTISFKDDGNYKTRAEVKFEGDTLVNRIELKGIDFKEDGNILGHKLEYNYNSHNVYITADKQKNGIKANFKIRHNIEDGSVQLADHYQQNTPIGDGPVLLPDNHYLSTQSALSKDPNEKRDHMVLLEFVTAAGITHGMDELYK*
}

func Example_cloning_out() {

	// In this example we're going to generate a primer library to clone out
	// every CDS in b.subtilis

	bsub := genbank.Read("./data/bsub.gbk")

	type cloneOut struct {
		CDS           poly.Feature
		Sequence      string
		FowardPrimer  string
		ReversePrimer string
	}

	var reactions []cloneOut // declaring our list of primers so we can append to it

	for _, feature := range bsub.Features {
		if feature.Type == "CDS" {
			var reaction cloneOut // initialize our reaction that will be appended

			// store the feature and its sequence in our reaction in case we need it later
			reaction.CDS = feature
			reaction.Sequence = feature.GetSequence()

			// generate forward and reverse primers and store it in our struct
			forward, reverse := pcr.DesignPrimers(reaction.Sequence, 56.0) // <- 56.0 is our melting temp. Again, don't hardcode at home kids.
			reaction.FowardPrimer = forward
			reaction.ReversePrimer = reverse

			// append our reaction to our reactions slice (slice is Go for essentially a list, or vector)
			reactions = append(reactions, reaction)
		}
	}
	fmt.Println(reactions[0].FowardPrimer)
	// Output: ATGGAAAATATATTAGACCTGTGGAACCA
}
