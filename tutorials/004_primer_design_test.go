package tutorials_test

import (
	"fmt"
	"log"
	"testing"

	"github.com/bebop/poly/io/genbank"
	"github.com/bebop/poly/primers/pcr"
)

/******************************************************************************
Sep, 12, 2022

== Designing Primers for Just About Anything ==

"Primers are short sequences of DNA that can be used to amplify DNA sequences
and they are the workhorse of modern molecular biology.

Essentially primers are short pieces of single stranded DNA that can
bind to a target sequence of single stranded DNA. These primers serve as a
marker for polymerases (the enzyme Poly is named after!) to bind and start adding
free floating nucleotides (ACTGs) to a single strand piece of DNA to form a
double stranded piece of DNA.

This is a crucial step in the process of PCR (polymerase chain reaction).
https://en.wikipedia.org/wiki/Polymerase_chain_reaction

You can read more about that at the link above but just know that an absolute huge
number of protocols from diagnostics to plasmid cloning use these primers so they're
super important."

- From the Poly's primer design package level documentation


Primers are the workhorse of modern molecular biology. They're involved in almost
every molecular biology experiment and are a crucial component in modern lab
diagnostics.

In This tutorial we're going to design a large set of primers to help us extract and
isolate every protein coding region in the bacillus subtilis genome so we can express
and characterize each protein's structure individually in vivo (in the lab).

We could also use these primers for RNA interference experiments to supress
protein expression in vivo to better understand how these proteins interact.

Point is that primers are incredibly versatile tools by automating their design
we can do some pretty interesting (and even potentially lucrative) experiments.

TTFN,
Tim
******************************************************************************/

// if you're using VS-CODE you should see a DEBUG TEST button right below this
// comment. Please set break points and use them early and often.
func TestPrimersTutorial(t *testing.T) {

	// This is a struct that we'll use to store the results of our primer designs for cloning out each gene
	type CloneOut struct {
		CDS           genbank.Feature
		Sequence      string
		ForwardPrimer string
		ReversePrimer string
	}

	var reactions []CloneOut // <- declaring our list of primers so we can append to it

	// First let's get our annotated bacillus subtillus genome
	bsub, err := genbank.Read("../data/bsub.gbk")

	if err != nil {
		log.Fatal(err)
	}

	// For each feature in the genome we're going to design a primer pair if the feature is a coding sequence
	for _, feature := range bsub.Features {
		if feature.Type == "CDS" { // CDS stands for coding sequence (which means that it codes for a protein in this case)

			var reaction CloneOut // initialize our reaction that will be appended

			// store the feature and its sequence in our reaction in case we need it later
			reaction.CDS = feature
			reaction.Sequence, _ = feature.GetSequence()

			// generate forward and reverse primers and store it in our struct
			forward, reverse := pcr.DesignPrimers(reaction.Sequence, 56.0) // <- 56.0 is our melting temp. Again. don't hardcode values like this in real life. Put it in a constant or something.
			reaction.ForwardPrimer = forward
			reaction.ReversePrimer = reverse

			// append our reaction to a our reactions slice (slice is essentially Go's version of a list, or vector)
			reactions = append(reactions, reaction)
		}
	}

	// Check to see if any primer pairs were created and report the total
	if len(reactions) < 1 {
		t.Errorf("no reactions were created")
	} else {
		fmt.Println("total number of reactions: ", len(reactions))
	}

}
