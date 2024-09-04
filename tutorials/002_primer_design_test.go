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

Now that you've learned what plasmids are you should probably know what primers
are as well.

"Primers are short sequences of DNA that can be used to amplify DNA sequences
and they are the workhorse of modern molecular biology.

Essentially primers are short pieces of single stranded DNA that can
bind to a target sequence of single stranded DNA. These primers serve as a
marker for polymerases (the enzyme Poly is named after!) to bind and start adding
free floating nucleotides (ACTGs) to a single strand piece of DNA to form a
double stranded piece of DNA.

This is a crucial step in the process of PCR (polymerase chain reaction).
https://en.wikipedia.org/wiki/Polymerase_chain_reaction

Here's also a video animation from Cold Spring Harbor's DNA Learning center explaining the process.
https://youtu.be/2KoLnIwoZKU?si=wqKs1NU5ZhU5O7Ui

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

We could also use these primers for RNA interference experiments to suppress
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
			forward, reverse := pcr.DesignPrimers(reaction.Sequence, 56.0) // <- 56.0 is our melting temp. The temperature at which we want our primers to bind to double stranded DNA. Again. don't hardcode values like this in real life. Put it in a constant or something.
			reaction.ForwardPrimer = forward
			reaction.ReversePrimer = reverse

			// append our reaction to a our reactions slice (slice is essentially Go's version of a list, or vector)
			reactions = append(reactions, reaction)
		}
	}

	fmt.Println("Total reactions:", len(reactions))
	// We've now just generated ~5000 primers.
	// Notice how the only numerical parameter we give was "melting temp" this is the temp at which they'll anneal to denatured DNA.
	// As mentioned in this Cold Spring Harbor video. PCR reactions are conducted in ~30 cycles of heating and cooling over 3 temperature stages.

	// Stage 1: We raise the temperature of our reaction to 95C to denature (split) our DNA such that our primers can bind to it.
	// Stage 2: We lower the temperature of our reaction to 55C so that our primers can bind to complementary regions of newly accessible single stranded DNA.
	// Stage 3: We raise the temperature of our reaction to ~72C (or whatever temperature is best for the polymerase we're using) to activate our polymerase
	// Stage 3 (cont): to bind to our primers and begin constructing a brand new second strand to our denatured DNA. Then we go back to Stage 1.

	// For each cycle of the above stage we end up doubling the number of copies of the gene we want so after we have n^30 copies of our desired region.

	// What we've done is design a gigantic set of primers that share the same melting temp. This makes it possible to run all of these reactions
	// concurrently within a single PCR run but we've ignored a lot of other design considerations.

	// This primers aren't particularly well designed. All pcr.DesignPrimers has done is figure out how long each primer should be so that they all bind
	// at a specific temperature while assuming that they should bind at the very beginning and end of each given sequence. This is an extremely common
	// use case but there are several caveats.

	// 1. The designed primers could be dimers (The primers could bind to each other and not to their target sequence)
	// 2. One primer in a pair could be a "hairpin" that binds to itself (which is something the fold package can help detect)
	// 3. Your primers may actually need a specific configuration to bind to the intended target (something the fold package may be able to help with)

	// Depending on your situation you may need to get creative. Lots of scientists design their primers to bind up or downstream of their gene of
	// interest to avoid primer dimers or hairpins. Some scientists don't care if they copy the whole gene but just want to copy enough of the
	// gene to verify that it's there. There's probably a million different way to design primers but I'd guess that poly itself covers about 95%
	// of what most scientists would need on a daily basis.
}
