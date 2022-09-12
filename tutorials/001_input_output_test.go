package tutorials_test

import (
	"fmt"
	"testing"

	"github.com/TimothyStiles/poly/io/genbank"
)

/******************************************************************************
Sep, 12, 2022

== Reading and Writing Biological Data ==

"There's 10,000 file formats for DNA and none of them are JSON."

- me

The bane of the bioinformatician's existence are the blasted file formats that
have been foisted upon us by our predecessors. Most are hard to understand and
even harder to parse programmatically.

My least favorite is the most common format for annotated sequence. The Genbank
format.

Made by some secret council of jerks in the deserts of New Mexico sometime in
1978 this Genbank format is used by the NIH for their massive genomic repository
that goes by the same name:

https://www.ncbi.nlm.nih.gov/genbank/

This format is gnarly and barely has any documentation or specification
associated with it. Figuring out how to parse it was a three week ordeal that
I've written about many times in many other places.

Point is. To start you're going to need to get data in and out of your computer
and Poly makes it easy. Follow along below and run each test through the
included debugger via gitpod if you really want to know what's going on.

TTFN,
Tim
******************************************************************************/

// if you're using VS-CODE you should see a DEBUG TEST button right below this
// comment. Please set break points and use it early and often.
func TestFileIOTutorial(t *testing.T) {

	// First we're going to read in a Genbank file for the well known plasmid
	// backbone puc19. Plasmids are super small rings of "Circular DNA" that are
	// between 1 and 10 kilobases in length.

	puc19, _ := genbank.Read("../data/puc19.gbk")

	// A plasmid backbone is an empty template designed so that we
	// can easily insert genes into it to be expressed in our organism of choice.
	// In this case the harmless e. coli strain DH5a which is not of the
	// dangerous Chipotle variety.

	// Think of plasmids like little power up tokens for cells.
	// In lieu of sexual reproduction bacteria use these plasmid
	// tokens as a method of maintaining genetic diversity and adapting to
	// adverse conditions. There are consensual and non-consensual methods
	// of plasmid exchange that bacteria participate in from formal conjugation
	// (pretty much bacteria sex) to simply eating the other bacteria to gain its
	// powers like Kirby or some weird anime character.

	// Since the 1970s people have been adapting these plasmids that are naturally
	// found in bacteria to work in other cells and cell types from all taxons
	// of life. There are human engineered plasmids for human cells, tobacco cells,
	// and even chinese hamster ovary (CHO) cells. As synthetic biologists we'll
	// often design plasmids to help us engineer organisms.

	// You can read more about plasmids here:
	// https://en.wikipedia.org/wiki/Plasmid

	// You can also check out this massive repository of plasmids:
	// https://www.addgene.org/

	// Now that we have our genbank file uploaded we can take a look at what's
	// inside. I HIGHLY SUGGEST opening this in the debugger provided so that you
	// can check out the entire structure of the data we're working with.

	// First lets check out the metadata of our puc19 backbone
	meta := puc19.Meta
	fmt.Println(meta)

	// The meta data locus is a genbank specific field
	metaName := meta.Locus.Name
	fmt.Println(metaName)

	sourceOrganism := meta.Source
	fmt.Println(sourceOrganism)

	// Next let's see what sort of features puc19 has
	// Features are the parts of the sequence that are annotated.

	for _, feature := range puc19.Features {
		fmt.Println(feature.Type)
	}

	// We'll go into more detail about features and DNA parts
	// in the next tutorial but for now know that we can also
	// get the sequence of each feature using the GetSequence method.

	randomFeature := puc19.Features[1]
	randomFeatureSequence, _ := randomFeature.GetSequence()
	fmt.Println(randomFeatureSequence)

	// we can also get the sequence of the entire plasmid
	sequence := puc19.Sequence
	fmt.Println(sequence)

	// As you've just seen Poly supports many different file formats for DNA most
	// of which share a general data structure containing these parts:

	// meta: A place where meta information is stored

	// features: a list of feature annotations that are associated with the
	// sequence including the sequence type (CDS, rRNA, etc), and where
	// the feature is located within the sequence.

	// sequence: the sequence string itself which is either composed of
	// nucleic acid notation or single letter amino acid notation.

	// Nucleic acid notation:
	// https://en.wikipedia.org/wiki/Nucleic_acid_notation

	// Amino acid notation:
	// https://en.wikipedia.org/wiki/Protein_primary_structure#Notation

}
