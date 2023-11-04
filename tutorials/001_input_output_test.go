package tutorials_test

import (
	"encoding/json"
	"fmt"
	"os"
	"path/filepath"
	"testing"

	"github.com/TimothyStiles/poly/bio"
	"github.com/TimothyStiles/poly/bio/genbank"
	"github.com/google/go-cmp/cmp"
	"github.com/google/go-cmp/cmp/cmpopts"
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

// Ignore ParentSequence as that's a pointer which can't be serialized.
func CmpOptions() []cmp.Option {
	return []cmp.Option{
		cmpopts.IgnoreFields(genbank.Feature{}, "ParentSequence"),
	}
}

// if you're using VS-CODE you should see a DEBUG TEST button right below this
// comment. Please set break points and use it early and often.
func TestFileIOTutorial(t *testing.T) {
	// First we're going to read in a Genbank file for the well known plasmid
	// backbone puc19. Plasmids are super small rings of "Circular DNA" that are
	// between 1 and 10 kilobases in length.

	file, _ := os.Open("../bio/genbank/data/puc19.gbk")
	defer file.Close()
	parser, _ := bio.NewGenbankParser(file)
	puc19, _ := parser.Next()

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

	// The metadata Locus is a genbank specific field
	metaName := meta.Locus.Name     // this is what the genbank file lists the name as
	expectedMetaName := "puc19.gbk" // this is what we expect the name of the record to be

	if metaName != expectedMetaName {
		t.Errorf("Expected puc19 to be named %s but got %s", expectedMetaName, metaName)
	}

	sourceOrganism := meta.Source                       // this is the source organism from the genbank file
	expectedSourceOrganism := "synthetic DNA construct" // this is what we expect the source organism to be

	if sourceOrganism != expectedSourceOrganism {
		t.Errorf("Expected puc19 to be %s but got %s", expectedSourceOrganism, sourceOrganism)
	}

	// Next let's see what sort of features puc19 has
	// Features are the parts of the sequence that are annotated.

	for _, feature := range puc19.Features {
		fmt.Println(feature.Type)
	}

	// We'll go into more detail about features and DNA parts
	// in the next tutorial but for now know that we can also
	// get the sequence of each feature using the GetSequence method.

	feature := puc19.Features[1]
	featureSequence, _ := feature.GetSequence()       // this is the sequence of the feature
	expectedFeatureSequence := "gggaaacgcctggtatcttt" // this is what we expect the sequence of the feature to be
	if featureSequence != expectedFeatureSequence {
		t.Errorf("Expected feature sequence to be %s but got %s", expectedFeatureSequence, featureSequence)
	}

	// we can also get the sequence of the entire plasmid
	plasmidSequence := puc19.Sequence // this is the sequence of the plasmid

	// this is what we expect the sequence of the plasmid to be
	expectedPlasmidSequence := "gagatacctacagcgtgagctatgagaaagcgccacgcttcccgaagggagaaaggcggacaggtatccggtaagcggcagggtcggaacaggagagcgcacgagggagcttccagggggaaacgcctggtatctttatagtcctgtcgggtttcgccacctctgacttgagcgtcgatttttgtgatgctcgtcaggggggcggagcctatggaaaaacgccagcaacgcggcctttttacggttcctggccttttgctggccttttgctcacatgttctttcctgcgttatcccctgattctgtggataaccgtattaccgcctttgagtgagctgataccgctcgccgcagccgaacgaccgagcgcagcgagtcagtgagcgaggaagcggaagagcgcccaatacgcaaaccgcctctccccgcgcgttggccgattcattaatgcagctggcacgacaggtttcccgactggaaagcgggcagtgagcgcaacgcaattaatgtgagttagctcactcattaggcaccccaggctttacactttatgcttccggctcgtatgttgtgtggaattgtgagcggataacaatttcacacaggaaacagctatgaccatgattacgccaagcttgcatgcctgcaggtcgactctagaggatccccgggtaccgagctcgaattcactggccgtcgttttacaacgtcgtgactgggaaaaccctggcgttacccaacttaatcgccttgcagcacatccccctttcgccagctggcgtaatagcgaagaggcccgcaccgatcgcccttcccaacagttgcgcagcctgaatggcgaatggcgcctgatgcggtattttctccttacgcatctgtgcggtatttcacaccgcatatggtgcactctcagtacaatctgctctgatgccgcatagttaagccagccccgacacccgccaacacccgctgacgcgccctgacgggcttgtctgctcccggcatccgcttacagacaagctgtgaccgtctccgggagctgcatgtgtcagaggttttcaccgtcatcaccgaaacgcgcgagacgaaagggcctcgtgatacgcctatttttataggttaatgtcatgataataatggtttcttagacgtcaggtggcacttttcggggaaatgtgcgcggaacccctatttgtttatttttctaaatacattcaaatatgtatccgctcatgagacaataaccctgataaatgcttcaataatattgaaaaaggaagagtatgagtattcaacatttccgtgtcgcccttattcccttttttgcggcattttgccttcctgtttttgctcacccagaaacgctggtgaaagtaaaagatgctgaagatcagttgggtgcacgagtgggttacatcgaactggatctcaacagcggtaagatccttgagagttttcgccccgaagaacgttttccaatgatgagcacttttaaagttctgctatgtggcgcggtattatcccgtattgacgccgggcaagagcaactcggtcgccgcatacactattctcagaatgacttggttgagtactcaccagtcacagaaaagcatcttacggatggcatgacagtaagagaattatgcagtgctgccataaccatgagtgataacactgcggccaacttacttctgacaacgatcggaggaccgaaggagctaaccgcttttttgcacaacatgggggatcatgtaactcgccttgatcgttgggaaccggagctgaatgaagccataccaaacgacgagcgtgacaccacgatgcctgtagcaatggcaacaacgttgcgcaaactattaactggcgaactacttactctagcttcccggcaacaattaatagactggatggaggcggataaagttgcaggaccacttctgcgctcggcccttccggctggctggtttattgctgataaatctggagccggtgagcgtgggtctcgcggtatcattgcagcactggggccagatggtaagccctcccgtatcgtagttatctacacgacggggagtcaggcaactatggatgaacgaaatagacagatcgctgagataggtgcctcactgattaagcattggtaactgtcagaccaagtttactcatatatactttagattgatttaaaacttcatttttaatttaaaaggatctaggtgaagatcctttttgataatctcatgaccaaaatcccttaacgtgagttttcgttccactgagcgtcagaccccgtagaaaagatcaaaggatcttcttgagatcctttttttctgcgcgtaatctgctgcttgcaaacaaaaaaaccaccgctaccagcggtggtttgtttgccggatcaagagctaccaactctttttccgaaggtaactggcttcagcagagcgcagataccaaatactgttcttctagtgtagccgtagttaggccaccacttcaagaactctgtagcaccgcctacatacctcgctctgctaatcctgttaccagtggctgctgccagtggcgataagtcgtgtcttaccgggttggactcaagacgatagttaccggataaggcgcagcggtcgggctgaacggggggttcgtgcacacagcccagcttggagcgaacgacctacaccgaact"

	if plasmidSequence != expectedPlasmidSequence {
		t.Errorf("Expected plasmid sequence to be %s but got %s", expectedPlasmidSequence, plasmidSequence)
	}

	// Now that we've explored the plasmid, let's change it a tiny bit and write it out to both
	// a GenBank file and a JSON file.

	// First, let's change the name of the plasmid to "pUC19_modified"
	puc19.Meta.Locus.Name = "pUC19_modified"

	// adding ourselves as the modification author
	var reference genbank.Reference
	reference.Authors = "Timothy Stiles"
	reference.Title = "Modified pUC19"
	reference.Journal = "Poly"
	reference.PubMed = "123456789"

	puc19.Meta.References = append(puc19.Meta.References, reference)

	// create a tempdir to write the files to. The tempdir will be deleted when the test is done.
	tmpDataDir, err := os.MkdirTemp("", "data-*")
	if err != nil {
		t.Error(err)
	}
	defer os.RemoveAll(tmpDataDir)

	// write the modified plasmid to a GenBank file
	puc19Path := filepath.Join(tmpDataDir, "pUC19_modified.gb")
	writeFile, _ := os.OpenFile(puc19Path, os.O_WRONLY|os.O_CREATE|os.O_TRUNC, 0644)
	defer writeFile.Close()
	_, err = puc19.WriteTo(writeFile)
	if err != nil {
		t.Error(err)
	}

	// read the plasmid back in from the GenBank file and make sure it's the same as the original
	fileCopy, _ := os.Open(puc19Path)
	defer fileCopy.Close()
	parser2, _ := bio.NewGenbankParser(fileCopy)
	puc19Copy, err := parser2.Next()
	if err != nil {
		t.Error(err)
	}

	// compare our read-in plasmid to the the one we wrote out.
	if diff := cmp.Diff(puc19, puc19Copy, CmpOptions()...); diff != "" {
		t.Errorf("Parsing the output of Build() does not produce the same output as parsing the original file, \"%s\", read with Read(). Got this diff:\n%s", filepath.Base(puc19Path), diff)
	}

	// write the modified plasmid to a JSON file
	puc19JSONPath := filepath.Join(tmpDataDir, "pUC19_modified.json")
	marshaledPuc19, err := json.MarshalIndent(*puc19, "", " ")
	if err != nil {
		t.Error(err)
	}

	// write the JSON file
	_ = os.WriteFile(puc19JSONPath, marshaledPuc19, 0644)

	// read the plasmid back in from the JSON file and make sure it's the same as the original
	jsonContent, err := os.ReadFile(puc19JSONPath)
	if err != nil {
		t.Error(err)
	}
	var unmarshaledPuc19 genbank.Genbank
	if err := json.Unmarshal(jsonContent, &unmarshaledPuc19); err != nil {
		t.Error(err)
	}
	if diff := cmp.Diff(puc19, &unmarshaledPuc19, CmpOptions()...); diff != "" {
		t.Errorf("Parsing the JSON does not produce the same output as parsing the original file, \"%s\", read with Read(). Got this diff:\n%s", filepath.Base(puc19Path), diff)
	}

	// Glossary of terms:

	// meta: A place where metadata information is stored

	// features: a list of feature annotations that are associated with the
	// sequence including the sequence type (CDS, rRNA, etc), and where
	// the feature is located within the sequence.

	// sequence: the sequence string itself which is either composed of
	// nucleic acid notation or single letter amino acid notation.

	// more references

	// Nucleic acid notation:
	// https://en.wikipedia.org/wiki/Nucleic_acid_notation

	// Amino acid notation:
	// https://en.wikipedia.org/wiki/Protein_primary_structure#Notation
}
