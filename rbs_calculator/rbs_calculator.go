package rbs_calculator

import (
	"fmt"
	"os"
	"regexp"
	"sort"

	"github.com/TimothyStiles/poly/rbs_calculator/model"
	rbs_model "github.com/TimothyStiles/poly/rbs_calculator/model/salis_lab_v2_1"
	"github.com/olekukonko/tablewriter"
)

/*****************************************************
June, 5, 2021
The Ribosome Binding Site Calculator

Translation is the process by which a messenger RNA (mRNA) is 'decoded' to a
protein by a ribosome.

When trying to synthesize proteins of our choice, we have to design a mRNA
sequence with a  5' untranslated region, the protein coding sequence and a
terminator sequence. The ribosome's interaction with the mRNA sequence
determines the amount of the protein that will be synthesized. Thus, we need
a model to understand these interactions and how they relate to the amount of
protein produced.

The ribosome binding site calculator returns a list of potential places (binding
sites) a ribosome might bind to a mRNA. Each binding site breaks down a mRNA
sequence into a 5' untranslated region and protein coding sequence, and includes
information about the predicted translation rate of that binding site. This is
useful because it allows one to examine where a ribosome is likely to bind to
the mRNA sequence by examining the relative differences between the translation
initiation rates of each binding site.

The other intended use of this calculator is to help one increase synthesis of
a desired protein by the user re-desinging their mRNA sequence to increase
translation initiation rate which should lead to an increase in the amount of
the protein produced in-vitro.

How the calculator works:
To calculate the translation initiation rate of a binding site, we need a model
of a ribosome's interactions with a mRNA strand and the relationship of the
interactions with the translation initiation rate of the coding sequence of the
mRNA strand.

The current model of this RBS calculator is based upon the research done by
The Salis Lab at The Penn State University. For more information on this model,
please have a look at the `properties.go` file in the
`rbs_calculator/model/salis_lab_v2_1` subpackage.

For information of how to develop a RBS Calculator of your own, please
have a look at the `model.go` file in the `rbs_calculator/model` subpackage.

*****************************************************/

// StartCodon specifies the start codon of the protein coding sequence of a
// mRNA strand.
type StartCodon string

const (
	AUG StartCodon = "AUG"
	GUG            = "GUG"
	UUG            = "UUG"
	CUG            = "CUG"
	AUC            = "AUC"
	AUA            = "AUA"
	AUU            = "AUU"
	CAU            = "CAU"
	GGA            = "GGA"
	UGC            = "UGC"
	CGC            = "CGC"
	UAG            = "UAG"
)

// Organism16SrRNAMap is a map of an organism to its 16S ribosomal RNA.
var Organism16SrRNAMap map[string]string = map[string]string{
	"Bacillus subtilis subsp. subtilis str. 168":                       "ACCUCCUUU",
	"Bacteroides thetaiotaomicron VPI-5482":                            "ACCUCCUUA",
	"Corynebacterium glutamicum B-2784":                                "ACCUCCUUU",
	"Escherichia coli BL21(DE3)":                                       "ACCUCCUUA",
	"Escherichia coli str. K-12 substr. DH10B":                         "ACCUCCUUA",
	"Escherichia coli str. K-12 substr. MG1655":                        "ACCUCCUUA",
	"Pseudomonas fluorescens A506":                                     "ACCUCCUUA",
	"Salmonella enterica subsp. enterica serovar Typhimurium str. LT2": "ACCUCCUUA",
}

// DefaultStartCodons are AUG, GUG and CUG
var DefaultStartCodons []StartCodon = []StartCodon{AUG, GUG, CUG}

// RibosomeBindingSites returns a list of ribosome binding sites for a given
// messenger RNA (mRNA) strand and 16s rRNA. Each site contains information
// about the five prime untranslated region, the protein coding sequence, the
// translation initiation rate as well as all the properties that were computed
// to figure out the translation initiation rate.
//
// The translation initiation rates of the binding sites can be compared with
// each other to help you figure out where a ribosome is likely to bind to
// your mRNA strand.
//
// The output can also be used to increase the amount of synthesis of your
// desired protein. Redesign your mRNA sequence to have higher translation
// initiation rates and this should correlate with an increase in the amount
// of protein synthesized by your mRNA strand in-vitro.
//
// Use the exported map `Organism16SrRNAMap` to find the 16s rRNA of the
// organism the mRNA strand is inserted into.
func RibosomeBindingSites(ribosomalRNA, mRNA string, temperatureInCelsius float64, startCodons []StartCodon) (ribosomeBindingSites []model.RibosomeBindingSite) {
	ribosomalRNA = model.CleanRNA(ribosomalRNA)
	mRNA = model.CleanRNA(mRNA)

	// for each start codon, find all occurrences of it in the mRNA sequence,
	// create a `model.RibosomeBindingSite` and compute the translation
	// initiation rate for each occurrence
	for _, startCondon := range startCodons {
		startCodonRegex := regexp.MustCompile(string(startCondon))
		matches := startCodonRegex.FindAllStringIndex(mRNA, -1)
		for _, match := range matches {
			startCodonIdx := match[0]
			fivePrimeUTR, codingSequence := mRNA[:startCodonIdx], mRNA[startCodonIdx:]

			_, ribosomeBindingSite := TranslationInitiationRate(fivePrimeUTR, codingSequence, ribosomalRNA, temperatureInCelsius)
			ribosomeBindingSites = append(ribosomeBindingSites, ribosomeBindingSite)
		}
	}

	// sort by length of 5' UTR in ascending order
	sort.Slice(ribosomeBindingSites, func(i, j int) bool {
		return ribosomeBindingSites[i].Properties["startPos"].(int) < ribosomeBindingSites[j].Properties["startPos"].(int)
	})
	return
}

// SortByTranslationInitiationRate sorts a list of binding sites by their
// translation initiation rates in descending order
func SortByTranslationInitiationRate(ribosomeBindingSites []model.RibosomeBindingSite) []model.RibosomeBindingSite {
	sort.Slice(ribosomeBindingSites, func(i, j int) bool {
		if ribosomeBindingSites[i].TranslationInitiationRate == ribosomeBindingSites[j].TranslationInitiationRate {
			return len(ribosomeBindingSites[i].FivePrimeUTR) < len(ribosomeBindingSites[j].FivePrimeUTR)
		}
		return ribosomeBindingSites[i].TranslationInitiationRate > ribosomeBindingSites[j].TranslationInitiationRate
	})
	return ribosomeBindingSites
}

// TranslationInitiationRate returns the the translation initiation rate of a
// ribsome binding site as well as the binding site with the properties
// computed to calculate the translation initiation rate
func TranslationInitiationRate(fivePrimeUTR, proteinCodingSequence, ribosomalRNA string, temperateureInCelsius float64) (translationInitiationRate float64, bindingSiteWithProperties model.RibosomeBindingSite) {
	ribosomeBindingSite := model.RibosomeBindingSite{
		FivePrimeUTR:          fivePrimeUTR,
		ProteinCodingSequence: proteinCodingSequence,
		Temperature:           temperateureInCelsius,
		RibosomalRNA:          ribosomalRNA,
		Properties:            make(map[string]interface{}),
		// we use inf as a sanity check to ensure `TranslationInitiationRate` has
		// been computed for the binding site
		TranslationInitiationRate: inf,
	}

	// compute the properties of the interactions between the mRNA strand
	// and ribosome required to calculate the translation initiation rate
	ribosomeBindingSite.ComputeProperties(rbs_model.PropertiesToCompute)

	// calculate the translation initiation rate
	translationInitiationRate = rbs_model.TranslationInitiationRate(ribosomeBindingSite)
	ribosomeBindingSite.TranslationInitiationRate = translationInitiationRate

	if ribosomeBindingSite.TranslationInitiationRate == inf {
		panic("failed to calculate the translation initiation rate for the mRNA sequence. ensure the `TranslationInitiationRate` field of the `model.RibosomeBindingSite` struct is set by the `rbs_model.TranslationInitiationRate` func.")
	}

	return translationInitiationRate, ribosomeBindingSite
}

// PrintBindingSites prints the important properties of a binding site.
// The optional argument `propertiesToPrint` specifies the computed properties
// of the binding site that will be printed. Note that the properties are
// case-sensitive and must match one of the properties available in the
// `model.RibosomeBindingSite.Properties` map after the properties of a binding
// site have been computed.
func PrintBindingSites(bindingSites []model.RibosomeBindingSite, includeSequences, includeStructures bool, propertiesToPrint ...string) {
	table := tablewriter.NewWriter(os.Stdout)
	table.SetRowLine(true)
	table.SetAutoFormatHeaders(false)

	// add the basic important properties of the binding site
	header := []string{
		"Start position",
		"Start codon",
		"TIR",
	}

	// add the additional properties
	header = append(header, propertiesToPrint...)

	if includeSequences {
		header = append(header, "5' Untranslated Region")
		header = append(header, "Protein Coding Sequence")
	}

	if includeStructures {
		header = append(header, "Initial state")

		// add the final state headers
		header = append(header, "Final state (pre ribosome)")
		header = append(header, "Final state (mRNA shine dalgarno binding site)")
		header = append(header, "Final state (spacing)")
		header = append(header, "Final state (ribosome footprint)")
		header = append(header, "Final state (post ribosome)")
		header = append(header, "Final state (16S rRNA shine dalgarno binding site)")
	}

	table.SetHeader(header)

	// add the relevant information for each binding site
	for _, bindingSite := range bindingSites {
		bindingSiteInfo := bindingSiteInformation(bindingSite, includeSequences, includeStructures, propertiesToPrint...)
		table.Append(bindingSiteInfo)
	}

	// finally, render the table
	table.Render()
}

// returns the information of a binding site
func bindingSiteInformation(bindingSite model.RibosomeBindingSite, includeSequences, includeStructures bool, propertiesToPrint ...string) (ret []string) {

	// add the basic important properties of the binding site
	data := []interface{}{
		bindingSite.Properties["startPos"],
		bindingSite.ProteinCodingSequence[:3],
		bindingSite.TranslationInitiationRate,
	}

	// add the additional properties
	for _, property := range propertiesToPrint {
		data = append(data, bindingSite.Properties[property])
	}

	if includeSequences {
		data = append(data, bindingSite.FivePrimeUTR)
		data = append(data, bindingSite.ProteinCodingSequence)
	}

	if includeStructures {
		initialState := bindingSite.Properties["usedMRNAStructure"]
		data = append(data, initialState)

		// add the final state structures
		sdBindingSiteMRNAStructure, sdBindingSiteRRNAStructure := bindingSite.Properties["sdBindingSiteMRNAStructure"], "&"+bindingSite.Properties["sdBindingSiteRRNAStructure"].(string)
		preRibosomeStructure, postRibosomeStructure := bindingSite.Properties["preRibosomeMRNAStructure"], bindingSite.Properties["postRibosomeMRNAStructure"]
		ribosomeFootprint := "............."
		nbSpacingNucleotides := bindingSite.Properties["lenFivePrimeUTR"].(int) - bindingSite.Properties["bindingSiteThreePrimeIdx"].(int)
		spacing := make([]rune, 0)
		for i := 0; i < nbSpacingNucleotides; i++ {
			spacing = append(spacing, '.')
		}
		finalState := []interface{}{preRibosomeStructure, sdBindingSiteMRNAStructure, string(spacing), ribosomeFootprint, postRibosomeStructure, sdBindingSiteRRNAStructure}
		data = append(data, finalState...)
	}

	ret = stringSlice(data)
	return
}

func stringSlice(interfaceSlice []interface{}) (ret []string) {
	ret = make([]string, len(interfaceSlice))
	for i, v := range interfaceSlice {
		ret[i] = fmt.Sprint(v)
	}
	return
}

var inf float64 = 1000000000.0
