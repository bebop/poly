package main

import (
	"bytes"
	"log"
	"regexp"
	"strconv"
	"strings"
)

var genbankDivisions = []string{
	"PRI", //primate sequences
	"ROD", //rodent sequences
	"MAM", //other mamallian sequences
	"VRT", //other vertebrate sequences
	"INV", //invertebrate sequences
	"PLN", //plant, fungal, and algal sequences
	"BCT", //bacterial sequences
	"VRL", //viral sequences
	"PHG", //bacteriophage sequences
	"SYN", //synthetic sequences
	"UNA", //unannotated sequences
	"EST", //EST sequences (expressed sequence tags)
	"PAT", //patent sequences
	"STS", //STS sequences (sequence tagged sites)
	"GSS", //GSS sequences (genome survey sequences)
	"HTG", //HTG sequences (high-throughput genomic sequences)
	"HTC", //unfinished high-throughput cDNA sequencing
	"ENV", //environmental sampling sequences
}

var genbankTopLevelFeatures = []string{
	"LOCUS",
	"DEFINITION",
	"ACCESSION",
	"VERSION",
	"KEYWORDS",
	"SOURCE",
	"REFERENCE",
	"FEATURES",
	"ORIGIN",
}

var genbankSubLevelFeatures = []string{
	"ORGANISM",
	"AUTHORS",
	"TITLE",
	"JOURNAL",
	"PUBMED",
	"REMARK",
}

func parseLocus(locusString string) Locus {
	locus := Locus{}
	locusSplit := strings.Split(strings.TrimSpace(locusString), " ")
	var filteredLocusSplit []string
	for i := range locusSplit {
		if locusSplit[i] != "" {
			filteredLocusSplit = append(filteredLocusSplit, locusSplit[i])
		}
	}
	locus.Name = filteredLocusSplit[1]
	locus.SequenceLength = strings.Join([]string{filteredLocusSplit[2], filteredLocusSplit[3]}, " ")
	locus.MoleculeType = filteredLocusSplit[4]
	if filteredLocusSplit[5] == "circular" || filteredLocusSplit[5] == "linear" {
		if filteredLocusSplit[5] == "circular" {
			locus.Circular = true
		} else {
			locus.Circular = false
		}
		locus.GenBankDivision = filteredLocusSplit[6]
		locus.ModDate = filteredLocusSplit[7]
	} else {
		locus.Circular = false
		locus.GenBankDivision = filteredLocusSplit[5]
		locus.ModDate = filteredLocusSplit[6]
	}
	return locus
}

func joinSubLines(splitLine, subLines []string) string {
	base := strings.TrimSpace(strings.Join(splitLine[1:], " "))
	for _, subLine := range subLines {
		if string(subLine[0]) == " " {
			base = strings.TrimSpace(strings.TrimSpace(base) + " " + strings.TrimSpace(subLine))
		} else {
			break
		}
	}
	return base
}

func getSourceOrganism(splitLine, subLines []string) (string, string) {
	source := strings.TrimSpace(strings.Join(splitLine[1:], " "))
	var organism string
	for numSubLine, subLine := range subLines {
		headString := strings.Split(strings.TrimSpace(subLine), " ")[0]
		if string(subLine[0]) == " " && headString != "ORGANISM" {
			source = strings.TrimSpace(strings.TrimSpace(source) + " " + strings.TrimSpace(subLine))
		} else {
			organismSubLines := subLines[numSubLine+1:]
			organismSplitLine := strings.Split(strings.TrimSpace(subLine), " ")
			organism = joinSubLines(organismSplitLine, organismSubLines)
			break
		}
	}
	return source, organism
}

func joinReferenceSubLines(splitLine, subLines []string) string {
	base := strings.TrimSpace(strings.Join(splitLine[1:], " "))

	for _, subLine := range subLines {
		featureSplitLine := strings.Split(strings.TrimSpace(subLine), " ")
		headString := featureSplitLine[0]
		breakFlag := false
		//breaks if there's a sub level feature.
		for _, subLevelFeature := range genbankSubLevelFeatures {
			if headString == subLevelFeature {
				breakFlag = true
			}
		}
		//breaks if there's a top level feature.
		for _, topLevelFeature := range genbankTopLevelFeatures {
			if headString == topLevelFeature {
				breakFlag = true
			}
		}
		//if there's no need to break append this line to the last and trim.
		if breakFlag != true {
			base = strings.TrimSpace(strings.TrimSpace(base) + " " + strings.TrimSpace(subLine))
		} else {
			break
		}
	}
	return base

}

func getReference(splitLine, subLines []string) Reference {
	base := strings.TrimSpace(strings.Join(splitLine[1:], " "))
	reference := Reference{}
	reference.Index = strings.Split(base, " ")[0]
	start, _ := strconv.Atoi(strings.TrimSpace(strings.Split(base, " ")[3]))
	reference.Start = start
	end := strings.TrimSpace(strings.Split(base, " ")[5])
	end = end[0 : len(end)-1]
	reference.End, _ = strconv.Atoi(end)

	for numSubLine, subLine := range subLines {
		featureSubLines := subLines[numSubLine+1:]
		featureSplitLine := strings.Split(strings.TrimSpace(subLine), " ")
		headString := featureSplitLine[0]
		// fmt.Println(headString)

		breakFlag := false
		for _, topLevelFeature := range genbankTopLevelFeatures {
			if headString == topLevelFeature {
				breakFlag = true
			}
		}
		if breakFlag {
			break
		}
		switch headString {
		case "AUTHORS":
			reference.Authors = joinReferenceSubLines(featureSplitLine, featureSubLines)
		case "TITLE":
			reference.Title = joinReferenceSubLines(featureSplitLine, featureSubLines)
		case "JOURNAL":
			reference.Journal = joinReferenceSubLines(featureSplitLine, featureSubLines)
		case "PUBMED":
			reference.PubMed = joinReferenceSubLines(featureSplitLine, featureSubLines)
		case "REMARK":
			reference.Remark = joinReferenceSubLines(featureSplitLine, featureSubLines)
		default:
			break
		}

	}
	return reference
}

// func getFeatures(splitLine, subLines []string) []Feature {

// }

func getSequence(subLines []string) Sequence {
	sequence := Sequence{}
	var sequenceBuffer bytes.Buffer
	reg, err := regexp.Compile("[^a-zA-Z]+")
	if err != nil {
		log.Fatal(err)
	}
	for _, subLine := range subLines {
		sequenceBuffer.WriteString(reg.ReplaceAllString(subLine, ""))
	}
	sequence.Sequence = sequenceBuffer.String()

	return sequence
}
