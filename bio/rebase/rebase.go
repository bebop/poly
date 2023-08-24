/*
Package rebase contains a rebase parser for rebase data dump #31.

In order to effectively simulate cloning reactions, we need to know how each
restriction enzyme in the reaction functions. This data can be derived, in
bulk, from the REBASE database.

REBASE is an amazing resource run by New England Biolabs listing essentially
every known restriction enzyme. In particular, this parser parses the REBASE
data dump format #31, which is what Bioperl uses.

https://bioperl.org/howtos/Restriction_Enzyme_Analysis_HOWTO.html
http://rebase.neb.com/rebase/rebase.f31.html

The actual data dump itself is linked here and updated once a month:
http://rebase.neb.com/rebase/link_withrefm

The header of this file gives a wonderful explanation of its structure. Here is the
header with the commercial suppliers format and an example enzyme.

```
REBASE version 104                                              withrefm.104

	=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=
	REBASE, The Restriction Enzyme Database   http://rebase.neb.com
	Copyright (c)  Dr. Richard J. Roberts, 2021.   All rights reserved.
	=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=

# Rich Roberts                                                    Mar 31 2021

<ENZYME NAME>   Restriction enzyme name.
<ISOSCHIZOMERS> Other enzymes with this specificity.
<RECOGNITION SEQUENCE>

	These are written from 5' to 3', only one strand being given.
	If the point of cleavage has been determined, the precise site
	is marked with ^.  For enzymes such as HgaI, MboII etc., which
	cleave away from their recognition sequence the cleavage sites
	are indicated in parentheses.

	For example HgaI GACGC (5/10) indicates cleavage as follows:
	                5' GACGCNNNNN^      3'
	                3' CTGCGNNNNNNNNNN^ 5'

	In all cases the recognition sequences are oriented so that
	the cleavage sites lie on their 3' side.

	REBASE Recognition sequences representations use the standard
	abbreviations (Eur. J. Biochem. 150: 1-5, 1985) to represent
	ambiguity.
	                R = G or A
	                Y = C or T
	                M = A or C
	                K = G or T
	                S = G or C
	                W = A or T
	                B = not A (C or G or T)
	                D = not C (A or G or T)
	                H = not G (A or C or T)
	                V = not T (A or C or G)
	                N = A or C or G or T



	ENZYMES WITH UNUSUAL CLEAVAGE PROPERTIES:

	Enzymes that cut on both sides of their recognition sequences,
	such as BcgI, Bsp24I, CjeI and CjePI, have 4 cleavage sites
	each instead of 2.

	Bsp24I
	          5'      ^NNNNNNNNGACNNNNNNTGGNNNNNNNNNNNN^   3'
	          3' ^NNNNNNNNNNNNNCTGNNNNNNACCNNNNNNN^        5'


	This will be described in some REBASE reports as:

	             Bsp24I (8/13)GACNNNNNNTGG(12/7)

<METHYLATION SITE>

	The site of methylation by the cognate methylase when known
	is indicated X(Y) or X,X2(Y,Y2), where X is the base within
	the recognition sequence that is modified.  A negative number
	indicates the complementary strand, numbered from the 5' base
	of that strand, and Y is the specific type of methylation
	involved:
	               (6) = N6-methyladenosine
	               (5) = 5-methylcytosine
	               (4) = N4-methylcytosine

	If the methylation information is different for the 3' strand,
	X2 and Y2 are given as well.

<MICROORGANISM> Organism from which this enzyme had been isolated.
<SOURCE>        Either an individual or a National Culture Collection.
<COMMERCIAL AVAILABILITY>

	Each commercial source of restriction enzymes and/or methylases
	listed in REBASE is assigned a single character abbreviation
	code.  For example:

	K        Takara (1/98)
	M        Boehringer Mannheim (10/97)
	N        New England Biolabs (4/98)

	The date in parentheses indicates the most recent update of
	that organization's listings in REBASE.

<REFERENCES>only the primary references for the isolation and/or purification
of the restriction enzyme or methylase, the determination of the recognition
sequence and cleavage site or the methylation specificity are given.

REBASE codes for commercial sources of enzymes

	B        Life Technologies (3/21)
	C        Minotech Biotechnology (3/21)
	E        Agilent Technologies (8/20)
	I        SibEnzyme Ltd. (3/21)
	J        Nippon Gene Co., Ltd. (3/21)
	K        Takara Bio Inc. (6/18)
	M        Roche Applied Science (4/18)
	N        New England Biolabs (3/21)
	O        Toyobo Biochemicals (8/14)
	Q        Molecular Biology Resources - CHIMERx (3/21)
	R        Promega Corporation (11/20)
	S        Sigma Chemical Corporation (3/21)
	V        Vivantis Technologies (1/18)
	X        EURx Ltd. (1/21)
	Y        SinaClon BioScience Co. (1/18)

<1>AaaI
<2>XmaIII,BseX3I,BsoDI,BstZI,EagI,EclXI,Eco52I,SenPT16I,TauII,Tsp504I
<3>C^GGCCG
<4>
<5>Acetobacter aceti ss aceti
<6>M. Fukaya
<7>
<8>Tagami, H., Tayama, K., Tohyama, T., Fukaya, M., Okumura, H., Kawamura, Y., Horinouchi, S., Beppu, T., (1988) FEMS Microbiol. Lett., vol. 56, pp. 161-166.

```
*/
package rebase

import (
	"encoding/json"
	"io"
	"os"
	"strings"
)

var (
	readAllFn  = io.ReadAll
	parseFn    = Parse
	marshallFn = json.Marshal
)

// Enzyme represents a single enzyme within the Rebase database
type Enzyme struct {
	Name                   string   `json:"name"`
	Isoschizomers          []string `json:"isoschizomers"`
	RecognitionSequence    string   `json:"recognitionSequence"`
	MethylationSite        string   `json:"methylationSite"`
	MicroOrganism          string   `json:"microorganism"`
	Source                 string   `json:"source"`
	CommercialAvailability []string `json:"commercialAvailability"`
	References             string   `json:"references"`
}

// Parse parses the Rebase database into a map of enzymes
func Parse(file io.Reader) (map[string]Enzyme, error) {
	fileBytes, err := readAllFn(file)
	if err != nil {
		return make(map[string]Enzyme), err
	}
	// Setup some variables
	var enzyme Enzyme
	enzymeMap := make(map[string]Enzyme)
	commercialSupplierMap := make(map[rune]string)

	// Get rebase as a large string
	rebase := string(fileBytes)

	// Split those strings into individual lines for parsing
	lines := strings.Split(rebase, "\n")

	commercialParsingLine := 0
	startCommercialParsing := false
	startReferenceParsing := false
	for _, line := range lines {
		// Parse commercial sources map
		if line == "REBASE codes for commercial sources of enzymes" {
			startCommercialParsing = true
		}
		// If startCommercialParsing is true, start building the commercial supplier map
		if startCommercialParsing {
			// if we start enzyme parsing, break the commercial supplier parsing
			if strings.Contains(line, "<1>") {
				commercialParsingLine = 0
				startCommercialParsing = false
			}

			// Skip two lines
			commercialParsingLine++
			if (commercialParsingLine > 3) && (len(strings.TrimLeft(line, "\t")) > 0) {
				// Trim indentation
				trimmedString := strings.TrimLeft(line, "\t")

				// The first letter of the trimmedString is the single letter code
				singleLetterCommercialCode := rune(trimmedString[0])

				// There are 8 spaces until the commercial companies's name. We are keeping the dates
				// attached, since it is additional information that could be useful for users down
				// the line
				commercialName := trimmedString[9:]

				// Add both to commercialSupplierMap
				commercialSupplierMap[singleLetterCommercialCode] = commercialName
			}
		}

		// If we are parsing references, continue appending to the current enzyme's references
		if startReferenceParsing && line != "" {
			// Break reference parsing if we encounter a new enzyme
			if strings.Contains(line, "<1>") {
				enzymeMap[enzyme.Name] = enzyme
				enzyme = Enzyme{}
				startReferenceParsing = false
			}

			enzyme.References += "\n" + line
		}

		// Normal enzyme parsing
		switch {
		case strings.Contains(line, "<1>"):
			enzyme.Name = line[3:]
		case strings.Contains(line, "<2>"):
			enzyme.Isoschizomers = strings.Split(line[3:], ",")
		case strings.Contains(line, "<3>"):
			enzyme.RecognitionSequence = line[3:]
		case strings.Contains(line, "<4>"):
			enzyme.MethylationSite = line[3:]
		case strings.Contains(line, "<5>"):
			enzyme.MicroOrganism = line[3:]
		case strings.Contains(line, "<6>"):
			enzyme.Source = line[3:]
		case strings.Contains(line, "<7>"):
			// We need to get a list of specific commercial suppliers from the commercialSupplierMap we previously made
			var commercialSuppliers []string
			for _, commercialLetter := range line[3:] {
				commercialSuppliers = append(commercialSuppliers, commercialSupplierMap[commercialLetter])
			}
			enzyme.CommercialAvailability = commercialSuppliers
		case strings.Contains(line, "<8>"):
			enzyme.References = line[3:]
			startReferenceParsing = true
		}
	}
	return enzymeMap, err
}

// Read returns an enzymeMap from a Rebase data dump
func Read(path string) (map[string]Enzyme, error) {
	file, err := os.Open(path)
	if err != nil {
		return map[string]Enzyme{}, err
	}
	enzymeMap, err := parseFn(file)
	if err != nil {
		return map[string]Enzyme{}, err
	}
	return enzymeMap, nil
}

// Export returns a json file of the Rebase database
func Export(enzymeMap map[string]Enzyme) ([]byte, error) {
	jsonRebase, err := marshallFn(enzymeMap)
	if err != nil {
		return []byte{}, err
	}
	return jsonRebase, nil
}
