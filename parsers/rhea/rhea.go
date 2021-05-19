package rhea

import (
	"bufio"
	"compress/gzip"
	"encoding/json"
	"encoding/xml"
	"errors"
	"io"
	"io/ioutil"
	"os"
	"strconv"
	"strings"
)

/******************************************************************************

Higher level structs

These structs are what you would put into a database or directly use. In order to
create a tree or insert into a normalized database, you would insert in the following
order:

	- Compound
	- ReactionSide
	- ReactionParticipant
	- Reaction
	- ReactionSide <-> Reaction

The entire structure of Rhea can simply be thought of as:

	1 There are Reactions. Those Reactions can have substrates and products, or substratesOrProducts
	  in the case that the reaction is bidirectional.
	2 There are ReactionSides. ReactionSides can be thought of as a many-to-many table between Reactions
	  and ReactionParticipants. It only serves as an inbetween, saying "this Reaction has these
	  ReactionParticipants on this side".
	3 There are ReactionParticipants. ReactionParticipants link ReactionSides with ReactiveParts and include
	  useful information like the number of ReactiveParts/Compounds (or chemicals) needed to do a certain Reaction.
	4 There are Compounds. These are chemicals. Sometimes they're complicated chemicals (like proteins) that
	  have ReactionSides, or portions of the larger chemical that are actually active in the reaction, though often
	  they're just small molecules.

With the "->" representing a "Reaction", and the right and left []ReactionParticipant each representing a "ReactionSide",
the entire Rhea database is basically a series of:

[]ReactionParticipant -> []ReactionParticipant

******************************************************************************/

// Rhea is a struct of the entire Rhea Database in a simplified higher-level way.
type Rhea struct {
	ReactionParticipants []ReactionParticipant `json:"reactionParticipants"`
	Compounds            []Compound            `json:"compounds"`
	Reactions            []Reaction            `json:"reactions"`
}

// Compound is a struct of a Rhea compound. These are chemicals - sometimes they are
// complicated chemicals like proteins or polymers, but often they are simple molecules.
// When Compounds are complicated molecules, like proteins, they can have reactiveParts
// that actually do the reactions. Compound, as defined here, is a pair containing a
// macro-molecule and a reactivePart. When building into a database, you will have to
// split this pairing for normalization.
type Compound struct {
	ID                  int    `json:"id" db:"id"`
	Accession           string `json:"accession" db:"accession"`
	Position            string `json:"position" db: "position"`
	Name                string `json:"name" db:"name"`
	HTMLName            string `json:"htmlName" db:"htmlname"`
	Formula             string `json:"formula" db:"formula"`
	Charge              string `json:"charge" db:"charge"`
	ChEBI               string `json:"chebi" db:"chebi"`
	SubclassOfChEBI     string `json:"subclassOfChEBI"`
	PolymerizationIndex string `json:"polymerizationIndex" db:"polymerizationindex"`
	CompoundID          int    `json:"id" db:"compoundid"`
	CompoundAccession   string `json:"accession" db:"compoundaccession"`
	CompoundName        string `json:"name" db:"compoundname"`
	CompoundHTMLName    string `json:"htmlName" db:"compoundhtmlname"`
	CompoundType        string `json:"compoundType" db:"compoundtype"`
}

// ReactionParticipant represents a Rhea ReactionParticipant. ReactionParticipants represent Compounds in Reactions. They
// can contain many of the same Compounds (including polymerized configurations) and are associated with a ReactionSide.
type ReactionParticipant struct {
	ReactionSide string `json:"reactionside" db:"reactionside"`
	Contains     int    `json:"contains" db:"contains"`
	ContainsN    bool   `json:"containsn" db:"containsn"`
	Minus        bool   `json:"minus" db:"minus"` // Only set to true if ContainsN == true to handle Nminus1
	Plus         bool   `json:"plus" db:"plus"`   // Only set to true if ContainsN == true to handle Nplus1
	Accession    string `json:"reactionParticipant" db:"ReactionParticipant"`
	Compound     string `json:"compound" db:"compound"`
}

// Reaction represents a Rhea reaction. Substrates, Products, and SubstrateOrProducts are all ReactionSide accession
// numbers, which can be linked to the ReactionParticipant's ReactionSide accession
type Reaction struct {
	ID                   int      `json:"id" db:"id"`
	Directional          bool     `json:"directional" db:"directional"`
	Accession            string   `json:"accession" db:"accession"`
	Status               string   `json:"status" db:"status"`
	Comment              string   `json:"comment" db:"comment"`
	Equation             string   `json:"equation" db:"equation"`
	HTMLEquation         string   `json:"htmlequation" db:"htmlequation"`
	IsChemicallyBalanced bool     `json:"ischemicallybalanced" db:"ischemicallybalanced"`
	IsTransport          bool     `json:"istransport" db:"istransport"`
	Ec                   string   `json:"ec" db:"ec"`
	Location             string   `json:"location" db:"location"`
	Citations            []string `json:"citations"`
	Substrates           []string `json:"substrates"`
	Products             []string `json:"products"`
	SubstrateOrProducts  []string `json:"substrateOrProducts"`
}

/******************************************************************************

Parse functions

These functions take in the rhea.rdf.gz dump file and return a Rhea struct,
which contains all of the higher level structs

******************************************************************************/

// NewReaction returns a Reaction.
func NewReaction(description Description, subclass Subclass) Reaction {
	return Reaction{
		ID:                   description.ID,
		Directional:          subclass.Resource == "http://rdf.rhea-db.org/DirectionalReaction",
		Accession:            description.Accession,
		Status:               description.Status.Resource,
		Comment:              description.Comment,
		Equation:             description.Equation,
		HTMLEquation:         description.HTMLEquation,
		IsChemicallyBalanced: description.IsChemicallyBalanced,
		IsTransport:          description.IsTransport,
		Ec:                   description.EC.Resource,
		Citations:            description.CitationStrings(),
		Substrates:           description.SubstrateAccessionIDs(),
		Products:             description.ProductAccessionIDs(),
		SubstrateOrProducts:  description.SubstrateOrProductAccessionIDs(),
		Location:             description.Location.Resource}
}

// NewCompound returns a Compound.
func NewCompound(description Description, subclass Subclass) Compound {
	var newCompound Compound
	compoundType := subclass.Resource[23:]
	switch subclass.Resource {
	case "http://rdf.rhea-db.org/SmallMolecule", "http://rdf.rhea-db.org/Polymer":
		newCompound = Compound{
			ID:        description.ID,
			Accession: description.About,
			Position:  description.Position,
			Name:      description.Name,
			HTMLName:  description.HTMLName,
			Formula:   description.Formula,
			Charge:    description.Charge,
			ChEBI:     description.ChEBI.Resource,

			CompoundID:        description.ID,
			CompoundAccession: description.Accession,
			CompoundName:      description.Name,
			CompoundHTMLName:  description.HTMLName,
			CompoundType:      compoundType}
		if compoundType == "Polymer" {
			newCompound.ChEBI = description.UnderlyingChEBI.Resource
		}
		// Add subclass ChEBI
		for _, sc := range description.Subclass {
			if strings.Contains(sc.Resource, "CHEBI") {
				newCompound.SubclassOfChEBI = sc.Resource
			}
		}
	case "http://rdf.rhea-db.org/GenericPolypeptide", "http://rdf.rhea-db.org/GenericPolynucleotide", "http://rdf.rhea-db.org/GenericHeteropolysaccharide":
		newCompound = Compound{
			Accession:        description.About,
			CompoundID:       description.ID,
			CompoundName:     description.Name,
			CompoundHTMLName: description.HTMLName,
			CompoundType:     compoundType}
	}
	return newCompound
}

// NewReactionParticipant returns a ReactionParticipant.
func NewReactionParticipant(description Description, containsx ContainsX, compoundParticipantMap map[string]string) (ReactionParticipant, error) {
	// Get reaction sides
	// gzip -d -k -c rhea.rdf.gz | grep -o -P '(?<=contains).*(?= rdf)' | tr ' ' '\n' | sort -u | tr '\n' ' '
	// The exceptions to numeric contains are 2n, N, Nminus1, and Nplus1
	var newReactionParticipant ReactionParticipant
	switch containsx.XMLName.Local {
	case "containsN":
		newReactionParticipant = ReactionParticipant{
			ReactionSide: description.About,
			Contains:     1,
			ContainsN:    true,
			Minus:        false,
			Plus:         false,
			Accession:    containsx.Content}
	case "contains2n":
		newReactionParticipant = ReactionParticipant{
			ReactionSide: description.About,
			Contains:     2,
			ContainsN:    true,
			Minus:        false,
			Plus:         false,
			Accession:    containsx.Content}
	case "containsNminus1":
		newReactionParticipant = ReactionParticipant{
			ReactionSide: description.About,
			Contains:     1,
			ContainsN:    true,
			Minus:        true,
			Plus:         false,
			Accession:    containsx.Content}
	case "containsNplus1":
		newReactionParticipant = ReactionParticipant{
			ReactionSide: description.About,
			Contains:     1,
			ContainsN:    true,
			Minus:        false,
			Plus:         true,
			Accession:    containsx.Content}
	default:
		i, err := strconv.Atoi(containsx.XMLName.Local[8:])
		if err != nil {
			return ReactionParticipant{}, err
		}
		newReactionParticipant = ReactionParticipant{
			ReactionSide: description.About,
			Contains:     i,
			ContainsN:    false,
			Minus:        false,
			Plus:         false,
			Accession:    containsx.Content}
	}
	newReactionParticipant.Compound = compoundParticipantMap[description.About]
	return newReactionParticipant, nil
}

// ParseRhea parses a list of bytes into a higher-level Rhea Struct.
func ParseRhea(rheaBytes []byte) (Rhea, error) {
	var err error
	// Read rheaBytes into a Rdf object
	var rdf Rdf
	err = xml.Unmarshal(rheaBytes, &rdf)
	if err != nil {
		return Rhea{}, err
	}

	// Initialize Rhea
	var rhea Rhea
	compoundParticipantMap := make(map[string]string)
	compoundMap := make(map[string]Compound)

	for _, description := range rdf.Descriptions {
		// Handle the case of a single compound -> reactive part, such as
		// <rdf:Description rdf:about="http://rdf.rhea-db.org/Compound_10594">
		// 	<rh:reactivePart rdf:resource="http://rdf.rhea-db.org/Compound_10594_rp2"/>
		// </rdf:Description>
		if (len(description.Subclass) == 0) && (description.ReactivePartXML.Resource != "") {
			compoundParticipantMap[description.ReactivePartXML.Resource] = description.About
		}
		if description.Compound.Resource != "" {
			compoundParticipantMap[description.About] = description.Compound.Resource
		}

		for _, subclass := range description.Subclass {
			switch subclass.Resource {
			case "http://rdf.rhea-db.org/BidirectionalReaction", "http://rdf.rhea-db.org/DirectionalReaction":
				newReaction := NewReaction(description, subclass)
				rhea.Reactions = append(rhea.Reactions, newReaction)
			case "http://rdf.rhea-db.org/SmallMolecule", "http://rdf.rhea-db.org/Polymer":
				newCompound := NewCompound(description, subclass)
				rhea.Compounds = append(rhea.Compounds, newCompound)
			case "http://rdf.rhea-db.org/GenericPolypeptide", "http://rdf.rhea-db.org/GenericPolynucleotide", "http://rdf.rhea-db.org/GenericHeteropolysaccharide":
				newCompound := NewCompound(description, subclass)
				compoundMap[description.About] = newCompound
				compoundParticipantMap[description.ReactivePartXML.Resource] = description.About
			}
		}
	}

	// Go back and get the ReactiveParts
	for _, description := range rdf.Descriptions {
		for _, subclass := range description.Subclass {
			if subclass.Resource == "http://rdf.rhea-db.org/ReactivePart" {
				newCompound, ok := compoundMap[compoundParticipantMap[description.About]]
				if ok != true {
					return Rhea{}, errors.New("Could not find " + description.About)
				}
				newCompound.ID = description.ID
				newCompound.CompoundAccession = description.About
				newCompound.Position = description.Position
				newCompound.Name = description.Name
				newCompound.HTMLName = description.HTMLName
				newCompound.Formula = description.Formula
				newCompound.Charge = description.Charge
				newCompound.ChEBI = description.ChEBI.Resource
				rhea.Compounds = append(rhea.Compounds, newCompound)
			}
		}
		for _, containsx := range description.ContainsX {
			if !strings.Contains(containsx.XMLName.Local, "contains") {
				continue
			}
			newReactionParticipant, err := NewReactionParticipant(description, containsx, compoundParticipantMap)
			if err != nil {
				return Rhea{}, err
			}
			rhea.ReactionParticipants = append(rhea.ReactionParticipants, newReactionParticipant)
		}
	}
	return rhea, nil
}

// ReadRhea reads in a a gzip'd Rhea dump (https://www.rhea-db.org/help/download) into bytes.
func ReadGzippedXml(gzipPath string) ([]byte, error) {
	// Get gz'd file bytes
	xmlFile, err := os.Open(gzipPath)
	if err != nil {
		return []byte{}, err
	}

	// Decompress gz'd file
	r, err := gzip.NewReader(xmlFile)
	if err != nil {
		return []byte{}, err
	}

	// Read decompressed gz'd file
	rheaBytes, err := ioutil.ReadAll(r)
	if err != nil {
		return []byte{}, err
	}
	return rheaBytes, nil
}

/******************************************************************************

JSON functions

These functions simply export the Rhea database as a JSON file for consumption
in different programs.

******************************************************************************/

// Export exports Rhea as a JSON file
func (rhea *Rhea) Export() ([]byte, error) {
	rheaJSON, err := json.Marshal(rhea)
	if err != nil {
		return []byte{}, err
	}
	return rheaJSON, nil
}

/******************************************************************************

Rhea2Uniprot tsv

Rhea conveniently provides a TSV list of reaction IDs to Uniprot accession numbers.
These can be used to map Rhea reactions to protein sequences, which is very useful
for actually building synthetic circuits. While these parsers are very basic TSV
parsing, they are convenient to use with the rhea2uniprot data dumps specifically.

The sprot and trembl data dumps are available at:
https://ftp.expasy.org/databases/rhea/tsv/rhea2uniprot_sprot.tsv
https://ftp.expasy.org/databases/rhea/tsv/rhea2uniprot_trembl.tsv.gz

The decompressed gzip side of the rhea2uniprot_trembl.tsv.gz file, as of Apr 30, 2021
is 678M. Because it is fairly large, the parser is implemented to stream in the file,
passing the parsed lines into a channel, to save on memory.

The TSV is structured like:

`
RHEA_ID DIRECTION       MASTER_ID       ID
10008   UN      10008   O17433
10008   UN      10008   O34564
...
`

******************************************************************************/

// Rhea2Uniprot represents a single line of the TSV file.
type Rhea2Uniprot struct {
	RheaID    int
	Direction string
	MasterID  int
	UniprotID string
}

// ParseRhea2UniprotTsv parses a rhea2uniprot TSV file and sends values to a channel.
func ParseRhea2UniprotTsv(r io.Reader, lines chan<- Rhea2Uniprot) {
	start := true
	scanner := bufio.NewScanner(r)
	for scanner.Scan() {
		// We skip the header line
		if start {
			start = false
			continue
		}

		// Get the line from the scanner
		line := scanner.Text()

		// Split the line between tabs
		lineSplit := strings.Split(line, "\t")

		// Send line to lines channel
		rheaId, err := strconv.Atoi(lineSplit[0])
		if err != nil {
			panic(err)
		}
		masterId, err := strconv.Atoi(lineSplit[2])
		if err != nil {
			panic(err)
		}
		lines <- Rhea2Uniprot{RheaID: rheaId, Direction: lineSplit[1], MasterID: masterId, UniprotID: lineSplit[3]}
	}
	close(lines)
}

// ReadRhea2UniprotSprot reads in the rhea2uniprot sprot TSV file (not gzipped) and sends values into a channel.
func ReadRhea2UniprotSprot(path string, lines chan Rhea2Uniprot) {
	file, _ := os.Open(path)
	go ParseRhea2UniprotTsv(file, lines)
}

// ReadRhea2UniprotTrembl reads in the rhea2uniprot trembl TSV file (gzipped) and sends values into channel.
func ReadRhea2UniprotTrembl(path string, lines chan Rhea2Uniprot) {
	file, _ := os.Open(path)
	r, _ := gzip.NewReader(file)
	go ParseRhea2UniprotTsv(r, lines)
}
