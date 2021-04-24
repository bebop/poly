package rhea

import (
	"compress/gzip"
	"database/sql"
	"encoding/xml"
	"errors"
	"io/ioutil"
	"os"
	"strconv"
	"strings"
)

/******************************************************************************

Lower level structs

These structs operate at the lowest level to parse the RheaRdf database dump
into structs that Golang can understand. Getting all of Rhea from Rdf->Golang
is quite verbose, and so most of the time you should not use these structs unless
you know what you are doing and it cannot be accomplished with higher level structs.

******************************************************************************/

// RheaRdf is the XML representation of the Rhea database.
type RheaRdf struct {
	XMLName      xml.Name      `xml:"RDF"`
	Descriptions []Description `xml:"Description"`
}

// Description is the XML representation of an object in the Rhea database.
type Description struct {
	// Reaction
	XMLName              xml.Name             `xml:"Description"`
	About                string               `xml:"about,attr"`
	Id                   int                  `xml:"id"`
	Accession            string               `xml:"accession"`
	Equation             string               `xml:"equation"`
	HtmlEquation         string               `xml:"htmlEquation"`
	IsChemicallyBalanced bool                 `xml:"isChemicallyBalanced"`
	IsTransport          bool                 `xml:"isTransport"`
	Citations            []Citation           `xml:"citation"`
	Substrates           []Substrate          `xml:"substrates"`
	Products             []Product            `xml:"products"`
	SubstrateOrProducts  []SubstrateOrProduct `xml:"substratesOrProducts"`
	Subclass             []Subclass           `xml:"subClassOf"`
	Comment              string               `xml:"comment"`
	EC                   EC                   `xml:"ec"`
	Status               Status               `xml:"status"`
	Compound             CompoundXml          `xml:"compound"`

	// ReactionSide / Reaction Participant
	BidirectionalReactions []BidirectionalReaction `xml:"bidirectionalReaction"`
	DirectionalReactions   []DirectionalReaction   `xml:"directionalReaction"`
	Side                   Side                    `xml:"side"`
	SeeAlsos               SeeAlso                 `xml:"seeAlso"`
	TransformableTo        TransformableTo         `xml:"transformableTo"`
	CuratedOrder           int                     `xml:"curatedOrder"`
	Contains               Contains                `xml:"contains"`

	// ContainsX contains all other name-attribute pairs, with names like "contains1" in mind
	ContainsX []ContainsX `xml:",any"`

	// Small Molecule tags
	Name     string   `xml:"name"`
	HtmlName string   `xml:"htmlName"`
	Formula  string   `xml:"formula"`
	Charge   string   `xml:"charge"`
	Chebi    ChebiXml `xml:"chebi"`

	// Generic Compound
	ReactivePartXml ReactivePartXml `xml:"reactivePart"`

	// ReactivePart
	Position string `xml:"position"`

	// Polymer
	UnderlyingChebi     UnderlyingChebi `xml:"underlyingChebi"`
	PolymerizationIndex string          `xml:"polymerizationIndex"`

	// Transport
	Location Location `xml:"location"`
}

// CitationStrings gives a list of citation strings from a description.
func (d *Description) CitationStrings() []string {
	var output []string
	for _, x := range d.Citations {
		output = append(output, x.Resource)
	}
	return output
}

// SubstrateStrings gives a list of substrate accessions from a description.
func (d *Description) SubstrateStrings() []string {
	var output []string
	for _, x := range d.Substrates {
		output = append(output, x.Resource)
	}
	return output
}

// ProductStrings gives a list of product accessions from a description.
func (d *Description) ProductStrings() []string {
	var output []string
	for _, x := range d.Products {
		output = append(output, x.Resource)
	}
	return output
}

// SubstrateOrProductStrings gives a list of substrateOrProduct accessions from a description.
func (d *Description) SubstrateOrProductStrings() []string {
	var output []string
	for _, x := range d.SubstrateOrProducts {
		output = append(output, x.Resource)
	}
	return output
}

// Citation is an XML representation of a citation of a description.
type Citation struct {
	XMLName  xml.Name `xml:"citation"`
	Resource string   `xml:"resource,attr"`
}

// Substrate is an XML representation of a substrate.
type Substrate struct {
	XMLName  xml.Name `xml:"substrates"`
	Resource string   `xml:"resource,attr"`
}

// Product is an XML representation of a product.
type Product struct {
	XMLName  xml.Name `xml:"products"`
	Resource string   `xml:"resource,attr"`
}

// SubstrateOrProduct is an XML representation of a SubstrateOrProduct.
type SubstrateOrProduct struct {
	XMLName  xml.Name `xml:"substratesOrProducts"`
	Resource string   `xml:"resource,attr"`
}

// Subclass is an XML representation of a subclass, which can mean many different things in Rhea.
type Subclass struct {
	XMLName  xml.Name `xml:"subClassOf"`
	Resource string   `xml:"resource,attr"`
}

// EC is an XML representation of an EC number (Enzyme Commission Number) of a description.
type EC struct {
	XMLName  xml.Name `xml:"ec"`
	Resource string   `xml:"resource,attr"`
}

// Status is an XML representation of the current status of an description.
type Status struct {
	XMLName  xml.Name `xml:"status"`
	Resource string   `xml:"resource,attr"`
}

// CompoundXml is an XML representation of a compound.
type CompoundXml struct {
	XMLName  xml.Name `xml:"compound"`
	Resource string   `xml:"resource,attr"`
}

// BidirectionalReaction is an XML representation of a Rhea bidirectional reaction.
type BidirectionalReaction struct {
	XMLName  xml.Name `xml:"bidirectionalReaction"`
	Resource string   `xml:"resource,attr"`
}

// DirectionalReaction is an XML representation of a Rhea directional reaction.
type DirectionalReaction struct {
	XMLName  xml.Name `xml:"directionalReaction"`
	Resource string   `xml:"resource,attr"`
}

// Side is an XML representation of a Rhea ReactionSide.
type Side struct {
	XMLName  xml.Name `xml:"side"`
	Resource string   `xml:"resource,attr"`
}

// SeeAlso is an XML representation of a SeeAlso XML in a description.
type SeeAlso struct {
	XMLName  xml.Name `xml:"seeAlso"`
	Resource string   `xml:"resource,attr"`
}

// TransformableTo is an XML representation of a transformableTo in a description. This essentially links two ReactionSides in Rhea.
type TransformableTo struct {
	XMLName  xml.Name `xml:"transformableTo"`
	Resource string   `xml:"resource,attr"`
}

// Contains is an XML representation of what Compound a ReactionParticipant contains.
type Contains struct {
	XMLName  xml.Name `xml:"contains"`
	Resource string   `xml:"resource,attr"`
}

// ContainsX is a catch-all XML representation of how many compounds a ReactionParticipant would use.
type ContainsX struct {
	XMLName xml.Name
	Content string `xml:"resource,attr"`
}

// ChebiXml is an XML representation of a Chebi.
type ChebiXml struct {
	XMLName  xml.Name `xml:"chebi"`
	Resource string   `xml:"resource,attr"`
}

// UnderlyingChebi is an XML representation of Chebi that builds a Polymer.
type UnderlyingChebi struct {
	XMLName  xml.Name `xml:"underlyingChebi"`
	Resource string   `xml:"resource,attr"`
}

// ReactivePartXml is an XML representation of a ReactivePart.
type ReactivePartXml struct {
	XMLName  xml.Name `xml:"reactivePart"`
	Resource string   `xml:"resource,attr"`
}

// Location is an XML representation of Locations in a cell, usually referring to transport enzymes.
type Location struct {
	XMLName  xml.Name `xml:"location"`
	Resource string   `xml:"resource,attr"`
}

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
type Compound struct {
	Id                  int    `json:"id" db:"id"`
	Accession           string `json:"accession" db:"accession"`
	Position            string `json:"position" db: "position"`
	Name                string `json:"name" db:"name"`
	HtmlName            string `json:"htmlName" db:"htmlname"`
	Formula             string `json:"formula" db:"formula"`
	Charge              string `json:"charge" db:"charge"`
	Chebi               string `json:"chebi" db:"chebi"`
	SubclassOfChebi     string `json:"subclassOfChebi"`
	PolymerizationIndex string `json:"polymerizationIndex" db:"polymerizationindex"`
	CompoundId          int    `json:"id" db:"compoundid"`
	CompoundAccession   string `json:"accession" db:"compoundaccession"`
	CompoundName        string `json:"name" db:"compoundname"`
	CompoundHtmlName    string `json:"htmlName" db:"compoundhtmlname"`
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
	Id                   int      `json:"id" db:"id"`
	Directional          bool     `json:"directional" db:"directional"`
	Accession            string   `json:"accession" db:"accession"`
	Status               string   `json:"status" db:"status"`
	Comment              string   `json:"comment" db:"comment"`
	Equation             string   `json:"equation" db:"equation"`
	HtmlEquation         string   `json:"htmlequation" db:"htmlequation"`
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

// ParseRhea parses a list of bytes into a higher-level Rhea Struct.
func ParseRhea(rheaBytes []byte) (Rhea, error) {
	var err error
	// Read rheaBytes into a RheaRdf object
	var rdf RheaRdf
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
		if (len(description.Subclass) == 0) && (description.ReactivePartXml.Resource != "") {
			compoundParticipantMap[description.ReactivePartXml.Resource] = description.About
		}
		if description.Compound.Resource != "" {
			compoundParticipantMap[description.About] = description.Compound.Resource
		}

		for _, subclass := range description.Subclass {
			switch subclass.Resource {
			case "http://rdf.rhea-db.org/DirectionalReaction":
				newReaction := Reaction{
					Id:                   description.Id,
					Directional:          true,
					Accession:            description.Accession,
					Status:               description.Status.Resource,
					Comment:              description.Comment,
					Equation:             description.Equation,
					HtmlEquation:         description.HtmlEquation,
					IsChemicallyBalanced: description.IsChemicallyBalanced,
					IsTransport:          description.IsTransport,
					Ec:                   description.EC.Resource,
					Citations:            description.CitationStrings(),
					Substrates:           description.SubstrateStrings(),
					Products:             description.ProductStrings(),
					SubstrateOrProducts:  description.SubstrateOrProductStrings(),
					Location:             description.Location.Resource}
				rhea.Reactions = append(rhea.Reactions, newReaction)
			case "http://rdf.rhea-db.org/BidirectionalReaction":
				newReaction := Reaction{
					Id:                   description.Id,
					Directional:          false,
					Accession:            description.Accession,
					Status:               description.Status.Resource,
					Comment:              description.Comment,
					Equation:             description.Equation,
					HtmlEquation:         description.HtmlEquation,
					IsChemicallyBalanced: description.IsChemicallyBalanced,
					IsTransport:          description.IsTransport,
					Ec:                   description.EC.Resource,
					Citations:            description.CitationStrings(),
					Substrates:           description.SubstrateStrings(),
					Products:             description.ProductStrings(),
					SubstrateOrProducts:  description.SubstrateOrProductStrings(),
					Location:             description.Location.Resource}
				rhea.Reactions = append(rhea.Reactions, newReaction)
			case "http://rdf.rhea-db.org/SmallMolecule", "http://rdf.rhea-db.org/Polymer":
				compoundType := subclass.Resource[23:]
				newCompound := Compound{
					Id:        description.Id,
					Accession: description.About,
					Position:  description.Position,
					Name:      description.Name,
					HtmlName:  description.HtmlName,
					Formula:   description.Formula,
					Charge:    description.Charge,
					Chebi:     description.Chebi.Resource,

					CompoundId:        description.Id,
					CompoundAccession: description.Accession,
					CompoundName:      description.Name,
					CompoundHtmlName:  description.HtmlName,
					CompoundType:      compoundType}
				if compoundType == "Polymer" {
					newCompound.Chebi = description.UnderlyingChebi.Resource
				}
				// Add subclass Chebi
				for _, sc := range description.Subclass {
					if strings.Contains(sc.Resource, "CHEBI") {
						newCompound.SubclassOfChebi = sc.Resource
					}
				}
				// Add new reactive parts and new compounds to rhea
				rhea.Compounds = append(rhea.Compounds, newCompound)
			case "http://rdf.rhea-db.org/GenericPolypeptide", "http://rdf.rhea-db.org/GenericPolynucleotide", "http://rdf.rhea-db.org/GenericHeteropolysaccharide":
				compoundType := subclass.Resource[23:]
				newCompound := Compound{
					Accession:        description.About,
					CompoundId:       description.Id,
					CompoundName:     description.Name,
					CompoundHtmlName: description.HtmlName,
					CompoundType:     compoundType}
				compoundMap[description.About] = newCompound
				compoundParticipantMap[description.ReactivePartXml.Resource] = description.About
			}
		}
	}

	// Go back and get the ReactiveParts
	for _, description := range rdf.Descriptions {
		for _, subclass := range description.Subclass {
			switch subclass.Resource {
			case "http://rdf.rhea-db.org/ReactivePart":
				newCompound, ok := compoundMap[compoundParticipantMap[description.About]]
				if ok != true {
					return Rhea{}, errors.New("Could not find " + description.About)
				}
				newCompound.Id = description.Id
				newCompound.CompoundAccession = description.About
				newCompound.Position = description.Position
				newCompound.Name = description.Name
				newCompound.HtmlName = description.HtmlName
				newCompound.Formula = description.Formula
				newCompound.Charge = description.Charge
				newCompound.Chebi = description.Chebi.Resource
				rhea.Compounds = append(rhea.Compounds, newCompound)
			}
		}
		for _, containsx := range description.ContainsX {
			if strings.Contains(containsx.XMLName.Local, "contains") {
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
						return Rhea{}, err
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
				rhea.ReactionParticipants = append(rhea.ReactionParticipants, newReactionParticipant)
			}
		}
	}
	return rhea, nil
}

// ReadRhea reads in a a gzip'd Rhea dump (https://www.rhea-db.org/help/download) into bytes.
func ReadRhea(gzipPath string) ([]byte, error) {
	// Get gz'd file bytes
	xmlFile, err := os.Open("data/rhea.rdf.gz")
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

Database functions

These functions take in the rhea.rdf.gz dump file and import them into an sqlite
database. This database can be used locally for mapping functions, but more importantly,
it is used because it can enforce the relationships between each part of Rhea and alert
the user if the parser is failing to pick up on any of the relationships in the Rhea
database.

******************************************************************************/

// CreateRheaSqlite initializes a SQLite database with the proper tables to have Rhea inserted using InsertRheaSqlite.
func CreateRheaSqlite(db *sql.DB) error {
	schema := `
	-- Rhea tables --
	CREATE TABLE chebi (
		accession TEXT PRIMARY KEY,
		subclassof TEXT REFERENCES chebi(accession)
	);
	
	CREATE TABLE compound (
		id INT,
		accession TEXT PRIMARY KEY,
		position TEXT,
		name TEXT,
		htmlname TEXT,
		formula TEXT,
		charge TEXT,
		chebi TEXT REFERENCES chebi(accession),
		polymerizationindex TEXT,
		compoundtype TEXT NOT NULL CHECK(compoundtype IN ('SmallMolecule', 'Polymer', 'GenericPolypeptide', 'GenericPolynucleotide', 'GenericHeteropolysaccharide'))
	);
	
	CREATE TABLE reactivepart (
		id INT,
		accession TEXT PRIMARY KEY,
		name TEXT,
		htmlname TEXT,
		compound TEXT NOT NULL REFERENCES compound(accession)
	);
	
	CREATE TABLE reaction (
		id INT,
		directional BOOL NOT NULL DEFAULT false,
		accession TEXT PRIMARY KEY,
		status TEXT,
		comment TEXT,
		equation TEXT,
		htmlequation TEXT,
		ischemicallybalanced BOOL NOT NULL DEFAULT true,
		istransport BOOL NOT NULL DEFAULT false,
		ec TEXT,
		location TEXT
	);
	
	CREATE TABLE reactionside (
		accession TEXT PRIMARY KEY
	);
	
	CREATE TABLE reactionsidereaction (
		reaction TEXT NOT NULL REFERENCES reaction(accession),
		reactionside TEXT NOT NULL REFERENCES reactionside(accession),
		substrateorproduct BOOL NOT NULL DEFAULT false, -- If set to true, disregard substrate below
	        substrate BOOL NOT NULL DEFAULT false -- If substrateOrProduct is false, "substrate = true" means this is a substrate and "substrate = false" means this is a product
	);
	
	CREATE TABLE reactionparticipant (
		compound TEXT REFERENCES compound(accession),
	        reactionside TEXT NOT NULL REFERENCES reactionside(accession),
	        contains INTEGER,
	        containsn BOOL NOT NULL DEFAULT false,
	        minus BOOL NOT NULL DEFAULT false,
	        plus BOOL NOT NULL DEFAULT false
	);
	`
	_, err := db.Exec(schema)
	if err != nil {
		return errors.New("Failed to exec schema. Failed on: " + err.Error())
	}
	return nil
}

// InsertRheaSqlite inserts the Rhea database into an SQLite database with proper normalization for advanced queries.
func InsertRheaSqlite(db *sql.DB, rheaPath string) error {

	// Read the Compressed Rhea XML to bytes
	rheaBytes, err := ReadRhea(rheaPath)
	if err != nil {
		return err
	}

	// Parse the Rhea bytes into the rhea struct
	rhea, err := ParseRhea(rheaBytes)
	if err != nil {
		return err
	}

	// Start transaction with database for insertion. This ensures if there are any problems, they are seamlessly rolled back
	tx, err := db.Begin()
	if err != nil {
		return err
	}

	// First, insert Chebis and Compounds
	compoundKeys := make(map[string]bool)
	for _, compound := range rhea.Compounds {
		// Insert root chebi. Ie, what this current compound's subclass is
		_, err = tx.Exec("INSERT OR IGNORE INTO chebi(accession) VALUES (?)", compound.SubclassOfChebi)
		if err != nil {
			tx.Rollback()
			return err
		}
		// Insert the chebi of the current compound. If it is already inserted, update the subclassification
		_, err = tx.Exec("INSERT INTO chebi(accession, subclassof) VALUES (?, ?) ON CONFLICT (accession) DO UPDATE SET subclassof = ?", compound.Chebi, compound.SubclassOfChebi, compound.SubclassOfChebi)
		if err != nil {
			tx.Rollback()
			return err
		}
		// Insert the compound itself
		_, err = tx.Exec("INSERT INTO compound(id, accession, position, name, htmlname, formula, charge, chebi, polymerizationindex, compoundtype) VALUES (?, ?, ?, ?, ?, ?, ?, ?, ?, ?) ON CONFLICT DO NOTHING", compound.Id, compound.Accession, compound.Position, compound.Name, compound.HtmlName, compound.Formula, compound.Charge, compound.Chebi, compound.PolymerizationIndex, compound.CompoundType)
		if err != nil {
			tx.Rollback()
			return err
		}
		// If the compound isn't a small molecule or polymer, that means it would be a reactive part of a larger compound. So we add it in
		if (compound.CompoundType != "SmallMolecule") && (compound.CompoundType != "Polymer") {
			_, err = tx.Exec("INSERT INTO reactivepart(id, accession, name, htmlname, compound) VALUES (?, ?, ?, ?, ?)", compound.CompoundId, compound.CompoundAccession, compound.CompoundName, compound.CompoundHtmlName, compound.Accession)
			if err != nil {
				tx.Rollback()
				return err
			}
		}
		// Add compound.Access to the compoundKeys map
		compoundKeys[compound.Accession] = true
	}

	// Next, insert the ReactionSides and ReactionParticipants
	for _, reactionParticipant := range rhea.ReactionParticipants {
		// Insert ReactionSide, which is needed to insert the ReactionParticipant
		_, err = tx.Exec("INSERT INTO reactionside(accession) VALUES (?) ON CONFLICT DO NOTHING", reactionParticipant.ReactionSide)
		if err != nil {
			tx.Rollback()
			return err
		}
		// Insert the ReactionParticipants
		_, err = tx.Exec("INSERT INTO reactionparticipant(reactionside, contains, containsn, minus, plus, compound) VALUES (?, ?, ?, ?, ?, ?)", reactionParticipant.ReactionSide, reactionParticipant.Contains, reactionParticipant.ContainsN, reactionParticipant.Minus, reactionParticipant.Plus, reactionParticipant.Compound)
		if err != nil {
			tx.Rollback()
			return err
		}
	}

	// Insert the Reactions themselves
	for _, reaction := range rhea.Reactions {
		// Insert Reaction
		_, err = tx.Exec("INSERT INTO reaction(id, directional, accession, status, comment, equation, htmlequation, ischemicallybalanced, istransport, ec, location) VALUES (?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?)", reaction.Id, reaction.Directional, reaction.Accession, reaction.Status, reaction.Comment, reaction.Equation, reaction.HtmlEquation, reaction.IsChemicallyBalanced, reaction.IsTransport, reaction.Ec, reaction.Location)
		if err != nil {
			tx.Rollback()
			return err
		}

		// Insert ReactionsideReaction. Basically, these represent the substrates, products, or substratesOrProducts of a given reaction
		for _, substrate := range reaction.Substrates {
			_, err = tx.Exec("INSERT INTO reactionsidereaction(reaction, reactionside, substrate) VALUES (?, ?, true)", reaction.Accession, substrate)
			if err != nil {
				tx.Rollback()
				return err
			}
		}
		for _, product := range reaction.Products {
			_, err = tx.Exec("INSERT INTO reactionsidereaction(reaction, reactionside) VALUES (?, ?)", reaction.Accession, product)
			if err != nil {
				tx.Rollback()
				return err
			}
		}
		for _, substrateorproduct := range reaction.SubstrateOrProducts {
			_, err = tx.Exec("INSERT INTO reactionsidereaction(reaction, reactionside, substrateorproduct) VALUES (?, ?, true)", reaction.Accession, substrateorproduct)
			if err != nil {
				tx.Rollback()
				return err
			}
		}

	}

	// Commit
	err = tx.Commit()
	if err != nil {
		return err
	}
	return nil
}
