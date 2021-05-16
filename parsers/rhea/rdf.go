package rhea

import (
	"encoding/xml"
)

/******************************************************************************

Lower level structs

These structs operate at the lowest level to parse the Rhea Rdf database dump
into structs that Golang can understand. Getting all of Rhea from Rdf->Golang
is quite verbose, and so most of the time you should not use these structs unless
you know what you are doing and it cannot be accomplished with higher level structs.

******************************************************************************/

// Rdf is the RDF XML representation of the Rhea database.
type Rdf struct {
	XMLName      xml.Name      `xml:"RDF"`
	Descriptions []Description `xml:"Description"`
}

// Description is the XML representation of an object in the Rhea database.
type Description struct {
	// Reaction
	XMLName              xml.Name             `xml:"Description"`
	About                string               `xml:"about,attr"`
	ID                   int                  `xml:"id"`
	Accession            string               `xml:"accession"` // Accession refers to a Rhea accession number
	Equation             string               `xml:"equation"`
	HTMLEquation         string               `xml:"htmlEquation"`
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
	Compound             CompoundXML          `xml:"compound"`

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
	HTMLName string   `xml:"htmlName"`
	Formula  string   `xml:"formula"`
	Charge   string   `xml:"charge"`
	ChEBI    ChEBIXML `xml:"chebi"`

	// Generic Compound
	ReactivePartXML ReactivePartXML `xml:"reactivePart"`

	// ReactivePart
	Position string `xml:"position"`

	// Polymer
	UnderlyingChEBI     UnderlyingChEBI `xml:"underlyingChEBI"`
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

// SubstrateAccessionIDs gives a list of substrate accessions from a description.
func (d *Description) SubstrateAccessionIDs() []string {
	var output []string
	for _, x := range d.Substrates {
		output = append(output, x.Resource)
	}
	return output
}

// ProductAccessionIDs gives a list of product accessions from a description.
func (d *Description) ProductAccessionIDs() []string {
	var output []string
	for _, x := range d.Products {
		output = append(output, x.Resource)
	}
	return output
}

// SubstrateOrProductAccessionIDs gives a list of substrateOrProduct accessions from a description.
func (d *Description) SubstrateOrProductAccessionIDs() []string {
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

// CompoundXML is an XML representation of a compound.
type CompoundXML struct {
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

// ChEBIXML is an XML representation of a ChEBI.
type ChEBIXML struct {
	XMLName  xml.Name `xml:"chebi"`
	Resource string   `xml:"resource,attr"`
}

// UnderlyingChEBI is an XML representation of ChEBI that builds a Polymer.
type UnderlyingChEBI struct {
	XMLName  xml.Name `xml:"underlyingChEBI"`
	Resource string   `xml:"resource,attr"`
}

// ReactivePartXML is an XML representation of a ReactivePart.
type ReactivePartXML struct {
	XMLName  xml.Name `xml:"reactivePart"`
	Resource string   `xml:"resource,attr"`
}

// Location is an XML representation of Locations in a cell, usually referring to transport enzymes.
type Location struct {
	XMLName  xml.Name `xml:"location"`
	Resource string   `xml:"resource,attr"`
}
