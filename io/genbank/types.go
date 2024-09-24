/*
Package genbank provides genbank parsers and writers.

GenBank is a flat text file format developed in the 1980s to annotate genetic
sequences, and has since become the standard for sharing annotated genetic
sequences. A full specification of the GenBank flatfile format can be found at
https://www.ncbi.nlm.nih.gov/genbank/release/current/.

This package provides a parser and writer to convert between the GenBank file
format and the more general GenBank struct.
*/
package genbank

import "time"

// Genbank is the main struct for the Genbank file format.
type Genbank struct {
	Header  Header
	Entries []Entry
}

// Header holds the information at the beginning of all
// Genbank files (see Genbank spec section 3.1)
type Header struct {
	FileName     string
	Title        string
	Date         time.Time
	MajorRelease int
	MinorRelease int
	NumEntries   int
	NumBases     int
	NumSequences int
}

// Strandedness represents whether a sequence is
// single-, double-, or mixed-stranded (see Genbank
// spec section 3.4.4.2).
type Strandedness string

const (
	DoubleStranded Strandedness = "ds"
	SingleStranded Strandedness = "ss"
	MixedStranded  Strandedness = "ms"
	Unknown        Strandedness = "unknown"
)

// MoleculeType represents what kind of molecule a sequence
// consists of (see Genbank spec section 3.4.4.2).
type MoleculeType string

const (
	NucleicAcid     MoleculeType = "NA"
	DNA             MoleculeType = "DNA"
	RNA             MoleculeType = "RNA"
	TransferRNA     MoleculeType = "tRNA"
	RibosomalRNA    MoleculeType = "rRNA"
	MessengerRNA    MoleculeType = "mRNA"
	SmallNuclearRNA MoleculeType = "uRNA"
	ViralCRNA       MoleculeType = "cRNA"
)

// MoleculeType represents the topology of a molecule
// (see Genbank spec section 3.4.4.2).
type MoleculeToplogy string

const (
	Circular MoleculeToplogy = "circular"
	Linear   MoleculeToplogy = "linear"
)

// DivisionCode represents which division of Genbank an entry
// belongs to (see Genbank spec section 3.4.4.2).
type DivisionCode string

const (
	PrimateDivision               DivisionCode = "PRI"
	RodentDivision                DivisionCode = "ROD"
	OtherMammalianDivision        DivisionCode = "MAM"
	OtherVertebrateDivision       DivisionCode = "VRT"
	InvertebrateDivision          DivisionCode = "INV"
	PlantFungalAlgalDivision      DivisionCode = "PLN"
	BacterialDivision             DivisionCode = "BCT"
	ViralDivision                 DivisionCode = "VRL"
	BacteriophageDivision         DivisionCode = "PHG"
	SyntheticDivision             DivisionCode = "SYN"
	UnnanotatedDivision           DivisionCode = "UNA"
	ExpressedSequenceTagDivision  DivisionCode = "EST"
	PatentDivision                DivisionCode = "PAT"
	SequenceTaggedSiteDivision    DivisionCode = "STS"
	GenomeSurveyDivision          DivisionCode = "GSS"
	HighThroughputGenomicDivision DivisionCode = "HTG"
	HighThroughputCDNADivision    DivisionCode = "HTC"
	EnvironmentalSamplingDivision DivisionCode = "ENV"
	ConstructedDivision           DivisionCode = "CON"
	TranscriptomeShotgunDivision  DivisionCode = "TSA"
)

type Range struct {
	From uint
	To   uint
}

// Reference is a reference to a publication (see
// Genbank spec section 3.4.11).
type Reference struct {
	Number      int
	BaseRange   Range
	Authors     string
	Title       string
	Medline     string
	PubMed      string
	PageNumbers string
	Remark      string
	Year        string
}

type Book struct {
	Reference
	Editors           []string
	BookTitle         string
	PublisherName     string
	PublisherLocation string
}

type Article struct {
	Reference
	Journal string
	Volume  string
}

// Source describes the source of an entry (see
// Genbank spec section 3.4.10).
type Source struct {
	Name           string
	ScientificName string
	Taxonomy       []string
}

// A Feature represents genes, gene products, and
// other areas of significance in an entry (see
// Genbank spec section 3.4.12).
type Feature struct {
	Key        string
	Location   Location
	Qualifiers map[string][]string
}

type Location string

// Entry represents a single entry in the Genbank file
// (see Genbank spec section 3.1).
type Entry struct {
	Name             string
	Length           int
	Strandedness     Strandedness
	MoleculeType     MoleculeType
	MoleculeToplogy  MoleculeToplogy
	DivisionCode     DivisionCode
	UpdateDate       time.Time
	Definition       string
	Accession        string
	AccessionVersion int
	CrossReferences  map[string][]string
	References       []Reference
	Source           Source
	Features         []Feature
	Origin           string
	Sequence         string
	DatabaseLinks    map[string][]string
	Keywords         []string
	Segment          int
	TotalSegments    int
}
