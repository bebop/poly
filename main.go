package main

import (
	"log"
	"os"

	"github.com/urfave/cli/v2"
)

// Meta Holds all the meta information of an AnnotatedSequence struct.
type Meta struct {
	// shared
	Name        string
	GffVersion  string
	RegionStart int
	RegionEnd   int
	// genbank specific
	Size            int
	Type            string
	GenbankDivision string
	Date            string
	Definition      string
	Accession       string
	Version         string
	Keywords        string
	Organism        string
	Source          string
	Origin          string
	Locus           Locus
	References      []Reference
	Primaries       []Primary
}

// Primary Holds all the Primary information of a Meta struct.
type Primary struct {
	RefSeq, PrimaryIdentifier, Primary_Span, Comp string
}

// genbank specific
// type Reference struct {
// 	Authors []string
// 	Title   string
// 	Journal string
// 	PubMed  string
// }

// Reference holds information one reference in a Meta struct.
type Reference struct {
	Index, Authors, Title, Journal, PubMed, Remark, Range string
}

// Locus holds Locus information in a Meta struct.
type Locus struct {
	Name, SequenceLength, MoleculeType, GenBankDivision, ModDate string
	Circular                                                     bool
}

// Feature holds a single annotation in a struct. from https://github.com/blachlylab/gff3/blob/master/gff3.go
type Feature struct {
	Name string //Seqid in gff, name in gbk
	//gff specific
	Source     string
	Type       string
	Start      int
	End        int
	Score      float64
	Strand     byte
	Phase      int
	Attributes map[string]string // Known as "qualifiers" for gbk, "attributes" for gff.
	//gbk specific
	Location string
	Sequence string
}

// Sequence holds raw sequence information in an AnnotatedSequence struct.
type Sequence struct {
	Description string
	Sequence    string
}

// AnnotatedSequence holds all sequence information in a single struct.
type AnnotatedSequence struct {
	Meta     Meta
	Features []Feature
	Sequence Sequence
}

func main() {
	app := &cli.App{
		Name:  "poly",
		Usage: "A command line utility for engineering organisms.",
		Flags: []cli.Flag{
			&cli.BoolFlag{
				Name:  "y",
				Usage: "Answers yes for all prompts.",
			},
		},
		Commands: []*cli.Command{
			{
				Name:    "convert",
				Aliases: []string{"c"},
				Usage:   "Convert a single file or set of files from one type to another. Genbank to Json, Gff to genbank, etc.",
				Action: func(c *cli.Context) error {
					return nil
				}},
		}}

	err := app.Run(os.Args)
	if err != nil {
		log.Fatal(err)
	}
}
