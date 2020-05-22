package main

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

type Reference struct {
	Index, Authors, Title, Journal, PubMed, Remark string
	Start, End                                     int
}

type Locus struct {
	Name, SequenceLength, MoleculeType, GenBankDivision, ModDate string
	Circular                                                     bool
}

// from https://github.com/blachlylab/gff3/blob/master/gff3.go
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

type Sequence struct {
	Description string
	Sequence    string
}

type AnnotatedSequence struct {
	Meta     Meta
	Features []Feature
	Sequence Sequence
}

func main() {

	// fmt.Println(parseGff("data/ecoli-mg1655.gff"))
	// parseGbk("data/addgene-plasmid-50005-sequence-74677.gbk")
	parseGff("data/ecoli-mg1655.gff")
	// parseGbk("data/test.gbk")
	// parseQualifiersList()
}
