package genbank

import (
	"strings"
	"testing"
	"time"

	"github.com/google/go-cmp/cmp"
)

func TestParseHeader(t *testing.T) {
	testcases := []struct {
		name    string
		data    string
		want    Header
		wantErr bool
	}{{
		name: "parses example header from spec",
		data: `GBBCT1.SEQ          Genetic Sequence Data Bank
                          October 15 2023

                NCBI-GenBank Flat File Release 258.0

                     Bacterial Sequences (Part 1)

   51396 loci,    92682287 bases, from    51396 reported sequences
`,
		want: Header{
			FileName:     "GBBCT1.SEQ",
			Date:         time.Date(2023, time.October, 15, 0, 0, 0, 0, time.UTC),
			Title:        "Bacterial Sequences (Part 1)",
			NumEntries:   51396,
			NumBases:     92682287,
			NumSequences: 51396,
			MajorRelease: 258,
			MinorRelease: 0,
		},
	}, {
		name: "fails on premature EOF",
		data: `GBBCT1.SEQ          Genetic Sequence Data Bank
                          October 15 2023

                NCBI-GenBank Flat File Release`,
		wantErr: true,
	}}

	for _, tc := range testcases {
		t.Run(tc.name, func(t *testing.T) {
			p := NewParser(strings.NewReader(tc.data))
			got, err := p.parseHeader()

			if gotErr := err != nil; gotErr != tc.wantErr {
				t.Fatalf("incorrect error returned from parseHeader, wantErr: %v, err: %v", tc.wantErr, err)
			}

			if diff := cmp.Diff(tc.want, got); !tc.wantErr && diff != "" {
				t.Fatalf("parseHeader returned incorrect header, (-want, +got): %v", diff)
			}

		})
	}
}

func TestParseEntry(t *testing.T) {
	testcases := []struct {
		name           string
		data           string
		want           Entry
		wantReachedEOF bool
		wantErr        bool
	}{{
		name:           "parses LOCUS keyword",
		data:           "LOCUS       CP032762             5868661 bp    DNA     circular BCT 15-OCT-2018",
		wantReachedEOF: true,
		want: Entry{
			Name:            "CP032762",
			Length:          5868661,
			Strandedness:    Unknown,
			MoleculeType:    DNA,
			MoleculeToplogy: Circular,
			DivisionCode:    BacterialDivision,
			UpdateDate:      time.Date(2018, time.October, 15, 0, 0, 0, 0, time.UTC),
		},
	}, {
		name:           "parses DEFINITION keyword",
		data:           "DEFINITION test",
		wantReachedEOF: true,
		want:           Entry{Definition: "test"},
	}, {
		name:           "parses multiline DEFINITION keyword",
		data:           "DEFINITION test\n          another line\n\n          yet another line",
		wantReachedEOF: true,
		want:           Entry{Definition: "test\nanother line\n\nyet another line"},
	}, {
		name:           "parses repeated DEFINITION keywords",
		data:           "DEFINITION test\nDEFINITION another line\n\n          yet another line",
		wantReachedEOF: true,
		want:           Entry{Definition: "test\nanother line\n\nyet another line"},
	}, {
		name: "parses only first entry in reader",
		data: `DEFINITION test
//
DEFINITION another test`,
		want: Entry{Definition: "test"},
	}, {
		name:           "parses VERSION keyword",
		data:           "VERSION ASDF.130",
		want:           Entry{AccessionVersion: 130},
		wantReachedEOF: true,
	}, {
		name:           "parses ACCESSION keyword",
		data:           "ACCESSION HIIII",
		want:           Entry{Accession: "HIIII"},
		wantReachedEOF: true,
	}, {
		name:           "parses simple DBLINK keyword",
		data:           "DBLINK SomeDB: someReference",
		want:           Entry{DatabaseLinks: map[string][]string{"SomeDB": {"someReference"}}},
		wantReachedEOF: true,
	}, {
		name:           "parses DBLINK keyword with multiple refs for a single DB",
		data:           "DBLINK SomeDB: someReference, anotherReference",
		want:           Entry{DatabaseLinks: map[string][]string{"SomeDB": {"someReference", "anotherReference"}}},
		wantReachedEOF: true,
	}, {
		name:           "parses DBLINK keyword with multiple DBs",
		data:           "DBLINK SomeDB: someReference, anotherReference\n          AnotherDB: yetAnotherReference",
		want:           Entry{DatabaseLinks: map[string][]string{"SomeDB": {"someReference", "anotherReference"}, "AnotherDB": {"yetAnotherReference"}}},
		wantReachedEOF: true,
	}, {
		name:           "parses KEYWORDS keyword",
		data:           "KEYWORDS first; second; last.",
		want:           Entry{Keywords: []string{"first", "second", "last"}},
		wantReachedEOF: true,
	}, {
		name:    "KEYWORDS line must end with period",
		data:    "KEYWORDS first; second; last",
		wantErr: true,
	}, {
		name: "parses SOURCE keyword",
		data: `SOURCE common name (more info) molecule type.
  ORGANISM Scientific name.
          Taxon 1; Taxon 2; Taxon 3;
          Taxon 4; Taxon 5.
`,
		want: Entry{Source: Source{
			Name:           "common name (more info) molecule type",
			ScientificName: "Scientific name",
			Taxonomy:       []string{"Taxon 1", "Taxon 2", "Taxon 3", "Taxon 4", "Taxon 5"},
		}},
		wantReachedEOF: true,
	}}

	for _, tc := range testcases {
		t.Run(tc.name, func(t *testing.T) {
			p := NewParser(strings.NewReader(tc.data))
			got, gotReachedEOF, err := p.parseEntry()

			if gotErr := err != nil; gotErr != tc.wantErr {
				t.Fatalf("incorrect error returned from parseEntry, wantErr: %v, err: %v", tc.wantErr, err)
			}

			if gotReachedEOF != tc.wantReachedEOF {
				t.Fatalf("incorrect reached EOF status returned from parseEntry, want: %v, got: %v", tc.wantReachedEOF, gotReachedEOF)
			}

			if diff := cmp.Diff(tc.want, got); !tc.wantErr && diff != "" {
				t.Fatalf("parseEntry returned incorrect Entry, (-want, +got): %v", diff)
			}

		})
	}
}
