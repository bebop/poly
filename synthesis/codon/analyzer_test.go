package codon

import (
	"testing"

	"github.com/TimothyStiles/poly/io/genbank"
	"github.com/google/go-cmp/cmp"
	"github.com/google/go-cmp/cmp/cmpopts"
)

func TestComputeStats(t *testing.T) {
	t.Parallel()

	sortStrings := cmpopts.SortSlices((func(a, b string) bool { return a < b }))

	tests := []struct {
		name string

		testDataPath string
		wantStats    *Stats

		wantErr error
	}{
		{
			name: "pUC19 plasmid",

			testDataPath: "../../data/puc19.gbk",

			wantStats: &Stats{
				StartCodonCount: map[string]int{"atg": 2},
			},
			wantErr: nil,
		},
	}

	for _, tt := range tests {
		var tt = tt
		t.Run(tt.name, func(t *testing.T) {
			t.Parallel()
			sequence, err := genbank.Read(tt.testDataPath)
			if err != nil {
				t.Fatal(err)
			}

			analyzer := &Analyzer{}

			stats, err := analyzer.ComputeStats(sequence)
			if err != nil {
				t.Fatal(err)
			}

			if diff := cmp.Diff(stats, tt.wantStats); diff != "" {
				t.Errorf("stats and tt.wantStats didn't match (-got,+want): %s", diff)
			}
		})
	}
}
