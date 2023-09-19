package codon

import "github.com/TimothyStiles/poly/io/genbank"

// Stats denotes a set of computed codon table statistics
type Stats struct {
	StartCodonCount map[string]int `json:"start_codons_counts"`
}

func NewStats() *Stats {
	return &Stats{
		StartCodonCount: map[string]int{},
	}
}

// Analyzer is used to compute codon usage statistics
type Analyzer struct{}

// ComputeStats returns a set of codon related statistics from a genebank sequence
func (a *Analyzer) ComputeStats(sequence genbank.Genbank) (*Stats, error) {
	stats := NewStats()

	// iterate through the features of the genbank file and if the feature is a coding region, get the start codon
	for _, feature := range sequence.Features {
		if feature.Type == "CDS" {
			sequence, err := feature.GetSequence()
			if err != nil {
				return nil, err
			}

			stats.StartCodonCount[sequence[:3]]++
		}
	}

	return stats, nil
}
