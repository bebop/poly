package bio_test

import (
	"testing"

	"github.com/TimothyStiles/poly/bio"
	"github.com/TimothyStiles/poly/bio/fasta"
)

func TestWriter(t *testing.T) {
	var _ bio.LowLevelParser[fasta.Fasta, fasta.Header] = &fasta.Parser{}
	var _ bio.Writer = &fasta.Fasta{}
}
