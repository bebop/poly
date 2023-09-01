package bio_test

import (
	"io"
	"testing"

	"github.com/TimothyStiles/poly/bio"
	"github.com/TimothyStiles/poly/bio/fasta"
)

func TestWriter(t *testing.T) {
	var _ bio.LowLevelParser[fasta.Record, fasta.Header] = &fasta.Parser{}
	var _ io.WriterTo = &fasta.Record{}
}
