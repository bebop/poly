package bio

import (
	"testing"

	"github.com/TimothyStiles/poly/bio/fasta"
	"github.com/TimothyStiles/poly/bio/fastq"
)

func TestWriter(t *testing.T) {
	var _ Writer = &fastq.Fastq{}
	var _ Writer = &fasta.Fasta{}
}
