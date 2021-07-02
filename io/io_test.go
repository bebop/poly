package io_test

import (
	"testing"

	"github.com/TimothyStiles/poly"
	"github.com/TimothyStiles/poly/io/genbank"
	"github.com/TimothyStiles/poly/io/polyjson"
	"github.com/google/go-cmp/cmp"
	"github.com/google/go-cmp/cmp/cmpopts"
)

func TestLocusParseRegression(t *testing.T) {
	gbk := genbank.Read("../../data/puc19.gbk").Meta.Locus
	json := polyjson.Read("../../data/puc19static.json").Meta.Locus

	if diff := cmp.Diff(gbk, json, cmpopts.IgnoreFields(poly.Feature{}, "ParentSequence")); diff != "" {
		t.Errorf("The meta parser has changed behaviour. Got this diff:\n%s", diff)
	}
}
