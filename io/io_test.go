package io

import (
	"testing"

	"github.com/TimothyStiles/poly/io/genbank"
	"github.com/TimothyStiles/poly/io/polyjson"
	"github.com/google/go-cmp/cmp"
	"github.com/google/go-cmp/cmp/cmpopts"
)

func TestLocusParseRegression(t *testing.T) {
	gbk, _ := genbank.Read("../../data/puc19.gbk")
	json, _ := polyjson.Read("../../data/puc19static.json")

	if diff := cmp.Diff(gbk, json, cmpopts.IgnoreFields(polyjson.Feature{}, "ParentSequence")); diff != "" {
		t.Errorf("The meta parser has changed behaviour. Got this diff:\n%s", diff)
	}
}
