package main

import (
	"testing"

	"github.com/google/go-cmp/cmp"
)

func testJSONIO(t *testing.T) {
	testSequence := parseGbk("data/test.gbk")
	writeJSON(testSequence, "data/test.json")
	readTestSequence := readJSON("data/test.json")
	if diff := cmp.Diff(testSequence, readTestSequence); diff != "" {
		t.Errorf("MakeGatewayInfo() mismatch (-want +got):\n%s", diff)
	}
}
