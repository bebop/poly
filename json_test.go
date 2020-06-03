package main

import (
	"testing"

	"github.com/google/go-cmp/cmp"
)

func testJSONIO(t *testing.T) {
	testSequence := ParseGbk("data/test.gbk")
	WriteJSON(testSequence, "data/test.json")
	readTestSequence := ReadJSON("data/test.json")
	if diff := cmp.Diff(testSequence, readTestSequence); diff != "" {
		t.Errorf("MakeGatewayInfo() mismatch (-want +got):\n%s", diff)
	}
}
