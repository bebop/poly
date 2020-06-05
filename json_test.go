package main

import (
	"testing"

	"github.com/google/go-cmp/cmp"
)

func TestJSONIO(t *testing.T) {
	testSequence := ReadGbk("data/bsub.gbk")
	WriteJSON(testSequence, "data/test.json")
	readTestSequence := ReadJSON("data/test.json")
	if diff := cmp.Diff(testSequence, readTestSequence); diff != "" {
		t.Errorf("MakeGatewayInfo() mismatch (-want +got):\n%s", diff)
	}
}
