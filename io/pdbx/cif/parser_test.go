package cif

import (
	"strings"
	"testing"

	"github.com/google/go-cmp/cmp"
)

func TestParser(t *testing.T) {
	testCases := []struct {
		name          string
		input         string
		want          CIF
		wantSyntaxErr bool
	}{{
		name:  "empty input yields empty CIF",
		input: "",
		want:  NewCIF(),
	}, {
		name:          "tag:value must be inside data block",
		input:         "_someTag someValue",
		wantSyntaxErr: true,
	}, {
		name:          "loop_ must be inside data block",
		input:         "loop_",
		wantSyntaxErr: true,
	}, {
		name:          "save frame header must be inside data block",
		input:         "save_asdf",
		wantSyntaxErr: true,
	}, {
		name: "data block header must have name",
		input: `
data_
_someTag someValue`,
		wantSyntaxErr: true,
	}, {
		name: "save frame header must have name",
		input: `
data_aBlock
save_
_someTag someValue
save_`,
		wantSyntaxErr: true,
	}, {
		name: "save frame must be terminated before EOF",
		input: `
data_aBlock
save_aFrame
_someTag someValue`,
		wantSyntaxErr: true,
	}, {
		name: "save frame must be terminated before another data block is opened",
		input: `
data_aBlock
save_aFrame
_someTag someValue
data_anotherBlock
save_`,
		wantSyntaxErr: true,
	}, {
		name: "loops must not be empty",
		input: `
data_aBlock
loop_
`,
		wantSyntaxErr: true,
	}, {
		name: "number of loop values must be a multiple of number of loop tags",
		input: `
data_aBlock
loop_
_tag1 _tag2
val1 val2 val3
`,
		wantSyntaxErr: true,
	}, {
		name: "single tag:value data item saved to appropriate data block",
		input: `
data_aBlock
_someTag someValue`,
		want: CIF{
			DataBlocks: map[string]DataBlock{
				"aBlock": {
					Name: "aBlock",
					DataItems: map[string]any{
						"_someTag": "someValue",
					},
					SaveFrames: map[string]SaveFrame{},
				},
			},
		},
	}, {
		name: "single tag:value data item saved to appropriate save frame",
		input: `
data_aBlock
save_aFrame
_someTag someValue
save_`,
		want: CIF{
			DataBlocks: map[string]DataBlock{
				"aBlock": {
					Name:      "aBlock",
					DataItems: map[string]any{},
					SaveFrames: map[string]SaveFrame{
						"aFrame": {
							Name: "aFrame",
							DataItems: map[string]any{
								"_someTag": "someValue",
							},
						},
					},
				},
			},
		},
	}, {
		name: "loop_ values saved to appropriate data block",
		input: `
data_aBlock
loop_
_tag1
_tag2
val1
2
val3
4.0
"val5"
0.6E1`,
		want: CIF{
			DataBlocks: map[string]DataBlock{
				"aBlock": {
					Name: "aBlock",
					DataItems: map[string]any{
						"_tag1": []any{"val1", "val3", "val5"},
						"_tag2": []any{int64(2), 4.0, 6.0},
					},
					SaveFrames: map[string]SaveFrame{},
				},
			},
		},
	}, {
		name: "loop_ values saved to appropriate save frame",
		input: `
data_aBlock
save_aFrame
loop_
_tag1
_tag2
val1
2
val3
4.0
"val5"
0.6E1
save_`,
		want: CIF{
			DataBlocks: map[string]DataBlock{
				"aBlock": {
					Name:      "aBlock",
					DataItems: map[string]any{},
					SaveFrames: map[string]SaveFrame{
						"aFrame": {
							Name: "aFrame",
							DataItems: map[string]any{
								"_tag1": []any{"val1", "val3", "val5"},
								"_tag2": []any{int64(2), 4.0, 6.0},
							},
						},
					},
				},
			},
		},
	}}

	for _, tc := range testCases {
		t.Run(tc.name, func(t *testing.T) {
			parser := NewParser(strings.NewReader(tc.input))

			got, err := parser.Parse()
			if err, ok := err.(CIFSyntaxError); ok && !tc.wantSyntaxErr {
				t.Fatalf("unexpected syntax error: %v", err)
			} else if !ok && tc.wantSyntaxErr {
				t.Fatalf("no syntax error returned when one was expected")
			} else if ok && tc.wantSyntaxErr {
				return
			}

			if err != nil {
				t.Fatalf("unexpected error: %v", err)
			}

			if diff := cmp.Diff(tc.want, got); diff != "" {
				t.Errorf("incorrect output CIF (-want,+got): %s", diff)
			}
		})
	}

}
