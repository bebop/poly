package io_test

import (
	"bytes"
	"io"
	"testing"

	pio "github.com/TimothyStiles/poly/io"
	"github.com/google/go-cmp/cmp"
)

func TestReplacingReader(t *testing.T) {
	testCases := []struct {
		name    string
		input   string
		replace map[string]string
		want    string
	}{{
		name:    "does nothing to an empty string",
		input:   "",
		replace: map[string]string{"hi": "hey"},
		want:    "",
	}, {
		name:    "performs a single replacement",
		input:   "replace me please",
		replace: map[string]string{"replace me": "replaced"},
		want:    "replaced please",
	}, {
		name:    "does nothing when no replacement provided",
		input:   "do not replace me please",
		replace: map[string]string{},
		want:    "do not replace me please",
	}, {
		name:    "makes multiple replacements",
		input:   "replace me please, replace me again",
		replace: map[string]string{"replace me": "replaced"},
		want:    "replaced please, replaced again",
	}, {
		name:    "matches greedily",
		input:   "replace me, no actually replace me instead",
		replace: map[string]string{"replace me": "short", "replace me instead": "long"},
		want:    "short, no actually long",
	}}

	for _, tc := range testCases {
		t.Run(tc.name, func(t *testing.T) {
			r := bytes.NewReader([]byte(tc.input))
			rr := pio.NewReplacingReader(r, tc.replace)

			res, err := io.ReadAll(rr)
			if err != nil {
				t.Fatalf("unexpected error: %v", err)
			}

			got := string(res)
			if diff := cmp.Diff(tc.want, got); diff != "" {
				t.Errorf("unexpected output (-want,+got): %v", diff)
			}
		})
	}
}

func TestNewlineNormalizingReader(t *testing.T) {
	testCases := []struct {
		name  string
		input string
		want  string
	}{{
		name:  "does nothing to an empty string",
		input: "",
		want:  "",
	}, {
		name:  "conserves \n newlines",
		input: "hi\ntest\nagain\n",
		want:  "hi\ntest\nagain\n",
	}, {
		name:  "normalizes \r newlines",
		input: "hi\ntest\ragain\n",
		want:  "hi\ntest\nagain\n",
	}, {
		name:  "normalizes \r\n newlines",
		input: "hi\ntest\r\nagain\n",
		want:  "hi\ntest\nagain\n",
	}}

	for _, tc := range testCases {
		t.Run(tc.name, func(t *testing.T) {
			r := bytes.NewReader([]byte(tc.input))
			rr := pio.NewNewlineNormalizingReader(r)

			res, err := io.ReadAll(rr)
			if err != nil {
				t.Fatalf("unexpected error: %v", err)
			}

			got := string(res)
			if diff := cmp.Diff(tc.want, got); diff != "" {
				t.Errorf("unexpected output (-want,+got): %v", diff)
			}
		})
	}

}
