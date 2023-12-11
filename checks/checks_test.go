package checks_test

import (
	"testing"

	"github.com/bebop/poly/checks"
)

// This also needs an example test.
func TestIsPalindromic(t *testing.T) {
	ecori := checks.IsPalindromic("GAATTC")
	if ecori != true {
		t.Errorf("IsPalindromic failed to call EcoRI a palindrome")
	}
	bsai := checks.IsPalindromic("GGTCTC")
	if bsai != false {
		t.Errorf("IsPalindromic failed call BsaI NOT a palindrome")
	}
}

func TestGcContent(t *testing.T) {
	content := checks.GcContent("GGTATC")
	if content != 0.5 {
		t.Errorf("GcContent did not properly calculate GC content")
	}
}

func TestIsDNA(t *testing.T) {
	tests := []struct {
		name string
		args string
		want bool
	}{
		{
			name: "Success",
			args: "GATTACA",
			want: true,
		},
		{
			name: "FailRNA",
			args: "GAUUACA",
			want: false,
		},
		{
			name: "FailUnknown",
			args: "RANDOM STRING",
			want: false,
		},
	}
	for _, tt := range tests {
		t.Run(tt.name, func(t *testing.T) {
			if got := checks.IsDNA(tt.args); got != tt.want {
				t.Errorf("IsDNA() = %v, want %v", got, tt.want)
			}
		})
	}
}

func TestIsRNA(t *testing.T) {
	tests := []struct {
		name string
		args string
		want bool
	}{
		{
			name: "Success",
			args: "GAUUACA",
			want: true,
		},
		{
			name: "FailDNA",
			args: "GATTACA",
			want: false,
		},
		{
			name: "FailUnknown",
			args: "RANDOM STRING",
			want: false,
		},
	}
	for _, tt := range tests {
		t.Run(tt.name, func(t *testing.T) {
			if got := checks.IsRNA(tt.args); got != tt.want {
				t.Errorf("IsRNA() = %v, want %v", got, tt.want)
			}
		})
	}
}
