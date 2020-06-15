package main

import "testing"

func TestHello(t *testing.T) {
	want := "4f0fc3582aeb0667d268bab48a8237036cc7334abb2269b42c6fa25fa4ce9938"
	if got := seqhash("Hello World", true); got != want {
		t.Errorf("Seqhash(\"Hello World\", true) = %q, want %q", got, want)
	}
}

