package io

import "io"

// trieNode is a simple implementation of a trie with values
// of a single type. Only supports []byte keys for simplicity.
type trieNode[V any] struct {
	hasValue bool
	value    V
	children map[byte]*trieNode[V]
}

func (tn *trieNode[V]) Set(key []byte, val V) {
	// Written iteratively rather than recursively for performance
	// reasons and to avoid OOM errors.
	currTn := tn
	var nextTn *trieNode[V]
	exists := false

	// Iterate down to the requested node, creating nodes
	// along the way if necessary.
	for i := 0; i < len(key); i++ {
		nextTn, exists = currTn.children[key[i]]
		if !exists {
			nextTn = newTrieNode[V]()
			currTn.children[key[i]] = nextTn
		}
		currTn = nextTn
	}

	currTn.hasValue = true
	currTn.value = val
}

func newTrieNode[V any]() *trieNode[V] {
	return &trieNode[V]{
		hasValue: false,
		children: make(map[byte]*trieNode[V]),
	}
}

// A ReplacingReader wraps an io.Reader, performing the requested
// replacements to its input.
//
// Matches the strings to be replaced in a greedy manner.
type ReplacingReader struct {
	innerR  io.Reader
	buf     []byte
	atEOF   bool
	replace *trieNode[[]byte]
}

// Wrap r, replacing all occurences of the keys in replace
// with their associated values.
func NewReplacingReader(r io.Reader, replace map[string]string) io.Reader {
	res := &ReplacingReader{
		innerR: r,
	}

	// Build the replacement trie.
	res.replace = newTrieNode[[]byte]()
	for key, val := range replace {
		res.replace.Set([]byte(key), []byte(val))
	}

	return res
}

// Wrap r, normalizing all newlines ("\n", "\r\n", or "\r") to '\n'.
func NewNewlineNormalizingReader(r io.Reader) io.Reader {
	return NewReplacingReader(r, map[string]string{
		"\n":   "\n",
		"\r":   "\n",
		"\r\n": "\n",
	})
}

// Reads data into p, (greedily) performing the requested replacements.
func (rr *ReplacingReader) Read(p []byte) (n int, err error) {
	// If we have any data in the buffer left over from the last call,
	// read it into p and return.
	if len(rr.buf) != 0 {
		n = copy(p, rr.buf)

		// Clear the buffer if we were able to read everything.
		// Otherwise, advance the buffer.
		if n < len(rr.buf) {
			rr.buf = rr.buf[n:]
		} else {
			rr.buf = nil
		}

		return
	}

	// If we've encountered an io.EOF, return an io.EOF.
	if rr.atEOF {
		return 0, io.EOF
	}

	// Read from the underlying reader.
	innerP := make([]byte, len(p))
	innerN, innerErr := rr.innerR.Read(innerP)
	if innerErr == io.EOF {
		rr.atEOF = true
	} else if innerErr != nil {
		return 0, innerErr
	}

	// Copy data into p making the necessary replacements.
	n = 0
	i := 0
	for n < len(p) && i < innerN {
		// Match the longest possible replacement.
		currTn := rr.replace
		exists := true
		var lastReplacement []byte
		matchSize := 0
		for j := 0; exists; j++ {
			if currTn.hasValue {
				lastReplacement = currTn.value
				matchSize = j
			}
			currTn, exists = currTn.children[innerP[i+j]]
		}

		// If we found a match, copy the replacement data into p
		// and skip what was matched. Otherwise, copy a single byte.
		if matchSize > 0 {
			for j, b := range lastReplacement {
				// If we run out of space in p to put the replacement into,
				// store the remainder in the buffer and return.
				if n == len(p) {
					rr.buf = lastReplacement[j:]
					return
				}

				p[n] = b
				n++
			}
			i = i + matchSize
		} else {
			p[n] = innerP[i]
			n++
			i++
		}
	}
	return
}
