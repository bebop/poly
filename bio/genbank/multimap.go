/*
This file provides utilities for working with a MultiMap,
which is simply a map which can store multiple values for a single key instead
of the usual one.

Useful for when we expect to encounter repeated keys but we want to store all pairs,
not just the last-inserted one. This implementation has the advantage of being compatible with
json.Marshal, cmp.Diff, pretty printing, and bracket indexing out of the box.
Does not make uniqueness quarantees for key value pairs.
Currently only used in genbank.Feature.Attributes, however his may end up
being useful for other parsers which allow for repeated keys, in which case
this should be made into its own module.
*/
package genbank

// MultiMap is defined as a simple type alias over a map of slices.
type MultiMap[K, V comparable] map[K][]V

// NewMultiMap creates a new empty multimap.
func NewMultiMap[K, V comparable]() MultiMap[K, V] {
	return make(map[K][]V)
}

// Put adds a key-value pair to the multimap.
func Put[K, V comparable](m MultiMap[K, V], k K, v ...V) {
	if _, ok := m[k]; !ok {
		m[k] = v
	} else {
		m[k] = append(m[k], v...)
	}
}

// ForEachKey iterates over the multimap, once for each key with all values passed as a slice.
// do is a callback that takes the key, values slice for that key
// This exists purely as a convenience function, if your use case would benefit from
// early break/return, it is recommended you do the usual range iteration operator.
func ForEachKey[K, V comparable](m MultiMap[K, V], do func(K, []V)) {
	for k, values := range m {
		do(k, values)
	}
}

// ForEachValue iterates over the multimap, once for each value
// do is a callback that takes and a key and value.
func ForEachValue[K, V comparable](m MultiMap[K, V], do func(K, V)) {
	ForEachKey(m, func(k K, values []V) {
		for _, v := range values {
			do(k, v)
		}
	})
}

// MapSlice efficiently applies a transformation to each element of a slice to create a new slice
func MapSlice[X any, Y any](slice []X, mapper func(X) Y) []Y {
	y := make([]Y, len(slice))
	for i, x := range slice {
		y[i] = mapper(x)
	}
	return y
}

// identity returns its input, useful for using MapSlice to do a shallow copy
func identity[X any](x X) X {
	return x
}
