/*
This file provides utilities for working with a MultiMap,
which is simply a map which can store multiple values for a single key instead
of the usual one.

Useful for when we expect to encounter repeated keys but we want to keep all pairs,
not just the latest one, while preserving O(1) time lookup cost.
Does not make uniqueness quarantees for key value pairs.
This may end up being useful for other parsers which allow for repeated keys, in which case
this should be made into its own module.
*/
package genbank

// defined as a simple type alias over a map of slices
// while not ideal (eg computing total number of items takes O(N))
// this has the advantage of being compatible with json.Marshal, cmp.Diff,
// pretty printing, and bracket indexing out of the box.
type MultiMap[K, V comparable] map[K][]V

// create a new empty multimap
func NewMultiMap[K, V comparable]() MultiMap[K, V] {
	return make(map[K][]V)
}

// adds a key-value pair to the multimap
func Put[K, V comparable](m MultiMap[K, V], k K, v ...V) {
	if _, ok := m[k]; !ok {
		m[k] = v
	} else {
		m[k] = append(m[k], v...)
	}
}

// iterates over the multimap, once for each key with all values passed as a slice
func ForEachKey[K, V comparable](m MultiMap[K, V], do func(K, []V)) {
	for k, values := range m {
		do(k, values)
	}
}

// iterates over the multimap, once for each value
func ForEachValue[K, V comparable](m MultiMap[K, V], do func(K, V)) {
	ForEachKey(m, func(k K, values []V) {
		for _, v := range values {
			do(k, v)
		}
	})
}

// efficiently apply a transformation to each element of a slice to create a new slice
func MapSlice[X any, Y any](slice []X, mapper func(X) Y) []Y {
	y := make([]Y, len(slice))
	for i, x := range slice {
		y[i] = mapper(x)
	}
	return y
}

// the identity function which returns its input, useful for using MapSlice to do a shallow copy
func identity[X any](x X) X {
	return x
}
