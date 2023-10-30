package genbank

// defines a new MultiMap type which can store multiple values for a single key
// useful for when we expect repeated keys, while preserving O(1) lookup
// does not make uniqueness guarantees for values (you can repeat key-value pairs)
type MultiMap[K, V comparable] map[K][]V

// create a new empty multimap
func NewMultiMap[K, V comparable]() MultiMap[K, V] {
	return make(map[K][]V)
}

// adds a key-value pair to the multimap
func Put[K, V comparable](m MultiMap[K, V], k K, v V) {
	if _, ok := m[k]; !ok {
		m[k] = []V{v}
	} else {
		m[k] = append(m[k], v)
	}
}

// iterates over the multimap, once for each key
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

// returns number of unique keys
func KeyCount[K, V comparable](m MultiMap[K, V]) int {
	return len(m)
}

// efficiently each element of a slice to create a new slice
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
