package seqhash

import (
	"math/rand"
)

// Item represents an item in a dataset
type Item struct {
	ID       int      // ID of the item
	Features []string // Features of the item
}

// LSH represents a locality-sensitive hashing algorithm
type LSH struct {
	numBuckets int           // Number of buckets to use for hashing
	hashFuncs  []hashFunc    // Hash functions to use for each bucket
	buckets    [][]*Item     // Buckets to store hashed items
	items      map[int]*Item // Map of items by ID
}

// hashFunc represents a hash function used by LSH
type hashFunc func(string) int

// NewLSH creates a new instance of LSH
func NewLSH(numBuckets int, numHashFuncs int) *LSH {
	lsh := &LSH{
		numBuckets: numBuckets,
		hashFuncs:  make([]hashFunc, numHashFuncs),
		buckets:    make([][]*Item, numBuckets),
		items:      make(map[int]*Item),
	}

	// Generate random hash functions
	for i := 0; i < numHashFuncs; i++ {
		lsh.hashFuncs[i] = func(s string) int {
			return rand.Intn(lsh.numBuckets)
		}
	}

	return lsh
}

// Add adds an item to LSH
func (lsh *LSH) Add(item *Item) {
	// Hash the item's features and add it to the corresponding bucket
	for _, hashFunc := range lsh.hashFuncs {
		for _, feature := range item.Features {
			bucketID := hashFunc(feature)
			lsh.buckets[bucketID] = append(lsh.buckets[bucketID], item)
		}
	}

	// Add the item to the map of items
	lsh.items[item.ID] = item
}

// Search searches for items similar to the given item
func (lsh *LSH) Search(item *Item) []*Item {
	// Hash the item's features and search for similar items in the corresponding buckets
	var similar []*Item
	for _, hashFunc := range lsh.hashFuncs {
		for _, feature := range item.Features {
			bucketID := hashFunc(feature)
			for _, i := range lsh.buckets[bucketID] {
				if i.ID != item.ID {
					similar = append(similar, i)
				}
			}
		}
	}

	return similar
}
