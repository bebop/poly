package seqhash

// import "testing"

// func TestLSH_Add(t *testing.T) {
// 	lsh := NewLSH(10, 2)

// 	// Add an item with 2 features
// 	item1 := &Item{
// 		ID:       1,
// 		Features: []string{"apple", "orange"},
// 	}
// 	lsh.Add(item1)

// 	// Check that the item was added to the map of items
// 	if lsh.items[1] != item1 {
// 		t.Errorf("Expected item 1 to be in map of items, got %v", lsh.items[1])
// 	}

// 	// Check that the item was added to the correct bucket
// 	if len(lsh.buckets[0]) != 1 || lsh.buckets[0][0] != item1 {
// 		t.Errorf("Expected item 1 to be in bucket 0, got %v", lsh.buckets[0])
// 	}

// }

// func TestLSH_Search(t *testing.T) {
// 	lsh := NewLSH(10, 2)

// 	// Add some items to LSH
// 	item1 := &Item{
// 		ID:       1,
// 		Features: []string{"apple", "orange"},
// 	}
// 	lsh.Add(item1)
// 	item2 := &Item{
// 		ID:       2,
// 		Features: []string{"apple", "banana"},
// 	}
// 	lsh.Add(item2)
// 	item3 := &Item{
// 		ID:       3,
// 		Features: []string{"grape", "strawberry"},
// 	}
// 	lsh.Add(item3)

// 	// Search for similar items to item 1
// 	similar := lsh.Search(item1)
// 	if len(similar) != 1 || similar[0] != item2 {
// 		t.Errorf("Expected item 2 to be similar to item 1, got %v", similar)
// 	}

// 	// Search for similar items to item 3
// 	similar = lsh.Search(item3)
// 	if len(similar) != 0 {
// 		t.Errorf("Expected no items to be similar to item 3, got %v", similar)
// 	}
// }
