/*
Package mash is for sketching sequence data to make it easier to search for and against.

The package is named mash after the mash sketching algorithm, which is based on the MinHash algorithm.

Mash: fast genome and metagenome distance estimation using MinHash.
Ondov, B.D., Treangen, T.J., Melsted, P. et al.
Genome Biol 17, 132 (2016).
https://doi.org/10.1186/s13059-016-0997-x

Mash Screen: high-throughput sequence containment estimation for genome discovery.
Ondov, B., Starrett, G., Sappington, A. et al.
Genome Biol 20, 232 (2019).
https://doi.org/10.1186/s13059-019-1841-x

The idea is that we can sketch a sequence of data, and then compare the sketch to other sketches to see how similar they are.
This saves a bunch of computation time and memory, because we don't have to compare the entire sequence to another sequence.

TTFN,
Tim
*/
package mash

import "github.com/spaolacci/murmur3" // murmur3 is a fast non-cryptographic hash algorithm that was also used in the original papers-> https://github.com/shenwei356/go-hashing-kmer-bench

// Mash is a collection of hashes of kmers from a given sequence.
type Mash struct {
	KmerSize   uint     // The kmer size is the size of the sliding window that is used to generate the hashes.
	SketchSize uint     // The sketch size is the number of hashes to store.
	Sketches   []uint32 // The sketches are the hashes of the kmers that we can compare to other sketches.
}

func NewMash(kmerSize uint, sketchSize uint) *Mash {
	return &Mash{
		KmerSize:   kmerSize,
		SketchSize: sketchSize,
		Sketches:   make([]uint32, sketchSize),
	}
}

func (m *Mash) Sketch(sequence string) {
	// slide a window of size k along the sequence
	for kmerStart := 0; kmerStart < len(sequence)-int(m.KmerSize); kmerStart++ {
		kmer := sequence[kmerStart : kmerStart+int(m.KmerSize)]
		// hash the kmer to a 32 bit number
		hash := murmur3.Sum32([]byte(kmer))
		// keep the minimum hash value of all the kmers in the window up to a given sketch size
		// the sketch is a vector of the minimum hash values
		var biggestHashIndex int

		// find the biggest hash value in the sketch
		for i := 0; i < len(m.Sketches); i++ {
			if m.Sketches[i] == 0 {
				biggestHashIndex = i
				break
			} else if m.Sketches[i] > m.Sketches[biggestHashIndex] {
				biggestHashIndex = i
			}
		}
		m.Sketches[biggestHashIndex] = hash
	}
}

func (m *Mash) Distance(other *Mash) float64 {
	var sameHashes int
	for i := 0; i < len(m.Sketches); i++ {
		for j := 0; j < len(other.Sketches); j++ {
			if m.Sketches[i] == other.Sketches[j] {
				sameHashes++
				break
			}
		}
	}
	return 1 - (float64(sameHashes) / float64(len(m.Sketches)))
}
