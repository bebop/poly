// TODO: create benchmark for space and time costs vs https://github.com/shenwei356/bio/tree/master/sketches
// TODO: create efficient storage for children + counts in top node via multidimensional array
// TODO: integrate with IO parsers
// TODO: enable reading and writing to binary format
// TODO: iterate over observed kmers (what is the appropriate iteration method?)
// TODO: comparison between two KmerCounter's
// TODO: initialize upper levels while respecting cache locality
// TODO: add counts from other KmerCounter, to enable reducing parallel counts
// TODO: lookup neighbors
// TODO: path compression for internal nodes
// TODO: have last succ point to next node from same level for faster iteration

package alphabet

import (
	"fmt"
)

type KmerCounter struct {
	alphabet    *Alphabet
	num_symbols uint8
	max_k       uint8
	total       uint64
	children    []Node
}

type Node struct {
	succ     *Node
	child    *Node
	encoding uint8
	count    uint32
}

func NewKmerCounter(a *Alphabet, max_k uint8) *KmerCounter {
	kc := new(KmerCounter)
	kc.alphabet = a
	kc.num_symbols = uint8(len(a.symbols))
	kc.max_k = max_k

	kc.children = make([]Node, kc.num_symbols)
	for i := uint8(1); i < kc.num_symbols; i++ {
		kc.children[i].encoding = i
		kc.children[i-1].succ = &kc.children[i]
	}
	return kc
}

func lookupChild(n *Node, encoding uint8) *Node {
	for c := n.child; c != nil && encoding <= c.encoding; c = c.succ {
		if encoding == c.encoding {
			return c
		}
	}
	return nil
}

func insertChild(p *Node, encoding uint8, n *Node) {
	if p.child == nil {
		p.child = n
	} else if n.encoding < p.child.encoding {
		n.succ = p.child
		p.child = n
	} else {
		for c := p.child; c.encoding < encoding; c = c.succ {
			if encoding < c.succ.encoding {
				n.succ = c.succ
				c.succ = n
			}
		}
	}
	n.count++ // not sure why this is needed
	n.encoding = encoding
}

func Observe(kc *KmerCounter, seq string) error {
	cbuf, err := kc.alphabet.EncodeAll(seq[:kc.max_k])
	if err != nil {
		return err
	}

	CreateChildren := func(index int, remaining int) *Node {
		if remaining == 0 {
			return nil
		} // this condition should never happen
		nodes := make([]Node, remaining)

		for j, n := range nodes {
			n.count = 1
			if j != 0 {
				nodes[j-1].child = &n
				n.encoding = cbuf[(index+j)%int(kc.max_k)]
			}
		}

		return &nodes[0]
	}

	UpdateBuffer := func(index int) (max_k int, err error) {
		max_k = int(kc.max_k)
		lookahead := index + max_k - 1
		var encoding uint8
		if lookahead < len(seq) {
			next := string(seq[lookahead])
			encoding, err = kc.alphabet.Encode(next)
			if err != nil {
				err = fmt.Errorf("in position %d: %w", index, err)
				return
			}
			cbuf[lookahead%int(kc.max_k)] = encoding
		} else {
			max_k = len(seq) - index
		}
		return
	}

	var encoding uint8
	for i := range seq {
		max_k, err := UpdateBuffer(i)
		if err != nil {
			return err
		}

		p := &kc.children[cbuf[i%int(kc.max_k)]]
		kc.total++
		for k := 1; k <= max_k; k++ {
			p.count++
			// fmt.Printf("%d %d %v\n", i, k, p)

			if k != max_k {
				encoding = cbuf[(i+k)%int(kc.max_k)]
				c := lookupChild(p, encoding)
				if c == nil {
					insertChild(p, encoding, CreateChildren(i+k, max_k-k))
					break // inserted nodes already have count = 1 added
				}
				p = c
			}
		}
	}
	return nil
}

func LookupCount(kc *KmerCounter, kmer string) (count uint32, err error) {
	if len(kmer) > int(kc.max_k) {
		err = fmt.Errorf("kmer_counter: attempted to lookup count of %d-mer which exceeds max supported length %d", len(kmer), kc.max_k)
		return
	}
	if len(kmer) == 0 {
		err = fmt.Errorf("kmer_counter: attempted to lookup count of 0mer")
		return
	}
	encoded, err := kc.alphabet.EncodeAll(kmer)
	if err != nil {
		return
	}

	var k int
	var p *Node
	for k, p = 1, &kc.children[encoded[0]]; k < len(kmer) && p != nil; k++ {
		p = lookupChild(p, encoded[k])
	}
	if p == nil {
		return
	}
	count = p.count
	return
}

func LookupFrequency(kc *KmerCounter, kmer string) (float64, error) {
	count, err := LookupCount(kc, kmer)
	if err != nil {
		return 0, err
	}
	return float64(count) / float64(kc.total-uint64(len(kmer)-1)), nil
}
