package poly

import (
	"bytes"
	"strings"
	"sync"
)

// complementBaseRuneMap provides 1:1 mapping between bases and their complements
var complementBaseRuneMap = map[rune]rune{
	65:  84,  // A -> T
	66:  86,  // B -> V
	67:  71,  // C -> G
	68:  72,  // D -> H
	71:  67,  // G -> C
	72:  68,  // H -> D
	75:  77,  // K -> M
	77:  75,  // M -> K
	78:  78,  // N -> N
	82:  89,  // R -> Y
	83:  83,  // S -> S
	84:  65,  // T -> A
	85:  65,  // U -> A
	86:  66,  // V -> B
	87:  87,  // W -> W
	89:  82,  // Y -> R
	97:  116, // a -> t
	98:  118, // b -> v
	99:  103, // c -> g
	100: 104, // d -> h
	103: 99,  // g -> a
	104: 100, // h -> d
	107: 109, // k -> m
	109: 107, // m -> k
	110: 110, // n -> n
	114: 121, // r -> y
	115: 115, // s -> s
	116: 97,  // t -> a
	117: 97,  // u -> a
	118: 98,  // v -> b
	119: 119, // w -> w
	121: 114, // y -> r
}

// GetSequence is a method to get the full sequence of an annotated sequence
func (sequence Sequence) GetSequence() string {
	return sequence.Sequence
}

// GetSequence is a method wrapper to get a Feature's sequence. Mutates with Sequence.
func (feature Feature) GetSequence() string {
	return getFeatureSequence(feature, feature.SequenceLocation)
}

// ReverseComplement takes the reverse complement of a sequence
func ReverseComplement(sequence string) string {
	complementString := strings.Map(ComplementBase, sequence)
	n := len(complementString)
	newString := make([]rune, n)
	for _, base := range complementString {
		n--
		newString[n] = base
	}
	return string(newString)
}

// ComplementBase accepts a base pair and returns its complement base pair
func ComplementBase(basePair rune) rune {
	return complementBaseRuneMap[basePair]
}

// getFeatureSequence takes a feature and location object and returns a sequence string.
func getFeatureSequence(feature Feature, location Location) string {
	var sequenceBuffer bytes.Buffer
	var sequenceString string
	parentSequence := feature.ParentSequence.Sequence

	if len(location.SubLocations) == 0 {
		sequenceBuffer.WriteString(parentSequence[location.Start:location.End])
	} else {

		for _, subLocation := range location.SubLocations {
			sequenceBuffer.WriteString(getFeatureSequence(feature, subLocation))
		}
	}

	// reverse complements resulting string if needed.
	if location.Complement {
		sequenceString = ReverseComplement(sequenceBuffer.String())
	} else {
		sequenceString = sequenceBuffer.String()
	}

	return sequenceString
}

// the following functions Iter and iter, derive from github.com/schwarmco/go-cartesian-product
// which uses interfaces, so I modified it to use runes
func Iter(params ...[]rune) chan []rune {
	// create channel
	c := make(chan []rune)
	// create waitgroup
	var wg sync.WaitGroup
	// call iterator
	wg.Add(1)
	iterate(&wg, c, []rune{}, params...)
	// call channel-closing go-func
	go func() { wg.Wait(); close(c) }()
	// return channel
	return c
}

// private, recursive Iteration-Function
func iterate(wg *sync.WaitGroup, channel chan []rune, result []rune, params ...[]rune) {
	// dec WaitGroup when finished
	defer wg.Done()
	// no more params left?
	if len(params) == 0 {
		// send result to channel
		channel <- result
		return
	}
	// shift first param
	p, params := params[0], params[1:]
	// iterate over it
	for i := 0; i < len(p); i++ {
		// inc WaitGroup
		wg.Add(1)
		// create copy of result
		resultCopy := append([]rune{}, result...)
		// call self with remaining params
		go iterate(wg, channel, append(resultCopy, p[i]), params...)
	}
}

func allVariantsIUPAC(seq string) []string {
	var allVariants = []string{}
	var iupacList = [][]rune{}
	iupac := map[rune][]rune{
		'G': []rune{'G'},
		'A': []rune{'A'},
		'T': []rune{'T'},
		'C': []rune{'C'},
		'R': []rune{'G', 'A'},
		'Y': []rune{'T', 'C'},
		'M': []rune{'A', 'C'},
		'K': []rune{'G', 'T'},
		'S': []rune{'G', 'C'},
		'W': []rune{'A', 'T'},
		'H': []rune{'A', 'C', 'T'},
		'B': []rune{'G', 'T', 'C'},
		'V': []rune{'G', 'C', 'A'},
		'D': []rune{'G', 'A', 'T'},
		'N': []rune{'G', 'A', 'T', 'C'},
	}
	for _, s := range seq {
		iupacList = append(iupacList, iupac[s])
	}

	cartesianProducts := Iter(iupacList...)
	for product := range cartesianProducts {
		allVariants = append(allVariants, string(product))
	}

	return allVariants

}
