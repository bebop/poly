package seqhash

import (
	_ "embed"
)

// All the modifications were pulled
//
// #/bin/bash
// curl -s https://www.idtdna.com/site/catalog/Modifications/GetAllMods | grep '</a>' | grep 'Modifications' | sed 's/.*>\(.*\)<.*/\1/' > modifications.txt
//
// #/bin/python
// import json
// f = open("modifications.txt", "r")
// currentDict = "five_prime"
// d={"five_prime":{},"internal":{},"three_prime":{}}
// i = 0
// for line in f:
//     line = line.strip()
//     d[currentDict][line] = i
//     i+=1
//     if line == "5SUN":
//         i = 0
//         currentDict = "internal"
//         d[currentDict]["DNA"] = i
//         i+=1
//         d[currentDict]["RNA"] = i
//         i+=1
//     if line == "i2MOErT":
//         i = 0
//         currentDict = "three_prime"
// with open("modifications_v2.json", "w") as outfile:
//     outfile.write(json.dumps(d))
// #/bin/bash
// rm modifications.txt

//go:embed modifications_v2.json
var modificationsJsonString string

// Rules:
// 1. Circular sequence cannot have 5' or 3' modifications
// 2. 5' and 3' modifications ADD to the last base pair, but do not affect indexing
// 3. DNA = 0 and RNA = 1. If either occur in a DNA or RNA sequence, they are removed
// 4. All internally methylated sequences are converted first to x
// 5. Indexing begins at zero. This is unlike many conventions in the DNA world - watch out
// 6. No modifications are allowed for proteins
