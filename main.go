package main

import (
	"bufio"
	"encoding/json"
	"fmt"
	"io"
	"io/ioutil"
	"os"
	"strings"
)

type Meta struct {
	// shared
	Name        string
	GffVersion  string
	RegionStart int
	RegionEnd   int
	// genbank specific
	Size            int
	Type            string
	GenbankDivision string
	Date            string
	Definition      string
	Accession       string
	Version         string
	Keywords        string
	Organism        string
	Origin          string
	Locus           Locus
	References      []Reference
	Primaries       []Primary
}

type Primary struct {
	RefSeq, PrimaryIdentifier, Primary_Span, Comp string
}

// genbank specific
// type Reference struct {
// 	Authors []string
// 	Title   string
// 	Journal string
// 	PubMed  string
// }

type Reference struct {
	Index, Authors, Title, Journal, PubMed, Remark string
}

type Locus struct {
	Name, SequenceLength, MoleculeType, GenBankDivision, ModDate string
	Circular                                                     bool
}

// from https://github.com/blachlylab/gff3/blob/master/gff3.go
type Feature struct {
	Name string //Seqid in gff, name in gbk
	//gff specific
	Source     string
	Type       string
	Start      int
	End        int
	Score      float64
	Strand     byte
	Phase      int
	Attributes map[string]string // Known as "qualifiers" for gbk, "attributes" for gff.
	//gbk specific
	Location string
	Sequence string
}

type Sequence struct {
	Description string
	Sequence    string
}

type AnnotatedSequence struct {
	Meta     Meta
	Features []Feature
	Sequence Sequence
}

// getReference loops through to get all references within genbank file
func getReference(buf int, lines []string) (Reference, int) {

	var Ref Reference

	Ref.Index = strings.TrimSpace(string(lines[buf])[12:14])
	for {
		if strings.Compare(strings.TrimSpace(string(lines[buf])[0:11]), "AUTHORS") == 0 {
			break
		}
		buf++
	}

	Ref.Authors = strings.TrimSpace(string(lines[buf])[12:])
	buf++
	for {
		if string(lines[buf][0:12]) == "            " {
			Ref.Authors += strings.TrimSpace(string(lines[buf])[12:]) + " "
		} else {
			break
		}
		buf++
	}

	for {
		if strings.Compare(strings.TrimSpace(string(lines[buf])[0:11]), "TITLE") == 0 {
			break
		}
		buf++
	}

	Ref.Title = strings.TrimSpace(string(lines[buf])[12:])
	buf++
	for {
		if string(lines[buf][0:12]) == "            " {
			Ref.Title += strings.TrimSpace(string(lines[buf])[12:]) + " "
		} else {
			break
		}
		buf++
	}

	for {
		if strings.Compare(strings.TrimSpace(string(lines[buf])[0:11]), "JOURNAL") == 0 {
			break
		}
		buf++
	}

	Ref.Journal = strings.TrimSpace(string(lines[buf])[12:])
	buf++
	for {
		if string(lines[buf][0:12]) == "            " {
			Ref.Journal += strings.TrimSpace(string(lines[buf])[12:]) + " "
		} else {
			break
		}
		buf++
	}

	for {
		if strings.Compare(strings.TrimSpace(string(lines[buf])[0:11]), "PUBMED") == 0 {
			break
		}
		buf++
	}

	Ref.PubMed = strings.TrimSpace(string(lines[buf])[12:])

	for {
		if strings.Compare(strings.TrimSpace(string(lines[buf])[0:11]), "REMARK") == 0 {
			break
		} else if strings.Compare(strings.TrimSpace(string(lines[buf])[0:11]), "REFERENCE") == 0 || string(lines[buf][0]) != " " {
			buf--
			return Ref, buf
		}
		buf++
	}

	Ref.Remark = strings.TrimSpace(string(lines[buf])[12:])
	buf++
	for {
		if string(lines[buf][0:12]) == "            " {
			Ref.Remark += strings.TrimSpace(string(lines[buf])[12:]) + " "
		} else {
			break
		}
		buf++
	}

	buf--
	return Ref, buf
}

func parseGbk(path string) {

	/*
	     // If missing filename, print usage
	   	if len(os.Args) < 2 {
	   		fmt.Println("Usage:")
	   		fmt.Println("    go run main.go {location/filename}")
	   		return ""
	   	}
	   	// Open file
	   	filename := os.Args[1]
	*/
	f, err := os.Open(path)
	if err != nil {
		fmt.Println(err)
	}

	bio := bufio.NewReader(f)

	var lines []string
	i := 0

	// Read all lines of the file into buffer
	for {
		line, err := bio.ReadBytes('\n')
		if err == io.EOF {
			break
		}
		if err != nil {
			panic(err)
		}
		sline := strings.TrimRight(string(line), `\n`)
		lines = append(lines, sline)
		i++
	}
	// End read of file into buffer

	// Create meta struct
	meta := Meta{}

	// Create features struct
	// features := []Feature{}

	// Create sequence struct
	// sequence := Sequence{}

	for numLine, line := range lines {
		// fmt.Print / ln(numLine)
		splitLine := strings.Split(line, " ")
		switch splitLine[0] {
		//case should just be alphanumeric character.
		case "":
			continue
		case "LOCUS":
			meta.Locus = parseLocus(line)
		case "DEFINITION":
			baseDefinition := strings.Join(splitLine[1:], " ")
			sublines := lines[numLine+1:]
			for _, subLine := range sublines {
				if string(subLine[0]) == " " {
					baseDefinition = strings.TrimSpace(strings.TrimSpace(baseDefinition) + " " + strings.TrimSpace(subLine))
				} else {
					break
				}
			}
			meta.Definition = baseDefinition
		case "ACCESSION":
			meta.Accession = strings.Join(splitLine[1:], " ")
		case "VERSION":
			meta.Version = strings.Join(splitLine[1:], " ")
		case "KEYWORDS":
			meta.Keywords = strings.Join(splitLine[1:], " ")
		case "SOURCE":
			continue
		case "REFERENCE":
			continue
		case "FEATURES":
			continue
		case "ORIGIN":
			continue
		default:
			fmt.Println(numLine)
			if numLine == 1000000 {
				continue
			}
		}

	}
	file, _ := json.MarshalIndent(meta, "", " ")

	_ = ioutil.WriteFile("test.json", file, 0644)

	//definition parsing

	// // Assign locus to meta struct

	// buf := 0

	// // Search for SOURCE
	// for {
	// 	if strings.Compare(strings.TrimSpace(string(lines[buf])[0:11]), "SOURCE") == 0 {
	// 		break
	// 	}
	// 	buf++
	// }

	// // Populate Source
	// meta.Name = strings.TrimSpace(string(lines[buf])[12:])
	// buf++
	// meta.Organism = strings.TrimSpace(string(lines[buf])[12:]) + " "
	// buf++
	// for {
	// 	if string(lines[buf][0]) == " " {
	// 		meta.Organism += strings.TrimSpace(string(lines[buf])[12:]) + " "
	// 	} else {
	// 		break
	// 	}
	// 	buf++
	// }

	// // Get all References
	// for {
	// 	if strings.Compare(strings.TrimSpace(string(lines[buf])[0:11]), "REFERENCE") == 0 {
	// 		Ref, newBuf := getReference(buf, lines)
	// 		buf = newBuf
	// 		// Append Reference to References array in GenBank struct
	// 		meta.References = append(meta.References, Ref)
	// 	} else {
	// 		break
	// 	}
	// 	// Once search reaches Comment, break from loop
	// 	buf++
	// }

	// buf--

	// // Search for PRIMARY
	// for {
	// 	if len(lines[buf]) >= 11 {
	// 		if strings.Compare(strings.TrimSpace(string(lines[buf])[0:11]), "PRIMARY") == 0 {
	// 			break
	// 		}
	// 	}

	// 	buf++
	// }

	// buf++

	// // Search for Primaries
	// for {
	// 	if string(lines[buf][0:12]) == "            " {
	// 		var P Primary
	// 		P.RefSeq = strings.TrimSpace(string(lines[buf])[12:23])
	// 		P.PrimaryIdentifier = strings.TrimSpace(string(lines[buf])[32:50])
	// 		if len(lines[buf]) > 73 {
	// 			P.Primary_Span = strings.TrimSpace(string(lines[buf])[51:66])
	// 			P.Comp = strings.TrimSpace(string(lines[buf])[72:73])
	// 		} else {
	// 			P.Primary_Span = strings.TrimSpace(string(lines[buf])[51:])
	// 		}
	// 		// Append each Primary to Primaries array in GenBank Struct
	// 		meta.Primaries = append(meta.Primaries, P)
	// 	}
	// 	if strings.Compare(strings.TrimSpace(string(lines[buf])[0:11]), "FEATURES") == 0 {
	// 		break
	// 	}
	// 	buf++
	// }

	// // Search for Features
	// for {
	// 	if strings.Compare(strings.TrimSpace(string(lines[buf])[0:11]), "FEATURES") == 0 {
	// 		feature := Feature{}
	// 		attributes := make(map[string]string)
	// 		if strings.Compare(strings.TrimSpace(string(lines[buf])[0:11]), "FEATURES") == 0 {
	// 			buf++
	// 		}
	// 		for {
	// 			if strings.Compare(strings.TrimSpace(string(lines[buf])[7:8]), "") != 0 {

	// 				// Creating a new feature
	// 				feature.Attributes = make(map[string]string)
	// 				attributes = make(map[string]string)
	// 				feature.Name = strings.TrimSpace(string(lines[buf])[5:21])
	// 				feature.Location = strings.TrimSpace(string(lines[buf])[21:])
	// 				buf++

	// 				for {
	// 					if strings.Compare(strings.TrimSpace(string(lines[buf])[21:22]), "/") == 0 { //parse qualifier
	// 						// Found qualifier
	// 						q := strings.TrimSpace(string(lines[buf])[22:])
	// 						if strings.Contains(q, "=") {
	// 							quarry := strings.Split(q, "=")
	// 							// Add qualifier
	// 							attributes[quarry[0]] = quarry[1]

	// 							// Handle attributes that take up multiple lines
	// 							for {
	// 								if len(lines[buf+1]) > 22 {
	// 									if strings.Compare(strings.TrimSpace(string(lines[buf+1])[21:22]), "/") != 0 && strings.Compare(strings.TrimSpace(string(lines[buf+1])[0:21]), "") == 0 {
	// 										if strings.Compare(strings.TrimSpace(string(lines[buf+1])[0:7]), "ORIGIN") != 0 {
	// 											if quarry[0] == "note" || quarry[0] == "experiment" {
	// 												attributes[quarry[0]] += " " + strings.TrimSpace(string(lines[buf+1])[21:])
	// 											} else {
	// 												attributes[quarry[0]] += strings.TrimSpace(string(lines[buf+1])[21:])
	// 											}
	// 										}
	// 									} else {
	// 										break
	// 									}

	// 								}
	// 								buf++
	// 							}

	// 						}
	// 					} else {
	// 						if attributes != nil {
	// 							feature.Attributes = attributes
	// 							buf--
	// 						}
	// 						break
	// 					}
	// 					buf++
	// 				}
	// 				// Assign each feature to features array in GenBank struct.
	// 				features = append(features, feature)
	// 				if strings.Compare(strings.TrimSpace(string(lines[buf])[0:7]), "ORIGIN") == 0 {
	// 					break
	// 				}

	// 			}
	// 			buf++
	// 		}

	// 	}
	// 	if strings.Compare(strings.TrimSpace(string(lines[buf])[0:7]), "ORIGIN") == 0 {
	// 		break
	// 	}
	// 	buf++
	// }
	// var Origin string
	// // Extract Origin by appending all lines and removing spaces and line information.
	// for {
	// 	if strings.Compare(strings.TrimSpace(string(lines[buf])[0:2]), "//") == 0 {
	// 		break
	// 	} else {
	// 		Origin += strings.TrimSpace(string(lines[buf])[10:])
	// 	}
	// 	buf++
	// }
	// Origin = strings.Replace(Origin, " ", "", -1)

	// // Assign Origin to GenBank struct
	// meta.Origin = Origin

	// // Loop back through all features, grabbing the location information and setting the associated sequence.
	// for i := range features {
	// 	if strings.Contains(features[i].Location, "JOIN") {
	// 		sublocation := strings.Replace(features[i].Location, "JOIN(", "", -1)
	// 		sublocation = strings.Replace(features[i].Location, ")", "", -1)
	// 		sublocations := strings.Split(sublocation, ",")
	// 		Seq := ""
	// 		for location := range sublocations {
	// 			if strings.Contains(string(location), "..") {
	// 				numbers := strings.Split(features[i].Location, "..")
	// 				start, err := strconv.Atoi(numbers[0])
	// 				start--
	// 				end, err := strconv.Atoi(numbers[1])
	// 				if err != nil {
	// 					fmt.Println(err)
	// 				}
	// 				Seq += Origin[start:end]
	// 			} else {
	// 				Seq += string(Origin[location] - 1)
	// 			}
	// 		}
	// 		// Assign associated sequence to GenBank's feature if a joined value
	// 		features[i].Sequence = Seq
	// 	} else if strings.Contains(features[i].Location, "..") {
	// 		numbers := strings.Split(features[i].Location, "..")
	// 		start, err := strconv.Atoi(numbers[0])
	// 		start--
	// 		end, err := strconv.Atoi(numbers[1])
	// 		// Assign associated sequence to GenBank's feature if a range
	// 		features[i].Sequence = Origin[start:end]
	// 		if err != nil {
	// 			fmt.Println(err)
	// 		}
	// 	} else {
	// 		k, err := strconv.Atoi(features[i].Location)
	// 		k--
	// 		if err != nil {
	// 			fmt.Println(err)
	// 		}
	// 		if k >= 0 {
	// 			// Assign associated sequence to GenBank's feature if a single value
	// 			features[i].Sequence = string(Origin[k])
	// 		}

	// 	}
	// }

	// annotatedSequence := AnnotatedSequence{Meta: meta, Features: features}

	// // convert GenBank struct into JSON (indented) format
	// out, err := json.MarshalIndent(annotatedSequence, "", "    ")
	// if err != nil {
	// 	panic(err)
	// }
	// fmt.Println(out)

}

func main() {

	// fmt.Println(parseGff("data/ecoli-mg1655.gff"))
	// parseGbk("data/addgene-plasmid-50005-sequence-74677.gbk")
	parseGbk("data/test.gbk")
}
