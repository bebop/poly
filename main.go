package main

import (
	"bufio"
	"bytes"
	"fmt"
	"io"
	"io/ioutil"
	"os"
	"strconv"
	"strings"
)

var genbankDivisions = []string{
	"PRI", //primate sequences
	"ROD", //rodent sequences
	"MAM", //other mamallian sequences
	"VRT", //other vertebrate sequences
	"INV", //invertebrate sequences
	"PLN", //plant, fungal, and algal sequences
	"BCT", //bacterial sequences
	"VRL", //viral sequences
	"PHG", //bacteriophage sequences
	"SYN", //synthetic sequences
	"UNA", //unannotated sequences
	"EST", //EST sequences (expressed sequence tags)
	"PAT", //patent sequences
	"STS", //STS sequences (sequence tagged sites)
	"GSS", //GSS sequences (genome survey sequences)
	"HTG", //HTG sequences (high-throughput genomic sequences)
	"HTC", //unfinished high-throughput cDNA sequencing
	"ENV", //environmental sampling sequences
}

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

func parseGff(path string) AnnotatedSequence {
	file, _ := ioutil.ReadFile(path)
	splitFile := strings.Split(string(file), "\n")
	metaString := splitFile[0:2]
	versionString := metaString[0]
	regionStringArray := strings.Split(metaString[1], " ")

	meta := Meta{}
	meta.GffVersion = strings.Split(versionString, " ")[1]
	meta.Name = regionStringArray[1] // Formally region name, but changed to name here for generality/interoperability.
	meta.RegionStart, _ = strconv.Atoi(regionStringArray[2])
	meta.RegionEnd, _ = strconv.Atoi(regionStringArray[3])
	meta.Size = meta.RegionEnd - meta.RegionStart

	records := []Feature{}
	sequence := Sequence{}
	var sequenceBuffer bytes.Buffer
	fastaFlag := false
	for _, line := range splitFile {
		if line == "##FASTA" {
			fastaFlag = true
		} else if len(line) == 0 {
			continue
		} else if line[0:2] == "##" {
			continue
		} else if fastaFlag == true && line[0:1] != ">" {
			// sequence.Sequence = sequence.Sequence + line
			sequenceBuffer.WriteString(line)
		} else if fastaFlag == true && line[0:1] == ">" {
			sequence.Description = line
		} else {
			record := Feature{}
			fields := strings.Split(line, "\t")
			record.Name = fields[0]
			record.Source = fields[1]
			record.Type = fields[2]
			record.Start, _ = strconv.Atoi(fields[3])
			record.End, _ = strconv.Atoi(fields[4])
			record.Score, _ = strconv.ParseFloat(fields[5], 64)
			record.Strand = fields[6][0]
			record.Phase, _ = strconv.Atoi(fields[7])
			record.Attributes = make(map[string]string)
			attributes := fields[8]
			var eqIndex int
			for i := strings.Index(attributes, ";"); i > 0; i = strings.Index(attributes, ";") {
				eqIndex = strings.Index(attributes[:i], "=")
				record.Attributes[attributes[:i][:eqIndex]] = attributes[:i][eqIndex+1:]
				attributes = attributes[i+1:]
			}
			records = append(records, record)
		}
	}
	sequence.Sequence = sequenceBuffer.String()
	annotatedSequence := AnnotatedSequence{}
	annotatedSequence.Meta = meta
	annotatedSequence.Features = records
	annotatedSequence.Sequence = sequence

	return annotatedSequence
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
	// meta := Meta{}

	// Create features struct
	// features := []Feature{}

	// Create sequence struct
	// sequence := Sequence{}

	// Populate Locus
	//locus parser
	locus := Locus{}

	locusString := strings.TrimSpace(lines[0])
	locusSplit := strings.Split(locusString, " ")
	var filteredLocusSplit []string
	for i := range locusSplit {
		if locusSplit[i] != "" {
			filteredLocusSplit = append(filteredLocusSplit, locusSplit[i])
		}
	}
	locus.Name = filteredLocusSplit[1]
	locus.SequenceLength = strings.Join([]string{filteredLocusSplit[2], filteredLocusSplit[3]}, " ")
	locus.MoleculeType = filteredLocusSplit[4]

	// locus.Name = strings.TrimSpace(string(lines[0])[12:25])
	// locus.SequenceLength = strings.TrimSpace(string(lines[0])[36:41])
	// locus.MoleculeType = strings.TrimSpace(string(lines[0])[41:43])
	// locus.GenBankDivision = strings.TrimSpace(string(lines[0])[47:51])
	// locus.ModDate = strings.TrimSpace(string(lines[0])[68:79])
	fmt.Println(locus)

	// // Assign locus to meta struct
	// meta.Locus = locus

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
	parseGbk("data/addgene-plasmid-50005-sequence-74677.gbk")
}
