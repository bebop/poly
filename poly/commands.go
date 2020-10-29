package main

import (
	"bufio"
	"bytes"
	"crypto"
	"encoding/json"
	"fmt"
	"io"
	"os"
	"path/filepath"
	"strconv"
	"strings"
	"sync"

	"github.com/TimothyStiles/poly"
	"github.com/urfave/cli/v2"
	"lukechampine.com/blake3"
)

/******************************************************************************
Oct, 15, 2020

File is structured as so:

	Top level commands:
		Convert
		Hash
		Translate
		Optimize

	Helper functions


This file contains a majority of the code that runs when command line routines
are run. The one exception is that argument flags and helper text for each
command are defined in main.go which then makes calls to their corresponding
function in this file. This keeps main.go clean and readable.

To ease development you'll notice there's a common command line template in
the section below along with helper functions for parsing stdin, getting blob
matches, etc which almost every command ends up using.

Each command must also have a test in commands_test.go that demonstrates its
correct usage by executing a bash subprocess. You can read more about that at
the top of commands_test.go itself.

Happy hacking,
Tim

https://github.com/urfave/cli/issues/731
******************************************************************************/

/******************************************************************************

convert currently has two modes. Pipe and fileio.

The function isPipe() detects if input is coming from a pipe like:

	cat data/bsub.gbk | poly c -i gbk -o json > test.json

In this case the output goes directly to standard out and can be redirected
into a file.

If not from a pipe convert checks args for file patterns to find, then iterates
over each matched file pattern to read in a file, then spit out the desired
output.

For example:

	poly c -o json *.gbk *.gff

will read all files in a directory with ".gbk" or ".gff" as their extension
parse them, and then spit out a similiarly named file with the .json extension.

******************************************************************************/

func convertCommand(c *cli.Context) error {
	if isPipe(c) {

		annotatedSequence := parseStdin(c)

		var output []byte

		// logic for chosing output format, then builds string to be output.
		if c.String("o") == "json" {
			output, _ = json.MarshalIndent(annotatedSequence, "", " ")
		} else if c.String("o") == "gff" {
			output = poly.BuildGff(annotatedSequence)
		} else if c.String("o") == "gbk" || c.String("o") == "gb" {
			output = poly.BuildGbk(annotatedSequence)
		}

		// output to stdout
		fmt.Fprint(c.App.Writer, string(output))

	} else {

		// gets glob pattern matches to determine which files to use.
		matches := getMatches(c)

		// TODO write basic check to see if input flag or all paths have accepted file extensions.

		// TODO write basic check for reduduncy. I.E converting gff to gff, etc.

		// declaring wait group outside loop
		var wg sync.WaitGroup

		// concurrently iterate through each pattern match, read the file, output to new format.
		for _, match := range matches {

			// incrementing wait group for Go routine
			wg.Add(1)

			// executing Go routine.
			go func(match string) {
				extension := filepath.Ext(match)
				annotatedSequence := fileParser(c, match)

				// determining output format and name, then writing out to name.
				outputPath := match[0 : len(match)-len(extension)]
				if c.String("o") == "json" {
					poly.WriteJSON(annotatedSequence, outputPath+".json")
				} else if c.String("o") == "gff" {
					poly.WriteGff(annotatedSequence, outputPath+".gff")
				} else if c.String("o") == "gbk" || c.String("o") == "gb" {
					poly.WriteGbk(annotatedSequence, outputPath+".gbk")
				}

				// decrementing wait group.
				wg.Done()

			}(match) // passing match to Go routine anonymous function.

			// TODO add delete flag with confirmation.

		}

		// waiting outside for loop for Go routines so they can run concurrently.
		wg.Wait()
	}

	return nil
}

/******************************************************************************

hash currently has two modes. Pipe and fileio.

The function isPipe() detects if input is coming from a pipe like:

	cat data/bsub.gbk | poly hash -i gbk -o json > test.json

In this case the output goes directly to standard out and can be redirected
into a file.

Without the "-o json" json flag only the hash is printed to stdout.

to force all output to go to stdout use --stdout.

If not from a pipe convert checks args for file patterns to find, then iterates
over each matched file pattern to read in a file, then spit out the desired
output.

For example:

	poly hash -o json *.gbk *.gff

will read all files in a directory with ".gbk" or ".gff" as their extension
parse them, and then spit out a similiarly named file with the .json extension along with their hashes.

******************************************************************************/
func hashCommand(c *cli.Context) error {

	if isPipe(c) {
		annotatedSequence := parseStdin(c)                   // get sequence from stdin
		sequenceHash := flagSwitchHash(c, annotatedSequence) // get hash include no-op which only rotates the sequence

		// handler for outputting String to stdout <- Default for pipes
		if c.String("o") == "string" {

			if strings.ToUpper(c.String("f")) == "NO" {
				fmt.Fprint(c.App.Writer, sequenceHash) // prints with no newline
			} else {
				fmt.Fprintln(c.App.Writer, sequenceHash) // prints with newline
			}

		}

		// handler for outputting JSON to stdout
		if c.String("o") == "json" {
			annotatedSequence.Sequence.SequenceHash = sequenceHash                           // adding hash to JSON
			annotatedSequence.Sequence.SequenceHashFunction = strings.ToUpper(c.String("f")) // adding hash type to JSON
			output, _ := json.MarshalIndent(annotatedSequence, "", " ")
			fmt.Fprint(c.App.Writer, string(output))
		}

	} else {

		// gets glob pattern matches to determine which files to use.
		matches := getMatches(c)

		// declaring wait group outside loop
		var wg sync.WaitGroup

		// concurrently iterate through each pattern match, read the file, output to new format.
		for _, match := range matches {

			// incrementing wait group for Go routine
			wg.Add(1)

			// executing Go routine.
			go func(match string) {
				extension := filepath.Ext(match)
				annotatedSequence := fileParser(c, match)
				sequenceHash := flagSwitchHash(c, annotatedSequence)

				// handler for outputting String <- Default
				if strings.ToLower(c.String("o")) == "string" {

					// if there's only one match then
					if len(matches) == 1 {

						if strings.ToUpper(c.String("f")) == "NO" {
							fmt.Fprint(c.App.Writer, sequenceHash)
						} else {
							// fmt.Println(sequenceHash)
							fmt.Fprintln(c.App.Writer, sequenceHash)
						}

					} else { // if there's more than one match then print each hash

						if strings.ToUpper(c.String("f")) == "NO" {
							fmt.Fprintln(c.App.Writer, sequenceHash) // just prints list of rotated, unhashed sequences
						} else {
							fmt.Fprintln(c.App.Writer, sequenceHash, " ", match) // prints list of hashes and corresponding filepaths.
						}

					}
				}

				// handler for outputting JSON.
				if strings.ToLower(c.String("o")) == "json" {
					annotatedSequence.Sequence.SequenceHash = sequenceHash
					annotatedSequence.Sequence.SequenceHashFunction = strings.ToUpper(c.String("f"))

					if c.Bool("--log") == true {
						output, _ := json.MarshalIndent(annotatedSequence, "", " ")
						fmt.Fprint(c.App.Writer, string(output))
					} else {
						outputPath := match[0 : len(match)-len(extension)]
						// should have way to support wildcard matches for varied output names.
						poly.WriteJSON(annotatedSequence, outputPath+".json")
					}

				}

				// decrementing wait group.
				wg.Done()

			}(match) // passing match to Go routine anonymous function.

		}

		// waiting outside for loop for Go routines so they can run concurrently.
		wg.Wait()

	}
	return nil
}

/******************************************************************************

poly optimize currently has one mode. Pipe io.

The function isPipe() detects if input is coming from a pipe like:

	cat data/DNAsequence.txt | poly optimize > test.txt

	poly optimize can also have it's codon table specified
	and files can be provided to add weight to the codon table

	cat data/DNAsequence.txt | poly optimize -ct 11 -wt data/puc19.gbk > test.txt

In the future there will be multi file support. Check our issues on github to see what's up.

******************************************************************************/
func optimizeCommand(c *cli.Context) error {

	// get appropriate codon table from flag
	var codonTable poly.CodonTable

	if isNumeric(c.String("ct")) {
		codonTableNumber, _ := strconv.Atoi(c.String("ct"))
		codonTable = poly.DefaultCodonTablesByNumber[codonTableNumber]
	} else {
		codonTable = poly.DefaultCodonTablesByName[c.String("ct")]
	}

	// if a file exists to weigh the table. Weigh it.
	if fileExists(c.String("wt")) {
		targetOrganism := fileParser(c, c.String("wt"))
		codonTable.CreateWeights(targetOrganism.Sequence.Sequence)
	}

	if isPipe(c) {

		// uncomment below to parse annotatedSequence from pipe
		annotatedSequence := parseStdin(c)
		var aminoAcids string

		if c.Bool("aa") {
			aminoAcids = annotatedSequence.Sequence.Sequence
		} else {
			aminoAcids = poly.Translate(annotatedSequence.Sequence.Sequence, codonTable)
		}

		optimizedSequence := poly.Optimize(aminoAcids, codonTable)
		fmt.Fprintln(c.App.Writer, optimizedSequence)

	}
	return nil
}

/******************************************************************************

poly translate currently has one mode. Pipe io.

The function isPipe() detects if input is coming from a pipe like:

	cat data/DNAsequence.txt | poly translate > test.txt

	poly optimize can also have it's codon table specified

	cat data/DNAsequence.txt | poly translate --codon-table Standard > test.txt

In the future there will be multi file support. Check our issues on github to see what's up.

******************************************************************************/

func translateCommand(c *cli.Context) error {

	if isPipe(c) {

		// get appropriate codon table from flag
		var codonTable poly.CodonTable

		if isNumeric(c.String("ct")) {
			codonTableNumber, _ := strconv.Atoi(c.String("ct"))
			codonTable = poly.DefaultCodonTablesByNumber[codonTableNumber]
		} else {
			codonTable = poly.DefaultCodonTablesByName[c.String("ct")]
		}

		// uncomment below to parse annotatedSequence from pipe
		annotatedSequence := parseStdin(c)

		aminoAcids := poly.Translate(annotatedSequence.Sequence.Sequence, codonTable)

		fmt.Fprintln(c.App.Writer, aminoAcids)

	}
	return nil
}

// a simple helper function to convert an *os.File type into a string.
func stdinToString(file io.Reader) string {
	var stringBuffer bytes.Buffer
	reader := bufio.NewReader(file)
	for {
		input, _, err := reader.ReadRune()
		if err != nil && err == io.EOF {
			break
		}
		stringBuffer.WriteRune(input)
	}
	return stringBuffer.String()
}

// a simple helper function to remove duplicates from a list of strings.
// Used to reduce reduncy in filepath pattern matching.
// from https://gist.github.com/johnwesonga/6301924
func uniqueNonEmptyElementsOf(s []string) []string {
	unique := make(map[string]bool, len(s))
	us := make([]string, len(unique))
	for _, elem := range s {
		if len(elem) != 0 {
			if !unique[elem] {
				us = append(us, elem)
				unique[elem] = true
			}
		}
	}

	return us

}

// a simple helper function to determine if there is input coming from stdin pipe.
func isPipe(c *cli.Context) bool {
	info, _ := os.Stdin.Stat()
	flag := false
	if info.Mode()&os.ModeNamedPipe != 0 {
		// we have a pipe input
		flag = true
	}
	if c.App.Reader != os.Stdin {
		flag = true
	}
	return flag
}

func isNumeric(s string) bool {
	_, err := strconv.ParseFloat(s, 64)
	return err == nil
}

// a simple helper function to take stdin from a pipe and parse it into an annotated sequence
func parseStdin(c *cli.Context) poly.AnnotatedSequence {
	var annotatedSequence poly.AnnotatedSequence
	// logic for determining input format, then parses accordingly.
	if c.String("i") == "json" {
		json.Unmarshal([]byte(stdinToString(c.App.Reader)), &annotatedSequence)
	} else if c.String("i") == "gbk" || c.String("i") == "gb" {
		annotatedSequence = poly.ParseGbk(stdinToString(c.App.Reader))
	} else if c.String("i") == "gff" {
		annotatedSequence = poly.ParseGff(stdinToString(c.App.Reader))
	} else if c.String("i") == "string" {
		annotatedSequence.Sequence.Sequence = stdinToString(c.App.Reader)
	}
	return annotatedSequence
}

// helper function to hash sequence based on flag using generic hash.
func flagSwitchHash(c *cli.Context, annotatedSequence poly.AnnotatedSequence) string {

	var hashString string
	switch strings.ToUpper(c.String("f")) {
	case "MD5":
		hashString = annotatedSequence.Hash(crypto.MD5.New())
	case "SHA1":
		hashString = annotatedSequence.Hash(crypto.SHA1.New())
	case "SHA244":
		hashString = annotatedSequence.Hash(crypto.SHA224.New())
	case "SHA256":
		hashString = annotatedSequence.Hash(crypto.SHA256.New())
	case "SHA384":
		hashString = annotatedSequence.Hash(crypto.SHA384.New())
	case "SHA512":
		hashString = annotatedSequence.Hash(crypto.SHA512.New())
	case "RIPEMD160":
		hashString = annotatedSequence.Hash(crypto.RIPEMD160.New())
	case "SHA3_224":
		hashString = annotatedSequence.Hash(crypto.SHA3_224.New())
	case "SHA3_256":
		hashString = annotatedSequence.Hash(crypto.SHA3_256.New())
	case "SHA3_384":
		hashString = annotatedSequence.Hash(crypto.SHA3_384.New())
	case "SHA3_512":
		hashString = annotatedSequence.Hash(crypto.SHA3_512.New())
	case "SHA512_224":
		hashString = annotatedSequence.Hash(crypto.SHA512_224.New())
	case "SHA512_256":
		hashString = annotatedSequence.Hash(crypto.SHA512_256.New())
	case "BLAKE2s_256":
		hashString = annotatedSequence.Hash(crypto.BLAKE2s_256.New())
	case "BLAKE2b_256":
		hashString = annotatedSequence.Hash(crypto.BLAKE2b_256.New())
	case "BLAKE2b_384":
		hashString = annotatedSequence.Hash(crypto.BLAKE2b_384.New())
	case "BLAKE2b_512":
		hashString = annotatedSequence.Hash(crypto.BLAKE2b_512.New())
	case "BLAKE3":
		hashString = annotatedSequence.Hash(blake3.New(32, nil))
		// hashString = annotatedSequence.Blake3Hash()
	case "NO":
		hashString = poly.RotateSequence(annotatedSequence.Sequence.Sequence)
	default:
		hashString = annotatedSequence.Hash(blake3.New(32, nil))
		break
	}
	return hashString
}

// helper function to get unique glob patterns from cli.context
func getMatches(c *cli.Context) []string {
	var matches []string

	//take all args and get their pattern matches.
	for argIndex := 0; argIndex < c.Args().Len(); argIndex++ {
		match, _ := filepath.Glob(c.Args().Get(argIndex))
		matches = append(matches, match...)
	}

	//filtering pattern matches for duplicates.
	matches = uniqueNonEmptyElementsOf(matches)

	return matches

}

// function to parse whatever file is at a matched path.
func fileParser(c *cli.Context, match string) poly.AnnotatedSequence {
	extension := filepath.Ext(match)
	var annotatedSequence poly.AnnotatedSequence

	// determining which reader to use and parse into AnnotatedSequence struct.
	if extension == ".gff" || c.String("i") == "gff" {
		annotatedSequence = poly.ReadGff(match)
	} else if extension == ".gbk" || extension == ".gb" || c.String("i") == "gbk" || c.String("i") == "gb" {
		annotatedSequence = poly.ReadGbk(match)
	} else if extension == ".json" || c.String("i") == "json" {
		annotatedSequence = poly.ReadJSON(match)
	} else {
		// TODO put default error handling here.
	}
	return annotatedSequence
}

// fileExists checks if a file exists and is not a directory before we
// try using it to prevent further errors.
// from https://golangcode.com/check-if-a-file-exists/

func fileExists(filename string) bool {
	info, err := os.Stat(filename)
	if os.IsNotExist(err) {
		return false
	}
	return !info.IsDir()
}
