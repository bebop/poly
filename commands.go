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
	"strings"
	"sync"

	"github.com/urfave/cli/v2"
)

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

func convert(c *cli.Context) error {
	if isPipe() {

		annotatedSequence := parseStdin(c)

		var output []byte

		// logic for chosing output format, then builds string to be output.
		if c.String("o") == "json" {
			output, _ = json.MarshalIndent(annotatedSequence, "", " ")
		} else if c.String("o") == "gff" {
			output = BuildGff(annotatedSequence)
		} else if c.String("o") == "gbk" || c.String("o") == "gb" {
			output = BuildGbk(annotatedSequence)
		}

		// output to stdout
		fmt.Print(string(output))

		//
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
					WriteJSON(annotatedSequence, outputPath+".json")
				} else if c.String("o") == "gff" {
					WriteGff(annotatedSequence, outputPath+".gff")
				} else if c.String("o") == "gbk" || c.String("o") == "gb" {
					WriteGbk(annotatedSequence, outputPath+".gbk")
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
func hash(c *cli.Context) {

	if isPipe() {
		annotatedSequence := parseStdin(c)
		sequenceHash := flagSwitchHash(c, annotatedSequence)
		if c.String("o") == "json" {
			annotatedSequence.Sequence.Hash = sequenceHash
			annotatedSequence.Sequence.HashFunction = strings.ToUpper(c.String("f"))
			output, _ := json.MarshalIndent(annotatedSequence, "", " ")
			fmt.Print(string(output))
		} else {
			if strings.ToUpper(c.String("f")) == "NO" {
				fmt.Print(sequenceHash)
			} else {
				fmt.Println(sequenceHash)
			}
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
				if c.String("o") == "json" {
					annotatedSequence.Sequence.Hash = sequenceHash
					annotatedSequence.Sequence.HashFunction = strings.ToUpper(c.String("f"))

					if c.Bool("stdout") == true {
						output, _ := json.MarshalIndent(annotatedSequence, "", " ")
						fmt.Print(string(output))
					} else {
						outputPath := match[0 : len(match)-len(extension)]
						// should have way to support wildcard matches for varied output names.
						WriteJSON(annotatedSequence, outputPath+".json")
					}

				} else {
					// should refactor before push
					if len(matches) == 1 {
						if strings.ToUpper(c.String("f")) == "NO" {
							fmt.Print(sequenceHash)
						} else {
							fmt.Println(sequenceHash)
						}
					} else {
						if strings.ToUpper(c.String("f")) == "NO" {
							fmt.Print(sequenceHash)
						} else {
							fmt.Println(sequenceHash, " ", match)
						}
					}
				}

				// decrementing wait group.
				wg.Done()

			}(match) // passing match to Go routine anonymous function.

		}

		// waiting outside for loop for Go routines so they can run concurrently.
		wg.Wait()

	}

}

// a simple helper function to convert an *os.File type into a string.
func stdinToString(file *os.File) string {
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
func isPipe() bool {
	info, _ := os.Stdin.Stat()
	flag := false
	if info.Mode()&os.ModeNamedPipe != 0 {
		// we have a pipe input
		flag = true
	}
	return flag
}

// a simple helper function to take stdin from a pipe and parse it into an annotated sequence
func parseStdin(c *cli.Context) AnnotatedSequence {
	var annotatedSequence AnnotatedSequence

	// logic for determining input format, then parses accordingly.
	if c.String("i") == "json" {
		json.Unmarshal([]byte(stdinToString(os.Stdin)), &annotatedSequence)
	} else if c.String("i") == "gbk" || c.String("i") == "gb" {
		annotatedSequence = ParseGbk(stdinToString(os.Stdin))
	} else if c.String("i") == "gff" {
		annotatedSequence = ParseGff(stdinToString(os.Stdin))
	}
	return annotatedSequence
}

// helper function to hash sequence based on flag using generic hash.
func flagSwitchHash(c *cli.Context, annotatedSequence AnnotatedSequence) string {

	var hashString string
	switch strings.ToUpper(c.String("f")) {
	case "MD5":
		hashString = annotatedSequence.hash(crypto.MD5)
	case "SHA1":
		hashString = annotatedSequence.hash(crypto.SHA1)
	case "SHA244":
		hashString = annotatedSequence.hash(crypto.SHA224)
	case "SHA256":
		hashString = annotatedSequence.hash(crypto.SHA256)
	case "SHA384":
		hashString = annotatedSequence.hash(crypto.SHA384)
	case "SHA512":
		hashString = annotatedSequence.hash(crypto.SHA512)
	case "RIPEMD160":
		hashString = annotatedSequence.hash(crypto.RIPEMD160)
	case "SHA3_224":
		hashString = annotatedSequence.hash(crypto.SHA3_224)
	case "SHA3_256":
		hashString = annotatedSequence.hash(crypto.SHA3_256)
	case "SHA3_384":
		hashString = annotatedSequence.hash(crypto.SHA3_384)
	case "SHA3_512":
		hashString = annotatedSequence.hash(crypto.SHA3_512)
	case "SHA512_224":
		hashString = annotatedSequence.hash(crypto.SHA512_224)
	case "SHA512_256":
		hashString = annotatedSequence.hash(crypto.SHA512_256)
	case "BLAKE2s_256":
		hashString = annotatedSequence.hash(crypto.BLAKE2s_256)
	case "BLAKE2b_256":
		hashString = annotatedSequence.hash(crypto.BLAKE2b_256)
	case "BLAKE2b_384":
		hashString = annotatedSequence.hash(crypto.BLAKE2b_384)
	case "BLAKE2b_512":
		hashString = annotatedSequence.hash(crypto.BLAKE2b_512)
	case "BLAKE3":
		hashString = annotatedSequence.blake3Hash()
	case "NO":
		hashString = RotateSequence(annotatedSequence.Sequence.Sequence)
	default:
		hashString = annotatedSequence.blake3Hash()
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
func fileParser(c *cli.Context, match string) AnnotatedSequence {
	extension := filepath.Ext(match)
	var annotatedSequence AnnotatedSequence

	// determining which reader to use and parse into AnnotatedSequence struct.
	if extension == ".gff" || c.String("i") == "gff" {
		annotatedSequence = ReadGff(match)
	} else if extension == ".gbk" || extension == ".gb" || c.String("i") == "gbk" || c.String("i") == "gb" {
		annotatedSequence = ReadGbk(match)
	} else if extension == ".json" || c.String("i") == "json" {
		annotatedSequence = ReadJSON(match)
	} else {
		// TODO put default error handling here.
	}
	return annotatedSequence
}
