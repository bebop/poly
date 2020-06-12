package main

import (
	"bufio"
	"bytes"
	"encoding/json"
	"fmt"
	"io"
	"os"
	"path/filepath"
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

		var annotatedSequence AnnotatedSequence

		// logic for determining input format, then parses accordingly.
		if c.String("i") == "json" {
			json.Unmarshal([]byte(stdinToString(os.Stdin)), &annotatedSequence)
		} else if c.String("i") == "gbk" || c.String("i") == "gb" {
			annotatedSequence = ParseGbk(stdinToString(os.Stdin))
		} else if c.String("i") == "gff" {
			annotatedSequence = ParseGff(stdinToString(os.Stdin))
		}

		var output []byte

		// logic for chosing output format, then builds string to be output.
		if c.String("o") == "json" {
			output, _ = json.MarshalIndent(annotatedSequence, "", " ")
		} else if c.String("o") == "gff" {
			output = BuildGff(annotatedSequence)
		}

		// output to stdout
		fmt.Print(string(output))

		//
	} else {

		var matches []string

		//take all args and get their pattern matches.
		for argIndex := 0; argIndex < c.Args().Len(); argIndex++ {
			match, _ := filepath.Glob(c.Args().Get(argIndex))
			matches = append(matches, match...)
		}

		//filtering pattern matches for duplicates.
		matches = uniqueNonEmptyElementsOf(matches)

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

				// determining output format and name, then writing out to name.
				outputPath := match[0 : len(match)-len(extension)]
				if c.String("o") == "json" {
					WriteJSON(annotatedSequence, outputPath+".json")
				} else if c.String("o") == "gff" {
					WriteGff(annotatedSequence, outputPath+".gff")
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
