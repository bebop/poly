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

	"github.com/TimothyStiles/poly"
	"github.com/urfave/cli/v2"
)

/******************************************************************************
Oct, 15, 2020

File is structured as so:

	Top level commands:
		Convert
		Hash

These commands are mainly for developer utility.
They are well tested but really, really shouldn't be used in production.
Unless there's a strong case for a new command line utility tool it's very unlikely that I'll
merge a new command so please ask before developing one with the intent of merging a pull request.
Too much effort for not enough pay off just to maintain it in most cases.

TTFN,
Tim

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

		sequence := parseStdin(c)

		output := buildStdOut(c, sequence)

		// logic for chosing output format, then builds string to be output.

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
				sequence := parseExt(match)
				writeFile(c, sequence, match)
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
		sequence := parseStdin(c) // get sequence from stdin
		hash, _ := sequence.Hash()
		printHash(c, hash, "-")

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
				sequence := parseExt(match)
				hash, _ := sequence.Hash()
				printHash(c, hash, match)

				// decrementing wait group.
				wg.Done()

			}(match) // passing match to Go routine anonymous function.

		}

		// waiting outside for loop for Go routines so they can run concurrently.
		wg.Wait()

	}
	return nil
}

// a simple helper function to convert an *os.File type into a string.
func stdinToBytes(file io.Reader) []byte {
	var stringBuffer bytes.Buffer
	reader := bufio.NewReader(file)
	for {
		input, _, err := reader.ReadRune()
		if err != nil && err == io.EOF {
			break
		}
		stringBuffer.WriteRune(input)
	}
	return stringBuffer.Bytes()
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

// a simple helper function to take stdin from a pipe and parse it into an Sequence
func parseStdin(c *cli.Context) poly.Sequence {
	var sequence poly.Sequence
	// logic for determining input format, then parses accordingly.
	if c.String("i") == "json" {
		_ = json.Unmarshal([]byte(stdinToBytes(c.App.Reader)), &sequence)
	} else if c.String("i") == "gbk" || c.String("i") == "gb" {
		sequence = poly.ParseGbk(stdinToBytes(c.App.Reader))
	} else if c.String("i") == "gff" {
		sequence = poly.ParseGff(stdinToBytes(c.App.Reader))
	} else if c.String("i") == "fasta" {
		sequence = poly.ParseFASTA(stdinToBytes(c.App.Reader))
	}
	return sequence
}

// a simple helper function to take stdin from a pipe and parse it into an Sequence
func parseFlag(file []byte, flag string) poly.Sequence {
	var sequence poly.Sequence
	// logic for determining input format, then parses accordingly.
	if flag == "json" {
		_ = json.Unmarshal(file, &sequence)
	} else if flag == "gbk" || flag == "gb" {
		sequence = poly.ParseGbk(file)
	} else if flag == "gff" {
		sequence = poly.ParseGff(file)
	} else if flag == "fasta" {
		sequence = poly.ParseFASTA(file)
	}
	return sequence
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

func buildStdOut(c *cli.Context, sequence poly.Sequence) []byte {
	var output []byte
	if c.String("o") == "json" {
		output, _ = json.MarshalIndent(sequence, "", " ")
	} else if c.String("o") == "gff" {
		output = poly.BuildGff(sequence)
	} else if c.String("o") == "gbk" || c.String("o") == "gb" {
		output = poly.BuildGbk(sequence)
	} else if c.String("o") == "fasta" {
		output = poly.BuildFASTA(sequence)
	}
	return output
}

func writeFile(c *cli.Context, sequence poly.Sequence, match string) {
	// determining output format and name, then writing out to name.
	inputExtension := filepath.Ext(match)

	var outputPath string
	var outputExtension string
	if filepath.Ext(c.String("o")) == "" {
		outputExtension = c.String("o")
		outputPath = match[0:len(match)-len(inputExtension)] + "." + outputExtension
	} else {
		outputPath = c.String("o")
		outputExtension = filepath.Ext(c.String("o"))[1:]
	}

	if match == outputPath {
		fmt.Fprintln(c.App.Writer, "WARNING: "+"input path and output path match. File: "+match+". Skipping to prevent possible data loss. Try providing a full path with a different name using the -o flag.")
	} else if outputExtension == "json" {
		poly.WriteJSON(sequence, outputPath)
	} else if outputExtension == "gff" {
		poly.WriteGff(sequence, outputPath)
	} else if outputExtension == "gbk" || c.String("o") == "gb" {
		poly.WriteGbk(sequence, outputPath)
	} else if outputExtension == "fasta" {
		poly.WriteFASTA(sequence, outputPath)
	}

}

func printHash(c *cli.Context, hash string, path string) {
	output := formatHashOutput(hash, path)
	fmt.Fprint(c.App.Writer, output)
}

func formatHashOutput(hash string, path string) string {
	return hash + "  " + path + "\n"
}

func parseExt(match string) poly.Sequence {
	extension := filepath.Ext(match)
	var sequence poly.Sequence

	// determining which reader to use and parse into Sequence struct.
	if extension == ".gff" {
		sequence = poly.ReadGff(match)
	} else if extension == ".gbk" || extension == ".gb" {
		sequence = poly.ReadGbk(match)
	} else if extension == ".json" {
		sequence = poly.ReadJSON(match)
	} else if extension == ".fasta" {
		sequence = poly.ReadFASTA(match)
	}
	return sequence
}
