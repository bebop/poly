package main

import (
	"bufio"
	"bytes"
	"crypto"
	"encoding/json"
	"fmt"
	"io"
	"log"
	"os"
	"path/filepath"
	"regexp"
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
				sequence := fileParser(c, match)
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
		sequence := parseStdin(c)                   // get sequence from stdin
		sequenceHash := flagSwitchHash(c, sequence) // get hash include no-op which only rotates the sequence
		printHash(c, sequenceHash, "-")

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
				sequence := fileParser(c, match)
				sequenceHash := flagSwitchHash(c, sequence)
				printHash(c, sequenceHash, match)

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
		codonTable = poly.GetCodonTable(codonTableNumber)
	}

	// if a file exists to weigh the table. Weigh it.
	if fileExists(c.String("wt")) {
		targetOrganism := fileParser(c, c.String("wt"))
		codonTable.OptimizeTable(targetOrganism.Sequence)
	}

	if isPipe(c) {

		// uncomment below to parse sequence from pipe
		sequence := parseStdin(c)
		var aminoAcids string

		if c.Bool("aa") {
			aminoAcids = sequence.Sequence
		} else {
			aminoAcids = poly.Translate(sequence.Sequence, codonTable)
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
			codonTable = poly.GetCodonTable(codonTableNumber)
		}

		sequence := parseStdin(c)

		aminoAcids := poly.Translate(sequence.Sequence, codonTable)

		fmt.Fprintln(c.App.Writer, aminoAcids)

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

func isNumeric(s string) bool {
	_, err := strconv.ParseFloat(s, 64)
	return err == nil
}

// a simple helper function to take stdin from a pipe and parse it into an Sequence
func parseStdin(c *cli.Context) poly.Sequence {
	var sequence poly.Sequence
	// logic for determining input format, then parses accordingly.
	if c.String("i") == "json" {
		err := json.Unmarshal([]byte(stdinToBytes(c.App.Reader)), &sequence)
		if err != nil {
			// do proper error handling here.
		}
	} else if c.String("i") == "gbk" || c.String("i") == "gb" {
		sequence = poly.ParseGbk(stdinToBytes(c.App.Reader))
	} else if c.String("i") == "gff" {
		sequence = poly.ParseGff(stdinToBytes(c.App.Reader))
	} else if c.String("i") == "fasta" {
		sequence = poly.ParseFASTA(stdinToBytes(c.App.Reader))
	} else if c.String("i") == "string" {
		sequence.Sequence = parseText(stdinToBytes(c.App.Reader))
	}
	return sequence
}

// a simple helper function to take stdin from a pipe and parse it into an Sequence
func parseFlag(file []byte, flag string) poly.Sequence {
	var sequence poly.Sequence
	// logic for determining input format, then parses accordingly.
	if flag == "json" {
		err := json.Unmarshal(file, &sequence)
		if err != nil {
			// todo proper error testing here.
		}
	} else if flag == "gbk" || flag == "gb" {
		sequence = poly.ParseGbk(file)
	} else if flag == "gff" {
		sequence = poly.ParseGff(file)
	} else if flag == "fasta" {
		sequence = poly.ParseFASTA(file)
	} else if flag == "string" {
		sequence.Sequence = string(file)
	}
	return sequence
}

// helper function to hash sequence based on flag using generic hash.
func flagSwitchHash(c *cli.Context, sequence poly.Sequence) string {

	var hashString string
	hashFunctionName := strings.ToUpper(c.String("f"))
	switch hashFunctionName {
	case "MD5":
		hashString = sequence.Hash(crypto.MD5.New())
	case "SHA1":
		hashString = sequence.Hash(crypto.SHA1.New())
	case "SHA244":
		hashString = sequence.Hash(crypto.SHA224.New())
	case "SHA256":
		hashString = sequence.Hash(crypto.SHA256.New())
	case "SHA384":
		hashString = sequence.Hash(crypto.SHA384.New())
	case "SHA512":
		hashString = sequence.Hash(crypto.SHA512.New())
	case "RIPEMD160":
		hashString = sequence.Hash(crypto.RIPEMD160.New())
	case "SHA3_224":
		hashString = sequence.Hash(crypto.SHA3_224.New())
	case "SHA3_256":
		hashString = sequence.Hash(crypto.SHA3_256.New())
	case "SHA3_384":
		hashString = sequence.Hash(crypto.SHA3_384.New())
	case "SHA3_512":
		hashString = sequence.Hash(crypto.SHA3_512.New())
	case "SHA512_224":
		hashString = sequence.Hash(crypto.SHA512_224.New())
	case "SHA512_256":
		hashString = sequence.Hash(crypto.SHA512_256.New())
	case "BLAKE2S_256":
		hashString = sequence.Hash(crypto.BLAKE2s_256.New())
	case "BLAKE2B_256":
		hashString = sequence.Hash(crypto.BLAKE2b_256.New())
	case "BLAKE2B_384":
		hashString = sequence.Hash(crypto.BLAKE2b_384.New())
	case "BLAKE2B_512":
		hashString = sequence.Hash(crypto.BLAKE2b_512.New())
	case "BLAKE3":
		hashString = sequence.Hash(blake3.New(32, nil))
	case "NO":
		hashString = poly.RotateSequence(sequence.Sequence)
	default:
		hashString = sequence.Hash(blake3.New(32, nil))
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

// function to parse whatever file is at a matched path.
func fileParser(c *cli.Context, match string) poly.Sequence {
	extension := filepath.Ext(match)
	var sequence poly.Sequence

	// determining which reader to use and parse into Sequence struct.
	if extension == ".gff" || c.String("i") == "gff" {
		sequence = poly.ReadGff(match)
	} else if extension == ".gbk" || extension == ".gb" || c.String("i") == "gbk" || c.String("i") == "gb" {
		sequence = poly.ReadGbk(match)
	} else if extension == ".json" || c.String("i") == "json" {
		sequence = poly.ReadJSON(match)
	} else if extension == ".fasta" || c.String("i") == "fasta" {
		sequence = poly.ReadFASTA(match)
	}
	// TODO put default error handling here.
	return sequence
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
	} else if c.String("o") == "string" || c.String("o") == "txt" {
		output = []byte(sequence.Sequence)
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

func parseText(file []byte) string {
	fileString := string(file)
	reg, err := regexp.Compile("\\*[^a-zA-Z]+")
	if err != nil {
		log.Fatal(err)
	}
	processedString := reg.ReplaceAllString(fileString, "")
	return processedString
}
