package main

import (
	"bufio"
	"bytes"
	"encoding/json"
	"fmt"
	"io"
	"os"
	"path/filepath"

	"github.com/urfave/cli/v2"
)

func convert(c *cli.Context) error {
	if isPipe() {

		var annotatedSequence AnnotatedSequence

		// logic for determining input format.
		if c.String("i") == "json" {
			json.Unmarshal([]byte(stdinToString(os.Stdin)), &annotatedSequence)
		} else if c.String("i") == "gbk" || c.String("i") == "gb" {
			annotatedSequence = ParseGbk(stdinToString(os.Stdin))
		} else if c.String("i") == "gff" {
			annotatedSequence = ParseGff(stdinToString(os.Stdin))
		}

		var output []byte

		// logic for chosing output format.
		if c.String("o") == "json" {
			output, _ = json.MarshalIndent(annotatedSequence, "", " ")
		} else if c.String("o") == "gff" {
			output = BuildGff(annotatedSequence)
		}

		// output to stdout
		fmt.Print(string(output))

	} else {

		var matches []string

		//take all args except for last and get their pattern matches.
		for argIndex := 0; argIndex < c.Args().Len()-1; argIndex++ {
			match, _ := filepath.Glob(c.Args().Get(argIndex))
			matches = append(matches, match...)
		}
		matches = uniqueNonEmptyElementsOf(matches)
	}

	return nil
}

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
	// fmt.Println("This is string buffer: " + stringBuffer.String())
	return stringBuffer.String()
}

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

func isPipe() bool {
	info, _ := os.Stdin.Stat()
	flag := false
	if info.Mode()&os.ModeNamedPipe != 0 {
		// we have a pipe input
		flag = true
	}
	return flag
}
