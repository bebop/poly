package main

import (
	"bytes"
	"io/ioutil"
	"os"
	"path/filepath"
	"testing"

	"github.com/TimothyStiles/poly"
	"github.com/TimothyStiles/poly/io/genbank"
	"github.com/TimothyStiles/poly/io/polyjson"
	"github.com/google/go-cmp/cmp"
	"github.com/google/go-cmp/cmp/cmpopts"
)

/******************************************************************************
Oct, 15, 2020

Testing command line utilities can be annoying.

The way poly does it is by using the cli.app  object to spoof input and output
via cli.app.reader and cli.app.writer. This is the only way to get true stack
traceable coverage.
******************************************************************************/

var testFilePaths = []string{"../../data/puc19.gbk", "../../data/ecoli-mg1655-short.gff", "../../data/sample.json"}

func TestMain(t *testing.T) {
	rescueStdout := os.Stdout
	_, w, _ := os.Pipe()
	os.Stdout = w

	arg := os.Args[0:1]
	os.Args = append(arg, "-h")
	main()
	os.Args = os.Args[0:1]
	w.Close()
	os.Stdout = rescueStdout
}
func TestConvertPipe(t *testing.T) {

	for _, match := range testFilePaths {
		extension := filepath.Ext(match)[1:]
		var writeBuffer bytes.Buffer
		app := application()
		app.Writer = &writeBuffer

		args := os.Args[0:1]                                    // Name of the program.
		args = append(args, "c", "-i", extension, "-o", "json") // Append a flag

		file, _ := ioutil.ReadFile(match)
		app.Reader = bytes.NewReader(file)

		err := app.Run(args)

		if err != nil {
			t.Fatalf("Run error: %s", err)
		}

		// getting test sequence from non-pipe io to compare against io to stdout
		baseTestSequence := parseExt(match)
		pipeOutputTestSequence := polyjson.Parse(writeBuffer.Bytes())

		if diff := cmp.Diff(baseTestSequence, pipeOutputTestSequence, cmpopts.IgnoreFields(poly.Feature{}, "ParentSequence")); diff != "" {
			t.Errorf(" mismatch converting from %q to json (-want +got):\n%s", extension, diff)
		}
	}

}

func TestConvertIO(t *testing.T) {
	for _, match := range testFilePaths {
		extension := filepath.Ext(match)[1:]
		var writeBuffer bytes.Buffer
		app := application()
		app.Writer = &writeBuffer

		args := os.Args[0:1] // Name of the program.
		args = append(args, "c", "-i", extension, "-o", extension)
		file, _ := ioutil.ReadFile(match)
		app.Reader = bytes.NewReader(file)

		err := app.Run(args)

		if err != nil {
			t.Fatalf("Run error: %s", err)
		}

		// getting test sequence from non-pipe io to compare against io to stdout
		baseTestSequence := parseExt(match)

		pipeOutputTestSequence := parseFlag(writeBuffer.Bytes(), extension)

		if diff := cmp.Diff(baseTestSequence, pipeOutputTestSequence, cmpopts.IgnoreFields(poly.Feature{}, "ParentSequence")); diff != "" {
			t.Errorf(" mismatch reading and writing %q (-want +got):\n%s", extension, diff)
		}
	}
}

func TestConvertWriteFile(t *testing.T) {

	for _, match := range testFilePaths {
		extension := filepath.Ext(match)
		testOutputPath := "../../data/test" + extension

		loopApp := application()

		args := os.Args[0:1] // Name of the program.
		args = append(args, "c", "-o", testOutputPath, match)

		err := loopApp.Run(args)

		if err != nil {
			t.Fatalf("Run error: %s", err)
		}

		// getting test sequence from non-pipe io to compare against io to stdout
		baseTestSequence := parseExt(match)

		outputSequence := parseExt(testOutputPath)
		os.Remove(testOutputPath)

		if diff := cmp.Diff(baseTestSequence, outputSequence, cmpopts.IgnoreFields(poly.Feature{}, "ParentSequence")); diff != "" {
			t.Errorf(" mismatch reading and writing %q (-want +got):\n%s", extension, diff)
		}
	}
}

// The test below passes on macosx but not ubuntu. Could someone try debugging it on a device running ubuntu?
func TestConvertWarning(t *testing.T) {

	testFilePath := "../../data/puc19.gbk"
	var writeBuffer bytes.Buffer

	app := application()
	app.Writer = &writeBuffer

	// testing file matching hash
	args := os.Args[0:1]                                       // Name of the program.
	args = append(args, "c", "-o", testFilePath, testFilePath) // Append a flag

	err := app.Run(args)

	if err != nil {
		t.Fatalf("Run error: %s", err)
	}

	if writeBuffer.Len() == 0 {
		t.Error("TestConvertWarning did not write output to desired writer.")
	}

	warningOutputString := writeBuffer.String()

	warningTestString := "WARNING: " + "input path and output path match. File: " + testFilePath + ". Skipping to prevent possible data loss. Try providing a full path with a different name using the -o flag.\n"
	if warningOutputString != warningTestString {
		t.Errorf("TestConvertWarning has failed. Returned %q, want %q", warningOutputString, warningTestString)
	}
}

func TestConvertFile(t *testing.T) {

	app := application()

	args := os.Args[0:1]                                                                      // Name of the program.
	args = append(args, "c", "-o", "json", "../../data/puc19.gbk", "../../data/t4_intron.gb") // Append a flag

	err := app.Run(args)

	if err != nil {
		t.Fatalf("Run error: %s", err)
	}

	puc19InputTestSequence := genbank.Read("../../data/puc19.gbk")
	puc19OutputTestSequence := polyjson.Read("../../data/puc19.json")

	//clearing test data.
	os.Remove("../../data/puc19.json")

	// compared input gff from resulting output json. Fail test and print diff if error.
	if diff := cmp.Diff(puc19InputTestSequence, puc19OutputTestSequence, cmpopts.IgnoreFields(poly.Feature{}, "ParentSequence")); diff != "" {
		t.Errorf(" mismatch from concurrent gbk input test (-want +got):\n%s", diff)
	}

	t4InputTestSequence := genbank.Read("../../data/t4_intron.gb")
	t4OutputTestSequence := polyjson.Read("../../data/t4_intron.json")

	// clearing test data.
	os.Remove("../../data/t4_intron.json")

	// compared input gbk from resulting output json. Fail test and print diff if error.
	if diff := cmp.Diff(t4InputTestSequence, t4OutputTestSequence, cmpopts.IgnoreFields(poly.Feature{}, "ParentSequence")); diff != "" {
		t.Errorf(" mismatch from concurrent gbk input test (-want +got):\n%s", diff)
	}
}

func TestHashFile(t *testing.T) {

	testFilePath := "../../data/puc19.gbk"
	puc19GbkBlake3Hash := "v1_DCD_4b0616d1b3fc632e42d78521deb38b44fba95cca9fde159e01cd567fa996ceb9"
	var writeBuffer bytes.Buffer

	app := application()
	app.Writer = &writeBuffer

	// testing file matching hash
	args := os.Args[0:1]                      // Name of the program.
	args = append(args, "hash", testFilePath) // Append a flag

	err := app.Run(args)

	if err != nil {
		t.Fatalf("Run error: %s", err)
	}

	if writeBuffer.Len() == 0 {
		t.Error("TestHash did not write output to desired writer.")
	}

	hashOutputString := writeBuffer.String()

	// use same function that outputs hash results to format test case.
	testHashString := formatHashOutput(puc19GbkBlake3Hash, testFilePath)
	if hashOutputString != testHashString {
		t.Errorf("TestHashFile has failed. Returned %q, want %q", hashOutputString, puc19GbkBlake3Hash)
	}

}

func TestHashPipe(t *testing.T) {

	puc19GbkBlake3Hash := "v1_DCD_4b0616d1b3fc632e42d78521deb38b44fba95cca9fde159e01cd567fa996ceb9"
	testFilePath := "../../data/puc19.gbk"

	var writeBuffer bytes.Buffer

	// create a mock application
	app := application()
	app.Writer = &writeBuffer
	file, _ := ioutil.ReadFile(testFilePath)
	app.Reader = bytes.NewReader(file)

	args := os.Args[0:1]                     // Name of the program.
	args = append(args, "hash", "-i", "gbk") // Append a flag

	err := app.Run(args)

	if err != nil {
		t.Fatalf("Run error: %s", err)
	}

	hashOutputString := writeBuffer.String()

	// use same function that outputs hash results to format test case.
	testHashString := formatHashOutput(puc19GbkBlake3Hash, "-")

	if hashOutputString != testHashString {
		t.Errorf("TestHashPipe has failed. Returned %q, want %q", hashOutputString, puc19GbkBlake3Hash)
	}
}
func TestInputParameterForConvertErrorsIfNotPipe(t *testing.T) {
	app := application()
	args := append(os.Args[0:1], "convert", "-i", "json", "../../data/puc19.gbk")
	err := app.Run(args)

	if err != errIllegalInputFlag {
		t.Fatal("passing '-i' should throw an error if we are not reading from a pipe")
	}
}

func TestInputParameterForHashErrorsIfNotPipe(t *testing.T) {
	app := application()
	args := append(os.Args[0:1], "hash", "-i", "json", "../../data/puc19.gbk")
	err := app.Run(args)

	if err != errIllegalInputFlag {
		t.Fatal("passing '-i' should throw an error if we are not reading from a pipe")
	}
}
