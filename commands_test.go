package main

import (
	"bytes"
	"fmt"
	"io"
	"log"
	"os"
	"sync"
	"testing"
)

/******************************************************************************
Oct, 15, 2020

Testing command line utilities via subroutines can be annoying so
if you're doing it from the commandline be sure to compile first.
From the project's root directory often use:

go build && go install && go test -v ./...

To accurately test your commands you MUST make sure to rebuild and reinstall
before you run your tests. Otherwise your system version will be out of
date and will give you results using an older build.

Happy hacking,
Tim


TODO:

write subtest to check for empty output before merge
******************************************************************************/

func TestPipe(t *testing.T) {

	r, stdinWriter, err := os.Pipe()
	if err != nil {
	}
	fmt.Println(isPipe())

	origStdin := os.Stdin
	_, err = stdinWriter.Write([]byte(BuildGbk(ReadGbk("data/puc19.gbk"))))

	if err != nil {
		stdinWriter.Close()
		os.Stdin = origStdin
	}

	os.Stdin = r
	stdinWriter.Close()

	args := os.Args[0:1]                   // Name of the program.
	args = append(args, "ha", "-i", "gbk") // Append a flag
	results := captureOutput(func() { run(args) })
	fmt.Print(results)
	r.Close()
	os.Stdin = origStdin
}

// func TestConvert(t *testing.T) {

// 	if runtime.GOOS == "windows" {
// 		fmt.Println("TestConvert was not run and autopassed. Currently Poly command line support is not available for windows. See https://github.com/TimothyStiles/poly/issues/16.")
// 	} else {

// 		// testing redirected pipe output
// 		command := "cat data/puc19.gbk | poly c -i gbk -o json > data/converttest.json"
// 		exec.Command("bash", "-c", command).Output()

// 		// getting test sequence from non-pipe io to compare against redirected io
// 		baseTestSequence := ReadGbk("data/puc19.gbk")
// 		outputTestSequence := ReadJSON("data/converttest.json")

// 		// cleaning up test data
// 		os.Remove("data/converttest.json")

// 		// diff original sequence and converted sequence reread back into
// 		// AnnotatedSequence format. If there's an error fail test and print diff.
// 		if diff := cmp.Diff(baseTestSequence, outputTestSequence, cmpopts.IgnoreFields(Feature{}, "ParentAnnotatedSequence")); diff != "" {
// 			t.Errorf(" mismatch from convert pipe input test (-want +got):\n%s", diff)
// 		}

// 		//clearing possibly still existent prior test data.
// 		os.Remove("data/ecoli-mg1655.json")
// 		os.Remove("data/puc19.json")

// 		// testing multithreaded non piped output
// 		command = "poly c -o json data/puc19.gbk data/ecoli-mg1655.gff"
// 		exec.Command("bash", "-c", command).Output()

// 		ecoliInputTestSequence := ReadGff("data/ecoli-mg1655.gff")
// 		ecoliOutputTestSequence := ReadJSON("data/ecoli-mg1655.json")

// 		//clearing test data.
// 		os.Remove("data/ecoli-mg1655.json")

// 		// compared input gff from resulting output json. Fail test and print diff if error.
// 		if diff := cmp.Diff(ecoliInputTestSequence, ecoliOutputTestSequence, cmpopts.IgnoreFields(Feature{}, "ParentAnnotatedSequence")); diff != "" {
// 			t.Errorf(" mismatch from concurrent gff input test (-want +got):\n%s", diff)
// 		}

// 		bsubInputTestSequence := ReadGbk("data/puc19.gbk")
// 		bsubOutputTestSequence := ReadJSON("data/puc19.json")

// 		// clearing test data.
// 		os.Remove("data/puc19.json")

// 		// compared input gbk from resulting output json. Fail test and print diff if error.
// 		if diff := cmp.Diff(bsubInputTestSequence, bsubOutputTestSequence, cmpopts.IgnoreFields(Feature{}, "ParentAnnotatedSequence")); diff != "" {
// 			t.Errorf(" mismatch from concurrent gbk input test (-want +got):\n%s", diff)
// 		}

// 	}
// }

// func TestHash(t *testing.T) {
// 	if runtime.GOOS == "windows" {
// 		fmt.Println("TestHash was not run and autopassed. Currently Poly command line support is not available for windows. See https://github.com/TimothyStiles/poly/issues/16.")
// 	} else {

// 		puc19GbkBlake3Hash := "4b0616d1b3fc632e42d78521deb38b44fba95cca9fde159e01cd567fa996ceb9"

// 		// testing pipe input
// 		command := "cat data/puc19.gbk | poly hash -i gbk"
// 		hashOutput, _ := exec.Command("bash", "-c", command).Output()
// 		hashOutputString := strings.TrimSpace(string(hashOutput))

// 		if hashOutputString != puc19GbkBlake3Hash {
// 			t.Errorf("TestHash for piped input has failed. Returned %q, want %q", hashOutputString, puc19GbkBlake3Hash)
// 		}

// 		// testing regular input
// 		command = "poly hash data/puc19.gbk"
// 		hashOutput, _ = exec.Command("bash", "-c", command).Output()
// 		hashOutputString = strings.TrimSpace(string(hashOutput))

// 		if hashOutputString != puc19GbkBlake3Hash {
// 			t.Errorf("TestHash for regular input has failed. Returned %q, want %q", hashOutputString, puc19GbkBlake3Hash)
// 		}

// 		// testing json write output
// 		command = "poly hash -o json data/puc19.gbk"
// 		exec.Command("bash", "-c", command).Output()
// 		hashOutputString = ReadJSON("data/puc19.json").Sequence.Hash
// 		os.Remove("data/puc19.json")

// 		if hashOutputString != puc19GbkBlake3Hash {
// 			t.Errorf("TestHash for json write output has failed. Returned %q, want %q", hashOutputString, puc19GbkBlake3Hash)
// 		}

// 	}
// }

// func TestOptimizeCommand(t *testing.T) {
// 	if runtime.GOOS == "windows" {
// 		fmt.Println("TestOptimize was not run and autopassed. Currently Poly command line support is not available for windows. See https://github.com/TimothyStiles/poly/issues/16.")
// 	} else {

// 		gfpTranslation := "MASKGEELFTGVVPILVELDGDVNGHKFSVSGEGEGDATYGKLTLKFICTTGKLPVPWPTLVTTFSYGVQCFSRYPDHMKRHDFFKSAMPEGYVQERTISFKDDGNYKTRAEVKFEGDTLVNRIELKGIDFKEDGNILGHKLEYNYNSHNVYITADKQKNGIKANFKIRHNIEDGSVQLADHYQQNTPIGDGPVLLPDNHYLSTQSALSKDPNEKRDHMVLLEFVTAAGITHGMDELYK*"

// 		command := "echo " + gfpTranslation + " |" + "poly op -aa"
// 		optimizeOutput, _ := exec.Command("bash", "-c", command).Output()
// 		optimizeOutputString := strings.TrimSpace(string(optimizeOutput))
// 		translation := Translate(optimizeOutputString, DefaultCodonTablesByName["Standard"])

// 		if translation != gfpTranslation {
// 			t.Errorf("TestOptimizeCommand for json write output has failed. Returned %q, want %q", translation, gfpTranslation)
// 		}

// 	}

// }

// spoofPipe is a helper function to spoof stdin for testing pipe workflows without a bunch of cmd.exec calls
func spoofPipe(bytes []byte, function func()) string {
	r, stdinWriter, err := os.Pipe()
	if err != nil {
	}
	fmt.Println(isPipe())

	origStdin := os.Stdin
	_, err = stdinWriter.Write(bytes)

	if err != nil {
		stdinWriter.Close()
		os.Stdin = origStdin
	}

	os.Stdin = r
	stdinWriter.Close()

	args := os.Args[0:1]                   // Name of the program.
	args = append(args, "ha", "-i", "gbk") // Append a flag
	results := captureOutput(func() { function() })
	fmt.Print(results)
	r.Close()
	os.Stdin = origStdin

	return results
}
func captureOutput(f func()) string {
	reader, writer, err := os.Pipe()
	if err != nil {
		panic(err)
	}
	stdout := os.Stdout
	stderr := os.Stderr
	defer func() {
		os.Stdout = stdout
		os.Stderr = stderr
		log.SetOutput(os.Stderr)
	}()
	os.Stdout = writer
	os.Stderr = writer
	log.SetOutput(writer)
	out := make(chan string)
	wg := new(sync.WaitGroup)
	wg.Add(1)
	go func() {
		var buf bytes.Buffer
		wg.Done()
		io.Copy(&buf, reader)
		out <- buf.String()
	}()
	wg.Wait()
	f()
	writer.Close()
	return <-out
}
