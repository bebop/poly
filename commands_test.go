package main

import (
	"os"
	"os/exec"
	"testing"

	"github.com/google/go-cmp/cmp"
)

func TestConvert(t *testing.T) {

	// testing redirected pipe output
	command := "cat data/bsub.gbk | poly c -i gbk -o json > data/converttest.json"
	exec.Command("bash", "-c", command).Output()

	// getting test sequence from non-pipe io to compare against redirected io
	baseTestSequence := ReadGbk("data/bsub.gbk")
	outPutTestSequence := ReadJSON("data/converttest.json")

	// cleaning up test data
	os.Remove("data/converttest.json")

	// diff original sequence and converted sequence reread back into
	// AnnotatedSequence format. If there's an error fail test and print diff.
	if diff := cmp.Diff(baseTestSequence, outPutTestSequence); diff != "" {
		t.Errorf(" mismatch from convert pipe input test (-want +got):\n%s", diff)
	}

	//clearing possibly still existent prior test data.
	os.Remove("data/ecoli-mg1655.json")
	os.Remove("data/bsub.json")

	// testing multithreaded non piped output
	command = "poly c -o json data/bsub.gbk data/ecoli-mg1655.gff"
	exec.Command("bash", "-c", command).Output()

	ecoliInputTestSequence := ReadGff("data/ecoli-mg1655.gff")
	ecoliOutPutTestSequence := ReadJSON("data/ecoli-mg1655.json")

	//clearing test data.
	os.Remove("data/ecoli-mg1655.json")

	// compared input gff from resulting output json. Fail test and print diff if error.
	if diff := cmp.Diff(ecoliInputTestSequence, ecoliOutPutTestSequence); diff != "" {
		t.Errorf(" mismatch from concurrent gff input test (-want +got):\n%s", diff)
	}

	bsubInputTestSequence := ReadGbk("data/bsub.gbk")
	bsubOutPutTestSequence := ReadJSON("data/bsub.json")

	// clearing test data.
	os.Remove("data/bsub.json")

	// compared input gbk from resulting output json. Fail test and print diff if error.
	if diff := cmp.Diff(bsubInputTestSequence, bsubOutPutTestSequence); diff != "" {
		t.Errorf(" mismatch from concurrent gbk input test (-want +got):\n%s", diff)
	}

}
