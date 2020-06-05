package main

import (
	"encoding/json"
	"io/ioutil"
)

// WriteJSON writes an AnnotatedSequence struct out to json.
func WriteJSON(annotatedSequence AnnotatedSequence, path string) {
	file, _ := json.MarshalIndent(annotatedSequence, "", " ")
	_ = ioutil.WriteFile(path, file, 0644)
}

// ReadJSON reads an AnnotatedSequence JSON file.
func ReadJSON(path string) AnnotatedSequence {
	file, err := ioutil.ReadFile(path)
	if err != nil {
		// return 0, fmt.Errorf("Failed to open file %s for unpack: %s", gzFilePath, err)
	}
	var annotatedSequence AnnotatedSequence
	json.Unmarshal([]byte(file), &annotatedSequence)
	return annotatedSequence
}
