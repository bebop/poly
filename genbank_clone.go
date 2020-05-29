package main

import (
	"fmt"
	"io"
	"log"
	"net/http"
	"os"
	"path/filepath"
	"regexp"

	"github.com/PuerkitoBio/goquery"
)

const genbankURL = "https://ftp.ncbi.nlm.nih.gov/genbank/"

func genbankClone() {
	// Request the HTML page.
	genbankDataPath := "data/genbank/"
	if _, err := os.Stat(genbankDataPath); os.IsNotExist(err) {
		os.MkdirAll("data/genbank", 0777)
	}

	res, err := http.Get(genbankURL)
	if err != nil {
		log.Fatal(err)
	}
	defer res.Body.Close()
	if res.StatusCode != 200 {
		log.Fatalf("status code error: %d %s", res.StatusCode, res.Status)
	}

	// Load the HTML document
	doc, err := goquery.NewDocumentFromReader(res.Body)
	if err != nil {
		log.Fatal(err)
	}

	reg, err := regexp.Compile("seq.gz$")
	// Find the review items

	doc.Find("a").Each(func(i int, s *goquery.Selection) {
		// For each item found, get the band and title

		genbankFileName, ok := s.Attr("href")
		if ok && reg.MatchString(genbankFileName) {
			genbankAddress := genbankURL + genbankFileName
			localDownloadPath := filepath.Join(genbankDataPath, genbankFileName)
			err := downloadFile(localDownloadPath, genbankAddress)
			if err != nil {
				panic(err)
			}
			fmt.Println("Downloaded: " + genbankAddress)
		}
	})
}

// taken from https://golangcode.com/download-a-file-from-a-url/
func downloadFile(filepath string, url string) error {

	// Get the data
	resp, err := http.Get(url)
	if err != nil {
		return err
	}
	defer resp.Body.Close()

	// Create the file
	out, err := os.Create(filepath)
	if err != nil {
		return err
	}
	defer out.Close()

	// Write the body to file
	_, err = io.Copy(out, resp.Body)
	return err
}
