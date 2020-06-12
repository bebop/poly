package main

import (
	"log"
	"os"

	"github.com/urfave/cli/v2"
)

func main() {
	app := &cli.App{
		Name:  "poly",
		Usage: "A command line utility for engineering organisms.",
		Flags: []cli.Flag{
			&cli.BoolFlag{
				Name:  "y",
				Usage: "Answers yes for all prompts.",
			},
		},
		Commands: []*cli.Command{
			{
				Name:    "convert",
				Aliases: []string{"c"},
				Usage:   "Convert a single file or set of files from one type to another. Genbank to Json, Json to Gff, etc.",
				Flags: []cli.Flag{
					&cli.StringFlag{
						Name:  "o",
						Value: "json",
						Usage: "Specify file output type. Options are Gff and json. Defaults to json.",
					},
					&cli.StringFlag{
						Name:  "i",
						Value: "",
						Usage: "Specify file input type. Options are Gff, gbk/gb, and json. Defaults to none.",
					},
				},
				Action: func(c *cli.Context) error {
					convert(c)
					return nil
				}},
		}}

	err := app.Run(os.Args)
	if err != nil {
		log.Fatal(err)
	}
}
