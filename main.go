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
						Usage: "Specify file output type. Options are Gff, gbk/gb, and json. Defaults to json.",
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
			{
				Name:    "hash",
				Aliases: []string{"ha"},
				Usage:   "Hash a sequence while accounting for circularity.",
				Flags: []cli.Flag{
					&cli.StringFlag{
						Name:  "f",
						Value: "blake3",
						Usage: "Specify hash function type. Has many options. Blake3 is probably fastest.",
					},
					&cli.StringFlag{
						Name:  "i",
						Value: "json",
						Usage: "Specify file input type. Options are Gff, gbk/gb, and json.",
					},
					&cli.StringFlag{
						Name:  "o",
						Value: "string",
						Usage: "Specify output type. Options are string and json. Defaults to string.",
					},
					&cli.BoolFlag{
						Name:  "stdout",
						Value: false,
						Usage: "Will write to standard out whenever applicable. Defaults to false.",
					},
				},
				Action: func(c *cli.Context) error {
					hash(c)
					return nil
				},
			},
		}}

	err := app.Run(os.Args)
	if err != nil {
		log.Fatal(err)
	}
}
