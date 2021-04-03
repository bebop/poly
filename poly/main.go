package main

import (
	"log"
	"os"

	"github.com/urfave/cli/v2"
)

/******************************************************************************
Oct, 15, 2020

This file is special because it is the entry point for our command line utility.
It also acts as a general template that outlines everything available to the user.

Initial argparsing and app definition is done entirely through
"github.com/urfave/cli/v2" for which you can find the docs here:

https://github.com/urfave/cli/blob/master/docs/v2/manual.md

Essentially poly's app is defined via the &cli.App{} struct which you initialize
with data needed to run your app. In our case we're providing it Name, Usage, Flags,
and Commands at the top level. Commands can also be nested to provide n-level sub commands.

When naming new flags please make sure they don't collide with already existent
flags and try to follow these naming conventions:

http://www.catb.org/~esr/writings/taoup/html/ch10s05.html

Happy hacking,
Tim

******************************************************************************/

// main is well... the main entry point for our command line app. We seperate it from the actual &cli.App to help with testing.
func main() {
	run(os.Args)
}

// run is seperated from main and application for debugging's sake.
func run(args []string) {

	app := application()
	err := app.Run(os.Args) // run app and log errors
	if err != nil {
		log.Fatal(err)
	}

}

// application is a function that defines instances of our app. it's where we template commands and where initial arg parsing occurs.
func application() *cli.App {

	app := &cli.App{
		Name:  "poly",
		Usage: "A command line utility for engineering organisms.",

		// This is where you define global flags. Each sub command can also have its own flags that overide globals
		Flags: []cli.Flag{

			&cli.BoolFlag{
				Name:  "y",
				Usage: "Answers yes for all confirmations before doing something possibly destructive.",
			},

			&cli.StringFlag{
				Name:  "i",
				Usage: "Specify file input type or input path.",
			},

			&cli.StringFlag{
				Name:  "o",
				Usage: "Specify file output type or output path.",
			},
		},

		// This is where you start defining subcommands there's a lot of spacing to enhance readability since these nested brackets can be a little much.
		Commands: []*cli.Command{

			// defining the kind of *cli.Context our subcommand will use.
			{
				Name:    "convert",
				Aliases: []string{"c"},
				Usage:   "Convert a single file or set of files from one type to another. Genbank to Json, Json to Gff, etc.",

				// defining flags for this specific sub command
				Flags: []cli.Flag{

					&cli.StringFlag{
						Name:  "o",
						Value: "json",
						Usage: "Specify file output type or path. Defaults to json.",
					},

					&cli.StringFlag{
						Name:  "i",
						Value: "",
						Usage: "Specify file input type. Options are Gff, gbk/gb, and json. Defaults to none.",
					},
				},
				// where we provide the actual function that is called by the subcommand.
				Action: func(c *cli.Context) error {
					convertCommand(c)
					return nil
				},
			},

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
						Usage: "Specify file input type.",
					},

					&cli.StringFlag{
						Name:  "o",
						Value: "string",
						Usage: "Specify output type. Defaults to string.",
					},
				},
				Action: func(c *cli.Context) error {
					hashCommand(c)
					return nil
				},
			},
			{
				Name:    "translate",
				Aliases: []string{"tr"},
				Usage:   "Translate a coding sequence into an amino acid string",

				Flags: []cli.Flag{
					&cli.StringFlag{
						Name:    "ct",
						Aliases: []string{"--codon-table"},
						Value:   "Standard",
						Usage:   "Specify codon table used for organism. Defaults to Standard NCBI table.",
					},
					&cli.StringFlag{
						Name:  "i",
						Value: "string",
						Usage: "Specify file input type.",
					},
				},
				Action: func(c *cli.Context) error {
					translateCommand(c)
					return nil
				},
			},
			{
				Name:    "optimize",
				Aliases: []string{"op"},
				Usage:   "Optimize a sequence for a specific organism.",

				Flags: []cli.Flag{
					&cli.StringFlag{
						Name:  "i",
						Value: "string",
						Usage: "Specify file input type.",
					},
					&cli.StringFlag{
						Name:    "ct",
						Aliases: []string{"--codon-table"},
						Value:   "1",
						Usage:   "Specify codon table used for organism. Defaults to Standard NCBI table.",
					},
					&cli.BoolFlag{
						Name:    "aa",
						Aliases: []string{"--amino-acid"},
						Value:   false,
						Usage:   "Specify that the input sequence is an amino acid sequence.",
					},
					&cli.StringFlag{
						Name:    "wt",
						Aliases: []string{"--weight-table"},
						Usage:   "Specify file to weigh table with",
					},
				},
				Action: func(c *cli.Context) error {
					optimizeCommand(c)
					return nil
				},
			},
      {
				Name:    "amplify",
				Aliases: []string{"amp"},
				Usage:   "Design primers to amplify a particular sequence of DNA.",

				Flags: []cli.Flag{
					&cli.StringFlag{
						Name:  "primer_for",
						Value: "string",
						Usage: "Name the output forward primer in primers.txt",
					},
          &cli.StringFlag{
						Name:  "primer_rev",
						Value: "string",
						Usage: "Name the output reverse primer in primers.txt",
					},
          &cli.StringFlag{
            Name:  "amplicon",
						Value: "string",
						Usage: "Amplify a particular string",
          },
          &cli.StringFlag{
            Name:  "range",
						Value: "string",
						Usage: "Amplify a particular range of sequence",
          },
          &cli.StringFlag{
            Name:  "no_amplify",
						Value: "string",
						Usage: "Prevent primers from amplifying a sequence from a different file",
          },
          &cli.StringFlag{
            Name:  "validate",
						Value: "string",
						Usage: "Specify file input type.",
          },
          &cli.StringFlag{
            Name:  "size",
						Value: "string",
						Usage: "Specify file input type.",
          },
          &cli.StringFlag{
            Name:  "coverage",
						Value: "string",
						Usage: "Specify file input type.",
          },
          &cli.StringFlag{
            Name:  "overlap",
						Value: "string",
						Usage: "Specify file input type.",
          },
				},
				Action: func(c *cli.Context) error {
					// amplifyCommand()
					return nil
				},
			},
		}, // subcommands list ends here
	} // app definition ends here

	return app
}
