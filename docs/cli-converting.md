---
id: cli-converting
title: Converting Sequence Files
---

One of Poly's surprisingly super useful features is being able to convert gbk and gff files to JSON and back again.

I know this sounds trivial but parsing and writing these legacy files isn't that simple and not many languages provide libraries for it. Our genbank parser is 500+ lines long and we're still adding to it!

## Examples
### Converting files to JSON
If you just want to turn the ugly gff and genbank files in a directory into JSON files so that you can push them to a database or use them with any language of your choice without having to write an annoying parser just run this inside any given directory.

```bash
poly convert -o json *.gbk *.gb *.gff
```
This will find every .gbk .gb and .gff file, parse them, and then write out a JSON file with a corresponding name.

### Converting pipes to JSON

Here's how you can pipe to standard input and then redirect to a new json file.

```bash
cat bsub.gbk | poly c -i gbk -o json > bsub.json
```

Notice the `-i` flag that specifies `gbk` input and the `-o` flag that specifies `json` output.

Convert output always defauts to `json` and can infer input from file paths but in the case of **piping input must always be specified using the `-i` flag**.

## More Info
For more info about poly convert and its usage try running `poly help convert` or `poly help c` from your command line.