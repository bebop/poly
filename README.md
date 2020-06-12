# (Poly)merase &middot; [![Contributor Covenant](https://img.shields.io/badge/Contributor%20Covenant-v2.0%20adopted-ff69b4.svg)](CODE_OF_CONDUCT.md)  [![GitHub license](https://img.shields.io/badge/license-MIT-blue.svg)](https://github.com/facebook/react/blob/master/LICENSE) ![Build Status](https://travis-ci.org/TimothyStiles/poly.svg?branch=master) ![PRs Welcome](https://img.shields.io/badge/PRs-welcome-brightgreen.svg) 
Poly is a Go library and command line utility for engineering organisms.

* **Fast:** Poly is fast and scalable.

* **Modern:** Poly tackles issues that other libraries and utilities just don't. From general codon optimization and primer design to circular sequence hashing. All written in a language that was designed to be fast, scalable, and easy to develop in and maintain. Did we say it was fast?

* **Reproducible:** Poly is well tested and designed to be used in industrial, academic, and hobbyist settings. No more copy and pasting strings into random websites to process the data you need.

* **Ambitious:** Poly's goal is to be the most complete, open, and well used collection of computational synthetic biology tools ever assembled. If you like our dream and want to support us please star this repo, request a feature, or open a pull request.

## Installation

```bash
go get github.com/TimothyStiles/poly
```

## Examples


### cli

At the moment Poly can be used to batch convert annotated sequence files like gff and genbank into JSON. You can do this by either providing wild card paths as arguments or piping to stdin.

Here's how you can pipe to standard input and then redirect to a new json file.

```bash
cat bsub.gbk | poly c -i gbk -o json > bsub.json
```

Here's how you can non-destructively copy and convert every genbank and gff file into JSON files. -o default to json and can also be used to specify gff as output.

```bash
poly c -o json *.gbk *.gb *gff
```

## Documentation

In progress, but most code is well commented.

You can also learn more about poly or a sub command like convert using the -h flag which will provide more documentation.

## Contributing

### Code of Conduct

Poly has adopted a [Code of Conduct](CODE_OF_CONDUCT.md). Please read the full text so you can understand what we're all about and please remember to be excellent to each other!

### Contributing Guide

Poly also has a [contributor's guide](CONTRIBUTING.md). Please read through it before you start hacking away and pushing contributions to this fine codebase.

### License
* [MIT](LICENSE)

* Copyright (c) 2020 Timothy Stiles
