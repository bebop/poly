# (Poly)merase &middot; [![Contributor Covenant](https://img.shields.io/badge/Contributor%20Covenant-v2.0%20adopted-ff69b4.svg)](CODE_OF_CONDUCT.md) [![GitHub license](https://img.shields.io/badge/license-MIT-blue.svg)](https://github.com/facebook/react/blob/master/LICENSE) ![Test](https://github.com/TimothyStiles/poly/workflows/Test/badge.svg) ![PRs Welcome](https://img.shields.io/badge/PRs-welcome-brightgreen.svg)

Poly is a Go library and command line utility for engineering organisms.

* **Fast:** Poly is fast and scalable.

* **Modern:** Poly tackles issues that other libraries and utilities just don't. From general codon optimization and primer design to circular sequence hashing. All written in a language that was designed to be fast, scalable, and easy to develop in and maintain. Did we say it was fast?

* **Reproducible:** Poly is well tested and designed to be used in industrial, academic, and hobbyist settings. No more copy and pasting strings into random websites to process the data you need.

* **Ambitious:** Poly's goal is to be the most complete, open, and well used collection of computational synthetic biology tools ever assembled. If you like our dream and want to support us please star this repo, request a feature, or open a pull request.

## Installation

### Installing Poly as a Go library and Command Line Utility
This assumes you already have a working Go environment, if not please see
[this page](https://golang.org/doc/install) first.

`go get` *will always pull the latest released version from the prime branch.*

```bash
go get github.com/TimothyStiles/poly
```

### Installing Poly as a Command Line Utility

Poly ships many binaries for many different operating systems and package managers thanks to the wonderful work of the [go releaser](https://goreleaser.com/) team. You can check out our [releases page](https://github.com/TimothyStiles/poly/releases) on github or install via package manager for your OS with the instructions below.

#### Mac OS

```bash
brew install timothystiles/poly/poly
```

#### Linux

```bash
brew install timothystiles/poly/poly
```

#### Windows

[Coming soon...](https://github.com/TimothyStiles/poly/issues/16)

## Building Poly from Scratch

This assumes you already have a working Go environment, if not please see
[this page](https://golang.org/doc/install) first.

```bash
git clone https://github.com/TimothyStiles/poly.git && cd poly && go build && go install
```

## Examples


### Command Line Interface

Converting a .gbk to .json using pipes. 

```bash
cat bsub.gbk | poly c -i gbk -o json > bsub.json
```

Here's how you can non-destructively copy and convert every genbank and gff file into JSON files. The -o flag defaults to json and can also be used to specify gff as output.

```bash
poly c -o json *.gbk *.gb *.gff
```

### Go Library

Here's how you'd read in a file from it's path:

```Go
bsubAnnotatedSequence := ReadGbk("data/bsub.gbk")
```



## Documentation

If you want to see a ton of great examples of what poly can do you should check out our [docs site](https://timothystiles.github.io/poly/docs/).

You can also learn more about poly or a sub command like convert using the -h flag which will provide more documentation.

## Contributing

### Code of Conduct

Poly has adopted a [Code of Conduct](CODE_OF_CONDUCT.md). Please read the full text so you can understand what we're all about and please remember to be excellent to each other!

### Contributing Guide

Poly also has a [contributor's guide](CONTRIBUTING.md). Please read through it before you start hacking away and pushing contributions to this fine codebase.

### License
* [MIT](LICENSE)

* Copyright (c) 2020 Timothy Stiles
