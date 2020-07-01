---
id: installation
title: Installing Poly
---

Poly can be used in two ways.

1. As a Go library where you have finer control and can make magical things happen.
2. As a command line utility where you can bash script your way to greatness and make DNA go vrroooom.

## Installing Poly as a Go library and Command Line Utility
This assumes you already have a working Go environment, if not please see
[this page](https://golang.org/doc/install) first.

`go get` *will always pull the latest released version from the prime branch.*

```bash
go get github.com/TimothyStiles/poly
```

## Installing Poly as a Command Line Utility

Poly ships many binaries for many different operating systems and package managers thanks to the wonderful work of the [go releaser](https://goreleaser.com/) team. You can check out our [releases page](https://github.com/TimothyStiles/poly/releases) on github or install via package manager for your OS with the instructions below.

### Mac OS
```bash
brew install timothystiles/poly/poly
```

### Linux
```bash
sudo snap install --classic gopoly
```

### Windows
[Coming soon...](https://github.com/TimothyStiles/poly/issues/16)

## Building Poly from Scratch
This assumes you already have a working Go environment, if not please see
[this page](https://golang.org/doc/install) first.

```bash
git clone https://github.com/TimothyStiles/poly.git && cd poly && go build && go install
```