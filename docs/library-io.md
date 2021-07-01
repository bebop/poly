---
id: library-io
title: Sequence Input Output
---

At the center of `poly`'s annotated sequence support is the `Sequence` struct. Structs are kind of Go's answer to objects in other languages. They provide a way of making custom datatypes and methods for developers to use. More on that [here](https://tour.golang.org/moretypes/2), [here](https://gobyexample.com/methods), and [here](https://www.golang-book.com/books/intro/9).

Anywho. `poly` centers around reading in various annotated sequence formats like genbank, or gff and parsing them into an `Sequence` to do stuff with them. Whether that's being written out to JSON or being used by `poly` itself. Here are some examples.

## Readers

For all supported file formats `poly` supports a reader. A reader is a function literally named `ReadJSON(path)`, `ReadGbk(path)`, or `ReadGff(path)` that takes one argument - a filepath where your file is located, and returns an `Sequence` struct.

```go
  bsubAnnotatedSequence := ReadGbk("data/bsub.gbk")
  ecoliAnnotatedSequence := ReadGff("data/ecoli-mg1655.gff")
  puc19AnnotatedSequence := ReadJSON("data/puc19static.json")
```

These Sequence structs contain all sorts of goodies but can be broken down into three sub main structs. `Sequence.Meta`, `Sequence.Features`, and `Sequence.Sequence`.

> Before we move on with the rest of IO I think it'd be good to go over these sub structs in the next section but of course you can skip to [writers](#writers) if you'd like.

## Sequence structs

Like I just said these Sequence structs contain all sorts of goodies but can be broken down into three main sub structs:

  * [Sequence.Meta](#annotatedsequencemeta)
  * [Sequence.Features](#annotatedsequencefeatures)
  * [Sequence.Sequence](#annotatedsequencesequence)

Here's how the Sequence struct is actually implemented as of [commit c4fc7e](https://github.com/TimothyStiles/poly/blob/c4fc7e6f6cdbd9e5ed2d8ffdbeb206d1d5a8d720/io.go#L108).

```go
  // Sequence holds all sequence information in a single struct.
  type Sequence struct {
    Meta     Meta
    Features []Feature
    Sequence Sequence
  }
```

> You can check out the original implementation [here](https://github.com/TimothyStiles/poly/blob/c4fc7e6f6cdbd9e5ed2d8ffdbeb206d1d5a8d720/io.go#L108) but I warn you that this is a snapshot and likely has been updated since last writing.

### Sequence.Meta

The Meta substruct contains various meta information about whatever record was parsed. Things like name, version, genbank references, etc.

So if I wanted to get something like the Genbank Accession number for a Sequence I'd get it like this:

```go
  bsubAnnotatedSequence := ReadGbk("data/bsub.gbk")
  bsubAccessionNumber := bsubAnnotatedSequence.Meta.Accession
```

Same goes for a lot of other stuff:

```go
  bsubAnnotatedSequence := ReadGbk("data/bsub.gbk")

  // Getting genbank division
  bSubGenbankDivision := bsubAnnotatedSequence.Meta.GenbankDivision

  // Getting genbank organism
  bsubGenbankOrganism := bsubAnnotatedSequence.Meta.Organism
```

Here's how the Meta struct is actually implemented in [commit c4fc7e](https://github.com/TimothyStiles/poly/blob/c4fc7e6f6cdbd9e5ed2d8ffdbeb206d1d5a8d720/io.go#L34) which is the latest as of writing.

```go
  // Meta Holds all the meta information of an Sequence struct.
  type Meta struct {
    Name            string
    GffVersion      string
    RegionStart     int
    RegionEnd       int
    Size            int
    Type            string
    GenbankDivision string
    Date            string
    Definition      string
    Accession       string
    Version         string
    Keywords        string
    Organism        string
    Source          string
    Origin          string
    Locus           Locus
    References      []Reference
    Primaries       []Primary
}
```

You'll notice that there are actually three more substructs towards the bottom. They hold extra genbank specific information that's handy to have grouped together. More about how genbank files are structered can be found [here](https://www.ncbi.nlm.nih.gov/Sitemap/samplerecord.html).

### Sequence.Features

The `Features` substruct is actually a slice (golang term for what is essentially a dynamic length list) of `Feature` structs that can be iterated through. For example if you wanted to iterate through an `Sequence`'s features and get their name (i.e GFP) and type (i.e CDS) you'd do it like this.

```go
  bsubAnnotatedSequence := ReadGbk("data/bsub.gbk")
  for feature in := range bsubAnnotatedSequence.Features {
    fmt.Println(feature.Name + " " + feature.Type)
  }
```

The `Feature` struct has about 10 or so fields which you can learn more about from this section in [commit c4fc7e](https://github.com/TimothyStiles/poly/blob/c4fc7e6f6cdbd9e5ed2d8ffdbeb206d1d5a8d720/io.go#L80).

### Sequence.Sequence

The Sequence Sequence substruct is by far the most basic and critical. Without it well, you ain't go no DNA. The substruct itself has 4 simple fields.

```go
  // Sequence holds raw sequence information in an Sequence struct.
  type Sequence struct {
    Description  string
    Hash         string
    HashFunction string
    Sequence     string
  }
```

The `Description`, `Hash`, and `HashFunction` are at all identifying fields of the Sequence string. The `Description` is the same kind of short description you'd find in a `fasta` or `fastq` file. The `Hash` and `HashFunction` are used to create a unique identifier specify to the sequence string which you'll learn more about in the next chapter on sequence hashing.

To get an Sequence sequence you can address it like so:

```go
  bsubAnnotatedSequence := ReadGbk("data/bsub.gbk")
  bsubSequenceString := bsubAnnotatedSequence.Sequence.Sequence
```

## Writers

 `poly` tries to supply a writer for all supported file formats that have a reader.

 Writers take two arguments. The first is an Sequence struct, the second is a path to write out to.

 ```go
  // getting Sequence(s) to write out again.
  bsubAnnotatedSequence := ReadGbk("data/bsub.gbk")

  // writing out gbk file input as json.
  WriteJSON(bsubAnnotatedSequence, "data/bsub.json")

  // writing out gbk file input as gff.
  WriteGff(bsubAnnotatedSequence, "data/bsub.gff")

  // writing out gbk file input as gff.
  WriteGff(bsubAnnotatedSequence, "data/bsub.gbk")

 ```

 Note that interop between file formats that aren't JSON<->* may be a little buggy. For example .gff doesn't contain a lot of sections that genbank has and genbank uses a system of sequence locations that isn't supported by the gff3 standard.

## Parsers

`poly` parsers are what actually parse input files from a string without any of the system IO. This is particularly useful if you're like me and have an old database holding genbank files as strings. You can take those strings from a database or whatever and just pass them to `ParseGbk()`, or `ParseGff()` and they'll convert them into Sequence structs.

```go
  puc19AnnotatedSequence := ParseGbk("imagine this is actually a gbk in string format.")
```

That's it. The reason we don't have a `ParseJSON()` is that golang, like almost all languages, already has a standard library to handle that.

## Builders

`poly` builders take Sequence structs and use them to build strings for different file formats. 

```go
  // generating an Sequence struct from a gff file.
  ecoliAnnotatedSequence := ReadGff("data/ecoli-mg1655.gff")

  // generating a gff string that then can be piped to stdout or written to a database.
  ecoliGffString := BuildGff(ecoliAnnotatedSequence)

  ecoliGbkString := Build(ecoliAnnotatedSequence)

```

## More Info

For more info about IO check out the most current [code](https://github.com/TimothyStiles/poly/blob/prime/io.go) and [tests](https://github.com/TimothyStiles/poly/blob/prime/io_test.go) and maybe open a [pull request](https://github.com/TimothyStiles/poly/issues/new/choose) for improvements.