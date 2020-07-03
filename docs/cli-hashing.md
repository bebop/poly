---
id: cli-hashing
title: Hashing Sequences
---

`poly` provides what is likely the only open source sequence hashing tool that handles circular sequences. By utilizing Booth's Least Rotation algorithm we're able to determistically rotate circular sequences to a fixed point which makes it possible to hash them. [More info on Booth's Least Rotation here](#more-info).

Hashes make incredibly powerful unique identifiers and with `poly hash` defaulting to the superfast blake3 hashing algorithm you can create them faster than ever before.

## Hashing a sequence from file

To hash a sequence from file all you have to do is this:

``` bash
poly hash bsub.gbk
```

`poly hash` will then parse the sequence string, rotate it to a deterministic point if it's circular, then hash it using the default blake3 algorithm returning something like this:

``` bash
949b2e18461fc354d989b14f8d4a58f710f3f46968b6bbffdbdc59a28ad77e83
```

## Hashing a sequence from a stream

`poly hash` can also accept input from streams to create a hash:

```bash
cat bsub.gbk | poly ha -i gbk
```

## Converting a sequence to JSON and hashing a sequence from a stream

```bash
cat bsub.gbk | poly ha -i gbk -o json > bsub.json
```

## Hashing multiple file inputs and writing strings to stdout

`poly hash` can also take all files in a directory and spit out their hashes along with their original file paths.

```bash
poly ha *.gbk *.gb *.gff
```

Result:

```bash
4031e1971acc8ff1bf0aa4ed623bc58beefc15e043075866a0854d592d80b28b   puc19.gbk
e35849d7d9d5476e84468f8527c1c8b8a0d4b6a2cf88d4329246b1cbba0920bc   sample.gbk
949b2e18461fc354d989b14f8d4a58f710f3f46968b6bbffdbdc59a28ad77e83   bsub.gbk
```

## Hashing multiple file inputs and writing out to JSON file

This is pretty much the same as `poly convert` but also hashes the sequence and stores the hash and meta info in the resulting jsons' Sequence struct.

```bash
poly ha -o json *.gbk *.gb *.gff
```

## Hashing multiple file inputs and streaming to stdout

I really woudn't recommend this but with the `--stdout` flag you can force all `json` output to be streamed to stdout. Useful if you want to start a bash process with hashes.

```bash
poly ha -o json --stdout *.gbk
```

## Hashing a sequence with a different hashing function

`poly hash` provides the option to use different hashing functions with the `-f` flag. Almost every hashing function available to Golang is available to `poly hash`. For example:

```bash
poly hash -f sha1 bsub.gbk
```

Will produce a sha1 hash. For a complete list of hashes and their flags check out the original source code [here](https://github.com/TimothyStiles/poly/blob/346e3eb58cdd74db14eba333ba428256f77c93b0/commands.go#L256). Hash flag values are case insensitive.

## Hashing with a system call

`poly hash` also provides the `no` argument to the function flag `-f`. This prints a rotated, unhashed sequence string to stdout to be piped as input into your system's hashing utility.

```bash
poly hash -f sha1 data/puc19.gbk
# returns: e5066a52a8b91eb8949b813347931f80e409b7c2

poly hash -f no data/puc19.gbk | shasum
# returns: e5066a52a8b91eb8949b813347931f80e409b7c2  -
```

## More Info

For more info about `poly hash` and its usage try running `poly help hash` or `poly help ha` from your command line.

For more info about circular sequence hashing and Booth's Least Rotation algorithm check out this [dev blog](https://www.ginkgobioworks.com/2020/04/20/fast-database-lookups-for-circular-dna-sequences/) by [Josh Hoffer](https://twitter.com/hofer) and this ridiculously hard to read [python implementation of Booth's Least Rotation on wikipedia](https://en.wikipedia.org/wiki/Lexicographically_minimal_string_rotation#Booth's_Algorithm).

For an easier to read implementation of Booth's Least Rotation you can also check out our [original implementation](https://github.com/TimothyStiles/poly/blob/346e3eb58cdd74db14eba333ba428256f77c93b0/hash.go#L40).