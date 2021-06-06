---
id: cli-hashing
title: Hashing Sequences
---

`poly` provides what is likely the only open source sequence hashing tool that handles circular sequences. By utilizing Booth's Least Rotation algorithm we're able to determistically rotate circular sequences to a fixed point which makes it possible to hash them. [More info on Booth's Least Rotation here](#more-info).

Hashes make incredibly powerful unique identifiers and with `poly hash` using our super fast seqhashV1 algorithm you can create them faster than ever before.

## Hashing a sequence from file

To hash a sequence from file all you have to do is this:

``` bash
poly hash bsub.gbk
```

`poly hash` will then parse the sequence string, rotate it to a deterministic point if it's circular, then hash it using seqhashV1 returning something like this:

``` bash
v1_DCD_e5bcb24c256aae63e4f22c948e06a38590dc9722de454b895f47dd19e40d5002
```

## Hashing a sequence from a stream

`poly hash` can also accept input from streams to create a hash:

```bash
cat bsub.gbk | poly ha -i gbk
```

## Hashing multiple file inputs and writing strings to stdout

`poly hash` can also take all files in a directory and spit out their hashes along with their original file paths.

```bash
poly ha *.gbk *.gb *.gff
```

Result:

```bash
v1_DCD_4b0616d1b3fc632e42d78521deb38b44fba95cca9fde159e01cd567fa996ceb9  puc19.gbk
v1_DLD_11e56bd6f159732e116b65f074768998ffea14688b45686b2fd71e5d13489d2d  sample.gbk
v1_DCD_e5bcb24c256aae63e4f22c948e06a38590dc9722de454b895f47dd19e40d5002  bsub.gbk
```

## More Info

For more info about `poly hash` and its usage try running `poly help hash` or `poly help ha` from your command line.

For more info about circular sequence hashing and Booth's Least Rotation algorithm check out this [dev blog](https://www.ginkgobioworks.com/2020/04/20/fast-database-lookups-for-circular-dna-sequences/) by [Josh Hoffer](https://twitter.com/hofer) and this ridiculously hard to read [python implementation of Booth's Least Rotation on wikipedia](https://en.wikipedia.org/wiki/Lexicographically_minimal_string_rotation#Booth's_Algorithm).

For an easier to read implementation of Booth's Least Rotation you can also check out our [original implementation](https://github.com/TimothyStiles/poly/blob/346e3eb58cdd74db14eba333ba428256f77c93b0/hash.go#L40).