---
id: library-hashing
title: Hashing Sequences
---

`poly` provides what are likely the only open source sequence hashing functions that handle circular sequences. By utilizing Booth's Least Rotation Algorithm we're able to determistically rotate circular sequences to a fixed point which makes it possible to hash them. [More info on Booth's Least Rotation here](#more-info).

Hashes make incredibly powerful unique identifiers and with a wide array of hash function options including the superfast blake3 poly has all your hashing needs covered.

## Blake3 Hashing

The golang team is currently figuring out the best way to implement blake3 into the standard lib but in the meantime `poly` provides this special function and method wrapper to hash sequences using blake3. This will eventually be deprecated in favor of only using the `GenericSequenceHash()` function and `.hash()` method wrapper.

```go
  // getting our example AnnotatedSequence struct
  puc19AnnotatedSequence := ReadJSON("data/puc19static.json")

  // there are two ways to use the blake3 Least Rotation hasher.

  // the first is with the method wrapper.
  puc19Blake3Hash := puc19AnnotatedSequence.blake3Hash()
  fmt.Println(puc19Blake3Hash)

  // the second is with the Blake3SequenceHash(annotatedSequence AnnotatedSequence) function.
  puc19Blake3Hash = puc19AnnotatedSequence.blake3Hash()
  fmt.Println(puc19Blake3Hash)
```

Again, this will be deprecated in favor of using generic hashing with blake3 in the future when available.

## Generic Hashing

`poly` also provides a generic hashing function and method wrapper for hashing sequences with arbitrary hashing functions that use the golang standard libraries' hash function interface. Check out this switch statement in the [hash command source code](https://github.com/TimothyStiles/poly/blob/f51ec1c08820394d7cab89a5a4af92d9b803f0a4/commands.go#L261) to see all that `poly` provides in the command line utility alone.

```go
  // getting our example AnnotatedSequence struct
  puc19AnnotatedSequence := ReadJSON("data/puc19static.json")

  // there are two ways to use the Least Rotation generic hasher.

  // the first is with the method wrapper where you pass your hashing function as an argument.
  puc19Sha1Hash := puc19AnnotatedSequence.hash(crypto.SHA1)
  fmt.Println(puc19Sha1Hash)

  // the second is with the GenericSequenceHash() function where you pass an AnnotatedSequence along with a hash function as arguments.
  puc19Sha1Hash = GenericSequenceHash(puc19AnnotatedSequence, crypto.SHA1)
  fmt.Println(puc19Sha1Hash)
```

## More Info
For more info about circular sequence hashing and Booth's Least Rotation algorithm check out this [dev blog](https://www.ginkgobioworks.com/2020/04/20/fast-database-lookups-for-circular-dna-sequences/) by [Josh Hoffer](https://twitter.com/hofer) and this ridiculously hard to read [python implementation of Booth's Least Rotation on wikipedia](https://en.wikipedia.org/wiki/Lexicographically_minimal_string_rotation#Booth's_Algorithm).

For an easier to read implementation of Booth's Least Rotation you can also check out our [original implementation](https://github.com/TimothyStiles/poly/blob/346e3eb58cdd74db14eba333ba428256f77c93b0/hash.go#L40).