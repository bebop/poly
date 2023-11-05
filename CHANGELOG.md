# Changelog

All notable changes to this project will be documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.0.0/),
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## [Unreleased]

### Added
- Alternative start codons can now be used in the `synthesis/codon` DNA -> protein translation package (#305)
- Added a parser and writer for the `pileup` sequence alignment format (#329)
- Added statistics to the `synthesis/codon` package (keeping track of the observed start codon occurrences in a translation table) (#350)

### Fixed
- `fastq` parser no longer becomes de-aligned when reading (#325)
- `fastq` now handles optionals correctly (#323)

## [0.26.0] - 2023-07-22
Oops, we weren't keeping a changelog before this tag!

[unreleased]: https://github.com/TimothyStiles/poly/compare/v0.26.0...main
[0.26.0]: https://github.com/TimothyStiles/poly/releases/tag/v0.26.0