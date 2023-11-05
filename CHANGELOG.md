# Changelog

All notable changes to this project will be documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.0.0/),
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## [Unreleased]

### Added
- Alternative start codons can now be used in the `synthesis/codon` DNA -> protein translation package (#305)
- Added a parser and writer for the `pileup` sequence alignment format (#329)
- Created copy methods for Feature and Location to address concerns raised by [(#342)](https://github.com/TimothyStiles/poly/issues/342)
- Created new methods to convert polyjson -> genbank.
- Created new `Feature.StoreSequence` method to enable [(#388)](https://github.com/TimothyStiles/poly/issues/388)

### Changed
- **Breaking**: Genbank parser uses new custom multimap for `Feature.Attributes`, which allows for duplicate keys. This changes the type of Features.Attributes from `map[string]string` to `MultiMap[string, string]`, an alias for `map[string]string` defined in `multimap.go`. [(#383)](https://github.com/TimothyStiles/poly/issues/383)
- Improves error reporting for genbank parse errors via a new `ParseError` struct.

### Fixed
- `fastq` parser no longer becomes de-aligned when reading (#325)
- `fastq` now handles optionals correctly (#323)
- Adds functional test and fix for [(#313)](https://github.com/TimothyStiles/poly/issues/313).
- In addition to expanding the set of genbank files which can be validly parsed, the parser is more vocal when it encounters unusual syntax in the "feature" section. This "fail fast" approach is better as there were cases where inputs triggered a codepath which would neither return a valid Genbank object nor an error, and should help with debugging.

## [0.26.0] - 2023-07-22
Oops, we weren't keeping a changelog before this tag!

[unreleased]: https://github.com/TimothyStiles/poly/compare/v0.26.0...main
[0.26.0]: https://github.com/TimothyStiles/poly/releases/tag/v0.26.0