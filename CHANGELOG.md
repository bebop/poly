# Changelog

All notable changes to this project will be documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.0.0/),
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## [Unreleased]

### Fixed
 - Made it possible to simulate primers shorter than design minimum.

## [0.31.1] - 2024-01-31

### Added
- Fixed package level doc strings for search and bwt packages.

[0.31.1]: https://github.com/TimothyStiles/poly/releases/tag/v0.31.0

## [0.31.0] - 2024-01-31

### Added
- Basic BWT for sub-sequence count and offset for sequence alignment. Only supports exact matches for now.
- Moved `BWT`, `align`, and `mash` packages to new `search` sub-directory.
- Implemented Run-Length Burrows Wheeler Transform.

[0.31.0]: https://github.com/TimothyStiles/poly/releases/tag/v0.31.0


## [0.30.0] - 2023-12-18
Oops, we weren't keeping a changelog before this tag!

### Fixed
-  Fixed bug that produced wrong overhang in linear, non-directional, single cut reactions. #408 

[0.30.0]: https://github.com/TimothyStiles/poly/releases/tag/v0.30.0
