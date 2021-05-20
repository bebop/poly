## This doc explains how MFE is calcuated

### Inputs

* sequence
* structure

### Pre-processing

* `fc.sequence_encoding[1:len(sequence)]` is the encoded version of `sequence` (mapping of characters to integers can be found at `vrna_nucleotide_encode`). The encoding is A -> 1, C -> 2, G -> 3, T / U -> 4.

*
