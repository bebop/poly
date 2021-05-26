The following is the modified version of the man pages from
[ViennaRNA](https://github.com/ViennaRNA/ViennaRNA/blob/master/man/RNAlib.texi)

Various loop parameters depend in general on the pairs closing the loops, 
as well as unpaired bases in the loops. Internally, the library
distinguishes 8 types of pairs (CG=1, GC=2, GU=3, UG=4, AU=5, UA=6, 
nonstandard=7, 0= no pair), and 5 types of bases (A=1, C=2, G=3, U=4 and 0
for anything else). Parameters belonging to pairs of type 0 are not listed
in the parameter files, but values for nonstandard pairs (type 7) and
nonstandard bases (type 0) are. Thus, a table for symmetric size 2 interior
loops would have `7*7*5*5` entries (2 pairs, two unpaired bases).

The order of entries always uses the closing pair or pairs as the first
indices followed by the unpaired bases in 5' to 3' direction.  To determine
the type of a pair consider the base at the (local) 5' end of each strand
first, i.e. use the pairs `(i, j)` and `(q, p)` for an interior
loop with `i<p<q<j` .
Consider the symmetric size 4 interior loop

```
		      5'-GAUA-3'
		      3'-CGCU-5'
```

which is equivalent to

```
					5'-UCGC-3'
					3'-AUAG-5'
```

(as the sequence is symmetric).
The first pair is GC, the second UA (not AU!), and the unpaired bases are (in 5'
to 3' direction, starting at the first pair) A U C G.
Thus the encoded base pair type array of the above sequence would be:

```
	pairs:                    GC UA  A  U  C  G
	encoded base pair array: [ 2  6  1  4  2  3]
```

The encoded base pair type array could also be described as:

```
	pairs:                    UA GC  C  G  A  U
	encoded base pair array: [ 6  2  2  3  1  4]
```

Be careful to preserve this symmetry when editing parameter tables!

### stack_energies

The list of free energies for stacked pairs, indexed by the two closing
pairs. The list should be formatted as symmetric an `7*7` matrix, conforming
to the order explained above. As an example the stacked pair

```
		      5'-GU-3'        5'-AC-3'
		      3'-CA-5'        3'-UG-5'
```

corresponds to the encoded base pair arrays

```

					GC AU         AU GC
				 [ 2  5]       [ 5  2]
```

### hairpin

Free energies of hairpin loops as a function of size. The list should
contain 31 entries on one or more lines. Since the minimum size of a
hairpin loop is 3 and we start counting with 0, the first three values
should be INF to indicate a forbidden value.

### bulge

Free energies of bulge loops. Should contain 31 entries, the first one
being INF.

### internal_loop

Free energies of internal loops. Should contain 31 entries, the first 4
being INF (since smaller loops are tabulated).

### mismatch_interior

Free energies for the interaction between the closing pair of an interior
loop and the two unpaired bases adjacent to the helix. This is a three
dimensional array indexed by the type of the closing pair (see above) and
the two unpaired bases. Since we distinguish 5 bases the list contains
`7*5*5` entries and should be formatted either as an `7*25` matrix or 7 `5*5`

matrices. The order is such that for example the mismatch

```
			       5'-CU-3'
			       3'-GC-5'
```

corresponds to entry [1 4 2] (CG=1, U=4, C=2), (in this notation
the first index runs from 1 to 7, second and third from 0 to 4)

### mismatch_hairpin

Same as above for hairpin loops.

### mismatch_enthalpies

Corresponding enthalpies for rescaling to temperatures other than 37C.

### int11_energies

Free energies for symmetric size 2 interior loops. `7*7*5*5` entries formatted
as 49 `5*5` matrices, or seven `7*25` matrices. Example:

```

			       5'-CUU-3'
			       3'-GCA-5'
```

corresponds to entry [1 5 4 2] (CG=1, AU=5, U=4, C=2), which should be
identical to [5 1 2 4] (AU=5, CG=1, C=2, U=4).

### int21_energies

Free energies for size 3 (2+1) interior loops. `7*7*5*5*5` entries formatted
in `5*5` or `5*25` matrices. The strand with a single unpaired base is listed
first, example:

```

			       5'-CU U-3'
			       3'-GCCA-5'
```

corresponds to entry [1 5 4 2 2] (CG=1, AU=5, U=4, C=2).

### int22_energies

Free energies for symmetric size 4 interior loops. To reduce the size of
parameter files this table only lists canonical bases (A, C, G, U) resulting in
a `7*7*4*4*4*4` table. See above for an example.

### dangle5

Energies for the interaction of an unpaired base on the 5' side and
adjacent to a helix in multiloops and free ends (the equivalent of mismatch
energies in interior and hairpin loops). The array is indexed by the type
of pair closing the helix and the unpaired base and, therefore, forms a `8*5`

matrix. For example the dangling base in

```
			       5'-GC-3'
			       3'- G-5'
```

corresponds to entry [1 3] (CG=1, G=3).

### dangle3

Same as above for bases on the 3' side of a helix.

```

			       5'- A-3'
			       3'-AU-5'

```

corresponds to entry [5 1] (AU=5, A=1)

### Tetraloops

Some tetraloops particularly stable tetraloops are assigned an energy
bonus. Up to forty tetraloops and their bonus energies can be listed
following the token, one sequence per line. For example:

```

       GAAA    -200

```

assigns a bonus energy of -2 kcal/mol to tetraloops containing
the sequence GAAA.
