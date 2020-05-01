# (poly)merase
### Transcribing between annotated genomic sequence file formats.


# Sequence objects

GFF <-> seq obj
GBK <-> seq obj


## struct types

* META
* ANNOTATIONS
* SEQUENCE


GFF3 Meta:

1. gff-version
2. sequence region
3. region start
4. region stop


Genbank meta (a ton of shit)~~~~
1. locus
2. size
3. type (DNA, RNA, ETC)
4. genbank division
5. date
6. definition
7. accession
8. version
9. keywords
10. source
* organism
7. Reference
* Authors
* title
* journal
* pubmed

GFF3 annotation :

1. seqid 
2. source - algo description?
3. type - CDS, gene, exon, etc
4. start - index of nucleotide location
5. end - index of nucleotide location
6. score - no score == "."
7. strand - + positive, - negative, ? unkown
8. phase - where feature begins with reference to the reading frame. 0, 1, 2
9. attributes - all the extra unstructure stuff as list of strings.
  
Genbank feature annotation struct:
1. source
2. 
## CLI

### Installation

## Library