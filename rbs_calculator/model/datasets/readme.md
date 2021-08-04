This folder contains the datasets used to create our RBS Calculator model.

`source` contains the 1014IC and Flow-Seq datasets. It is not used, but is
included incase we'd like to regenerate the `train` dataset.

`train` - To create an accurate model for the RBS calculator, we need high
quality data which includes the delimitation of the 5' UTR and CDS regions
of the mRNA sequences. 67 sequences from the 1014IC dataset didn't contain
this information so they were removed. Apart from this, `train` is a copy of
`source`.

`salis_v2_1.xlsx` contians the predicted output of the sequences in the `train`
dataset from Salis' RBS Calculator v2.1. This dataset is included to aid in the
initial development of our RBS Calculator. The excel file also includes a
sheet to compare the output from our RBS calculator with the results from
the RBS Calculator v2.1.

1014IC:
This dataset consists of 1014 individually characterized mRNA sequences. This
dataset has been taken from Supplementary Data S1 (sb0c00394_si_002.xlsx)
available at:
https://pubs.acs.org/doi/10.1021/acssynbio.0c00394
