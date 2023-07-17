# dce-extractor

This is a small package for extracting contextual information about dual-coding exons
highlighted by Mirage2 during intra-species multiple sequence alignment.

## Usage

When you first download this package, you'll need to run the pseudo-makefile
`lazy_maker.sh` within the `hsi` subdirectory.  Then, you can use the following
to run the extractor:

    ./DCE-Extractor.pl [Mirage2-Results] [Species-Guide]

Where `Mirage2-Results` is the output directory produced by a run of Mirage2, and
`Species-Guide` is the guide file used to associate species with their genomes
during that run of Mirage2.

## Output

The results of the extractor will be placed in a directory named `Extracted-DCEs/`
Within this directory, there will be a CSV-formatted file for each species that
lists each collection of left- and right-hand genomic windows of 32 nucleotides
(centered around the splice-in and splice-out sites that define the dual-coding
region), organized by "groups" within gene families.

In addition to the species-specific CSV files, there is also a subdirectory for
each species where the detailed protein-to-genome alignments are depicted.

Finally, there may be a file called `Gappy-Sequences.out` which lists any sequences
for which the mapping of amino acids to the genomic region relevant to a dual-coding
region highlighted by Mirage2 is gappy (suggestive of low-quality alignment and less
likely to be genuinely dual-coding).

