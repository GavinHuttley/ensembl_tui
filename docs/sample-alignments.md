# Subsampling Whole Genome Alignments

The alignments function returns whole genome alignments for a given set of coordinates from a reference species. The default behavior is to focus on sampling genomic coordinates corresponding to protein coding genes from the nominated reference species. The following command samples protein coding genes from human chromosome 22 using ensembl alignments that have the word "primate" in their name (see `eti compara-summary` to list the installed). The other options are explained below.

```console exec="1" source="console" result="ansi" workdir="./docs"
$ eti alignments -i data/apes-114 --align_name "*primates*" --outdir apes_aligns --ref human --coord_names 22 --limit 5 --mask "cds,dust"
```

To sample other types of genes, use the `--ref_genes` argument, providing a csv or tsv file with a "stableid" column. To sample non-genic regions by providing explicit genomic coordinates as a delimited file with columns for species, seqid, start, stop and strand to `--ref_coords`. ([See the `eti dump-genes` command](genome.md#export-genes).)

Identifying introns is achieved using `--mask cds`. Masking follows the [cogent3](https://cogent3.org/doc/cookbook/features.html#masking-annotated-regions-on-a-sequence) approach. Nucleotide characters are replaced with the most general IUPAC ambiguity character ('?'). To select everything except an annotated region, use the `--mask_shadow` option (see cogent3 docs on [mask shadow](https://cogent3.org/doc/cookbook/features.html#how-to-get-the-shadow-of-a-feature)). To limit the masking to the reference species use the `--mask_ref` flag.

This above command command produces a directory called `apes_aligns` and individual alignments are at the top level in this directory ending in `.fa`. We show the outcome of the mask option for one of the alignment files produced.

```python exec="on" result="ansi" workdir="./docs" source="above"
import cogent3

loader = cogent3.get_app("load_aligned", moltype="dna")
align_dir = cogent3.open_data_store("apes_aligns", suffix="fa")
aln = loader(align_dir[1])
# pretty print the first 200 bases
print(aln[:200].to_pretty(wrap=60))
```
We use `cogent3` to remove the masked columns from the alignment.

```python exec="on" result="ansi" workdir="./docs" source="above"
import cogent3

loader = cogent3.get_app("load_aligned", moltype="dna")
align_dir = cogent3.open_data_store("apes_aligns", suffix="fa")
aln = loader(align_dir[0])
# remove the alignment columns containing any degenerate character
# by default this includes the gap character
aln2 = aln.no_degenerates()
```
