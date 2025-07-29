# Subsampling whole genome alignments

The alignments function returns whole genome alignments for a given set of coordinates from a reference species.

```console exec="1" source="console" result="ansi" workdir="./docs"
$ eti alignments -i data/apes-114 --align_name "*primates*" --outdir apes_aligns --ref human --coord_names 22 --limit 5 --mask "cds,dust"
```

This command produces a directory called `apes_aligns` and individual alignments are at the top level in this directory ending in `.fa`. The mask option results in bases being replaced by the '?' for sequence segments annotated as either "cds" (exons) or "dust" (a repeat classification). This is shown for one of the alignment files produced.

```python exec="on" result="ansi" workdir="./docs" source="above"
import cogent3

loader = cogent3.get_app("load_aligned", moltype="dna")
align_dir = cogent3.open_data_store("apes_aligns", suffix="fa")
aln = loader(align_dir[1])
# pretty print the first 200 bases
print(aln[:200].to_pretty(wrap=60))
```
