# Subsampling whole genome alignments

The alignments function returns whole genome alignments for a given set of coordinates from a reference species.

```console exec="1" source="console" result="ansi" workdir="./docs"
$ eti alignments -i data/apes-114 --align_name "*primates*" --outdir apes_aligns --ref human --coord_names 22 --limit 5 --mask "cds,dust"
```
