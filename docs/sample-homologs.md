# Selecting Homologs

Homologs are related genes and the output is the raw sequence which is unaligned.

```console exec="1" source="console" result="ansi" workdir="./docs"
$ eti homologs -i data/apes-115 --outdir apes_homologs --ref human --coord_names 22 --limit 5
```

```console exec="1" source="console" result="ansi" workdir="./docs"
$ ls apes_homologs
```

By default: `"protein_coding"` genes are selected and `--homology_type` ([see compara summary](summary-compara.md#summary-compara)) is set `ortholog_one2one`. You can specify different gene biotypes ([see species summary](genome.md#summary-species)) by providing a delimited file to `--ref_genes`. This file must contain a "stableid" column where the values are the Ensembl stable IDs. To the get the full gene list for the reference species [see the `eti dump-genes` command](genome.md#export-genes).
