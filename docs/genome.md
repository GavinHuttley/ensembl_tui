# Querying genomes

## Summary of this installation

```console exec="1" source="console" result="ansi" workdir="./docs"
$ eti installed -i data/apes-114
```
> **Note**
> :material-download: [Download all the data](ensembl_tui_data.zip) (zip, ~196 MB).

## Summary for a species

```console exec="1" source="console" result="ansi" workdir="./docs"
$ eti species-summary -i data/apes-114 --species human
```

## Export gene meta-data for a species

> **Note**
> The list of data from this query only covers human chromosome 22 because we are using a custom subset of the original Ensembl data.

```console exec="1" source="console" result="ansi" workdir="./docs"
$ eti dump-genes -i data/apes-114 --species human -od human_data
```

```console exec="1" source="console" result="ansi" workdir="./docs"
$ head human_data/homo_sapiens-114-gene_metadata.tsv
```

```bash exec="1"
rm -rf human_data # markdown-exec: hide
```
