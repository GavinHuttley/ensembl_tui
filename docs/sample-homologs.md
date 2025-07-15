# Selecting homologs

Homologs are related genes and the output is the raw sequence which is unaligned.

```console exec="1" source="console" result="ansi" workdir="./docs"
$ eti homologs -i data/apes-114 --outdir apes_homologs --ref human --coord_names 22 --limit 5
```

```console exec="1" source="console" result="ansi" workdir="./docs"
$ ls apes_homologs
```
