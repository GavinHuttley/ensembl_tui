# Quickstart

Installation of `ensembl-tui` creates a command line tool `eti` which contains a number of subcommands that allow you to acquire and then sample from Ensembl genomic datasets. The general workflow is:

1. Create a demo config file
1. Download raw data from Ensembl
1. Make the local installation
1. Show summaries of the installation
1. Export gene meta-data
1. Export homology data

## Create A Demo Config

You specify the genomic resources that you want from Ensembl using a config file. `ensembl-tui` comes with an example file with comments describing what each of the components of that file are.

```console exec="1" source="console" result="ansi" workdir="./docs"
$ eti demo-config -o demo
```

!!! note
    The config specifies download the genomes and annotations from Ensembl release 115 of *Saccharomyces cerevisiae* and *Caenorhabditis elegans* and gene homology data. It also specifies the path to write the downloaded files and where to install them.

!!! warning
    Edit this file before using! It also specifies primate whole genome alignments -- which are large!

## Download The Specified Data

We use a custom config file which specifies just "yeast" and "worm". (You can do this yourself by :material-download: [downloading the small.cfg](data/small.cfg) and executing the following command

```
$ eti download -c <path to>/small.cfg
```

The data will be downloaded to `staging_path` specified in `small.cfg`, which is interpreted relative to the directory in which you executed the command.

!!! note
    If a download is interrupted and restarted, `eti` resumes downloads from where it stopped.

## Make The Local Installation

```console exec="1" source="console" result="ansi" workdir="./docs"
$ eti install -d data/small-download
```

## Show Summaries Of The Installation

### The Top Level

```console exec="1" source="console" result="ansi" workdir="./docs"
$ eti installed -i data/small-install
```

### Summary Of A Species

```console exec="1" source="console" result="ansi" workdir="./docs"
$ eti species-summary -i data/small-install --species saccharomyces_cerevisiae
```

### Summary Of Compara

```console exec="1" source="console" result="ansi" workdir="./docs"
$ eti compara-summary -i data/small-install
```

This shows the relationships between the species installed.

## Export Gene meta-data

```console exec="1" source="console" result="ansi" workdir="./docs"
$ eti dump-genes -i data/small-install --species saccharomyces_cerevisiae --outdir yeast
```

Show the first five lines of the output file.

```console exec="1" source="console" result="ansi" workdir="./docs"
$ head -n 5 yeast/saccharomyces_cerevisiae*.tsv
```

## Export Homology Data

```console exec="1" source="console" result="ansi" workdir="./docs"
$ eti homologs -i data/small-install --ref caenorhabditis_elegans --outdir worm_yeast --homology_type ortholog_one2one --limit 5
```

Listing the files that are written into the specified `worm_yeast` directory.

```console exec="1" source="console" result="ansi" workdir="./docs"
$ ls worm_yeast
```

!!! note
    The `not_completed` directory will contain any errors that occurred during the homolog command.
