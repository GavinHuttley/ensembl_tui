# Quickstart

Installation of `ensembl-tui` creates a command line tool `eti` which contains a number of subcommands that allow you to acquire and then sample from Ensembl genomic datasets. The general workflow is:

In this quick start, we'll perform the following steps:

1. Create a demo config file
2. Download raw data from Ensembl
3. Make the local installation
4. Show summaries of the installation
5. Export gene meta-data
6. Export homology data

## Create a demo config

You specify the genomic resources that you want from Ensembl using a config file. `ensembl-tui` comes with an example file with comments describing what each of the components of that file are.

```
$ eti demo-config -o demo
```

> **Note**
> The config specifies download the genomes and annotations from Ensembl release 114 of  *Saccharomyces cerevisiae* and *Caenorhabditis elegans* and gene homology data. It also specifies the path to write the downloaded files and where to install them.

## Download the specified data

```
$ eti download -c demo/sample.cfg
```

> **Note**
> Downloads can be interrupted.

## Make the local installation

```
$ eti install -c ensembl_download_114
```

## Show summaries of the installation

### The top level

```
$ eti installed -i ensembl_install_114
```

### Summary of a species

```
$ eti species-summary -i ensembl_install_114 --species saccharomyces_cerevisiae
```

### Summary of compara

```
$ eti compara-summary -i ensembl_install_114
```
 This shows the relationships between the species installed.

## Export gene meta-data

```bash exec="1"
rm -rf worm # markdown-exec: hide
rm -rf yeast # markdown-exec: hide
```

```
$ eti dump-genes -i ensembl_install_114 --species saccharomyces_cerevisiae --outdir yeast
```

Show the first five lines of the output file.

```
$ head -n 5 yeast/saccharomyces_cerevisiae*.tsv
```

## Export homology data

```bash exec="1"
rm -rf worm_yeast # markdown-exec: hide
```

```
$ eti homologs -i ensembl_install_114 --ref caenorhabditis_elegans --outdir worm_yeast --homology_type ortholog_one2one --limit 5
```
Listing the files that are written into the specified `worm_yeast` directory.

```
$ ls worm_yeast
```
