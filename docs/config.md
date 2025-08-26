# Selecting Ensembl Data

`ensembl-tui` requires a config file to specify the data that you want to select from Ensembl.

## Create a template config file to edit

`eti` can write a template config file to a directory you specify.

```console exec="1" source="console" result="ansi" workdir="./docs"
$ eti demo-config --outpath demo
```

The config template file is written to the specified directory along with a `species.tsv` file which includes a listing of Ensembl species.

```console exec="1" source="console" result="ansi" workdir="./docs"
$ ls demo/
```

### The species contents

```console exec="1" source="console" result="ansi" workdir="./docs"
$ head -n 5 demo/species.tsv
```

### The config contents

```console exec="1" source="console" result="ansi" workdir="./docs"
$ head -n 15 demo/sample.cfg
```

## The config format

This is a `.ini` format. A section is denoted by square brackets surrounding the section name, e.g. `[remote path]`. Variables within a section are denoted by the name followed by an `=`.

### `[remote path]`

This is the section defining which Ensembl FTP server hosts the data and the remote path. At present, we are only supporting the primary Ensembl Server, so this section should be left as is.

### `[local path]`

Specify where to download the data to (`staging_path`) on your machine and where to write your installation (`install_path`).

### `[release]`

Specify the ensemble release.

### `[<species name>]`

Selecting a species is done by providing a section with the species name as a section. Available species can be found in the species.tsv file which is written out by the `demo-config` command.

For now, you must include `db=core` under the species section.

### `[compara]`

Sources from the Compara database are indicated here. Including `homologies =` (without a value) indicates that you want to get the homology information for all the selected species.

You can indicate the alignments you want as comma separated values assigned to the variable `align_names`.

!!! note
    To select the alignments that you want, you will need to navigate the Ensembl FTP site for the release you are interested in.

## Implicit selection of genomes

If you specify a whole genome alignment and do not specify any specific species, then all of the species present in that whole genome alignment will be downloaded.

For example, the following config would download 10 primate genomes along with whole genome alignments and homology data.

```ini
[remote path]
host=ftp.ensembl.org
path=pub
[local path]
staging_path=download_114
install_path=install_114
[release]
release=114
[compara]
align_names=10_primates.epo
homologies =
```

```bash exec="1" workdir="./docs"
rm -rf demo  # markdown-exec: hide
```
