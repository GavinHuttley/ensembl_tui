# Installing Ensembl Data

The install step converts the downloaded data into more efficient data structures. These are written to the `install_path` as specified in the config file.

The install command requires the path to the download directory.

```
$ eti install -d <dirname>
```

!!! note
    You can utilize multiple processes on your machine for this installation step with the `-np #` argument. We recommend specifying the same number of processes as the number of genomes, e.g. `-np 10` for ten genomes.

## What Is Installed

```console exec="1" source="console" result="ansi" workdir="./docs"
$ ls data/apes-115
```

The `installed.cfg` file also specifies the Ensembl release, the software versions used during installation, and under the section `species_map`, the mapping of genome names to different "names", such as the abbreviation. This is just a plain text file which you can edit.

!!! note
    Changing the abbreviations to something that you find easier to type can be useful as these are employed in the command line interface.

```console exec="1" source="console" result="ansi" workdir="./docs"
$ cat data/apes-115/installed.cfg
```

## Check Your Installation

Once you have finished your installation, you can check its contents using the `installed` command. This includes the listing of software versions at the time of the installation (useful for troubleshooting) plus species names, abbreviations etc..

```console exec="1" source="console" result="ansi" workdir="./docs"
$ eti installed -i data/apes-115
```

!!! note
    Here we start specifying the installation directory using the `-i` option. This is required for all commands that reference an installation.

!!! warning
    At present, installation is not interruptible. If you need to reinstall, you will need to force overwriting of the current installation using the `--force_overwrite` argument.


