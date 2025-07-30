# Installing Ensembl Data

The install step converts the downloaded data into more efficient data structures. These are written to the `install_path` as specified in the config file. 

The install command requires the path to the download directory.

```
$ eti install -d <dirname>
```

> **Note**
> You can utilize multiple processes on your machine for this installation step with the `-np #` argument. We recommend specifying the same number of processes as the number of genomes, e.g. `-np 10` for ten genomes.

## Check Your Installation

Once you have finished your installation, you can check its contents using the `installed` command.

```console exec="1" source="console" result="ansi" workdir="./docs"
$ eti installed -i data/apes-114
```

> **Note**
> Here we start specifying the installation directory using the `-i` option. This is required for all commands that reference an installation.

> **Warning**
> At, present installation is not interruptible. If you need to reinstall, you will need to force overwriting of the current installation using the `--force_overwrite`` argument.
