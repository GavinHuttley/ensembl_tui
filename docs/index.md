`ensembl-tui` provides a textual interface to localised ensembl genomic data.

## Installing The Software

`ensembl-tui` can be installed from PyPI as follows

```
pip install ensembl-tui
```

!!! note
    If you experience any errors during installation, we recommend using [uv pip](https://docs.astral.sh/uv/). This command provides much better error messages than the standard `pip` command. If you cannot resolve the installation problem, please [open an issue](https://github.com/cogent3/ensembl_tui/issues).

!!! note
    The first usage of `ensembl-tui` is slow because both `cogent3` and `ensembl-tui` are compiling some functions. This is a one-time cost.

## Installation And Usage With `uv`

Speaking of `uv`, it provides a simplified approach to install `eti` as a command-line only tool as

```
uv tool install ensembl-tui
```

For the examples in this documentation, you can access the `eti` command as

```
uvx --from ensembl-tui eti
```

without having to activate a virtual environment.

## Getting Help

To see the options for `eti` or one of its subcommands, just enter the expression in the terminal and press return. For example,

```console exec="1" source="console" result="ansi" workdir="./docs" returncode="2"
$ eti
```

lists all of the subcommands. While

```console exec="1" source="console" result="ansi" workdir="./docs" returncode="2"
$ eti download
```

shows the options for the download command.
