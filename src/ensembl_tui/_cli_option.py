import pathlib
import sys

import click
from click import Context, Option
from cogent3 import load_table

from ensembl_tui import _config as eti_config
from ensembl_tui import _genome as eti_genome
from ensembl_tui import _species as eti_species
from ensembl_tui import _util as eti_util


def stableids_from_tsv(
    ctx: "Context",
    param: "Option",
    tsv_file: pathlib.Path,
) -> list[str] | None:
    if not tsv_file:
        return None
    table = load_table(tsv_file)
    if "stableid" not in table.columns:
        eti_util.print_colour(
            text=f"'stableid' column missing from {str(tsv_file)!r}",
            colour="red",
        )
        sys.exit(1)
    return table.columns["stableid"].tolist()


def values_from_csv_or_file(
    ctx: "Context",  # noqa: ARG001
    param: "Option",  # noqa: ARG001
    value: str | None,
) -> list[str] | None:
    """extract values from command line or a file

    Notes
    -----
    converts either comma separated values or a file with one value per line
    into values
    """
    if not value:
        return None

    path = pathlib.Path(value)
    if path.is_file():
        return [l.strip() for l in path.read_text().splitlines()]

    return [f.strip() for f in value.split(",")]


def get_installed_config_path(
    ctx: "Context",  # noqa: ARG001
    param: "Option",  # noqa: ARG001
    path: pathlib.Path | str | None,
) -> pathlib.Path:
    """path to installed.cfg"""
    path = pathlib.Path(path or ".")
    if path.name == eti_config.INSTALLED_CONFIG_NAME:
        return path

    path = path / eti_config.INSTALLED_CONFIG_NAME
    if not path.exists():
        eti_util.print_colour(text=f"{path!s} missing", colour="red")
        sys.exit(1)
    return path


def species_names_from_csv(
    ctx: "Context",
    param: "Option",
    species: str,
) -> list[str] | None:
    """returns species names"""
    species_names = values_from_csv_or_file(ctx, param, species)
    species_names = None if species_names == [""] else species_names
    if species_names is None:
        return None

    db_names = []
    for name in species_names:
        try:
            db_name = eti_species.Species.get_ensembl_db_prefix(name)
        except ValueError:
            eti_util.print_colour(text=f"ERROR: unknown species {name!r}", colour="red")
            sys.exit(1)

        db_names.append(db_name)

    return db_names


def genome_coords_from_tsv(
    ctx: "Context",  # noqa: ARG001
    param: "Option",  # noqa: ARG001
    tsv_file: pathlib.Path | None,
) -> list[eti_genome.genome_segment]:
    """reads a tsv file containing genomic coordinates and converts each
    line into a genome segment instance

    Notes
    -----
    A tsv file with the following column headings
    species seqid start stop strand
    """
    if not tsv_file:
        return None

    if not tsv_file.exists():
        eti_util.print_colour(
            text=f"ERROR: file {str(tsv_file)!r} does not exist",
            colour="red",
        )
        sys.exit(1)

    try:
        table = load_table(tsv_file, sep="\t")
    except Exception as e:  # noqa: BLE001
        eti_util.print_colour(
            text=f"ERROR: failed to load file '{tsv_file}'\n{e}",
            colour="red",
        )
        sys.exit(1)

    if table.shape[1] != 5:
        msg = f"ERROR: genome coord tsv must have 5 columns, got {table.shape[1]}"
        eti_util.print_colour(
            text=msg,
            colour="red",
        )
        sys.exit(1)

    columns = "species", "seqid", "start", "stop", "strand"
    required_columns = set(columns)
    header = set(table.header)
    if (header & required_columns) != required_columns:
        eti_util.print_colour(
            text=f"ERROR: genome coord tsv missing required columns {required_columns - header}",
            colour="red",
        )
        sys.exit(1)

    for col in ["start", "stop"]:
        if not table.columns[col].dtype.name.startswith("int"):
            eti_util.print_colour(
                text=f"ERROR: all values of {col!r} must be integers",
                colour="red",
            )
            sys.exit(1)

    return [
        eti_genome.genome_segment(
            species=str(eti_species.Species.get_ensembl_db_prefix(sp)),
            seqid=str(seqid),
            start=int(start),
            stop=int(stop),
            strand=str(strand),
        )
        for sp, seqid, start, stop, strand in table.to_list(columns=columns)
    ]


csv_or_file_help = "(comma separated or a path to file of names, one per line)"

species = click.option(
    "--species",
    required=True,
    callback=species_names_from_csv,
    help="Single species name or multiple (comma separated).",
)
mask = click.option(
    "--mask",
    callback=values_from_csv_or_file,
    help=f"mask the specified biotypes {csv_or_file_help}.",
)
mask_shadow = click.option(
    "--mask_shadow",
    callback=values_from_csv_or_file,
    help=f"mask everything but the specified biotypes {csv_or_file_help}.",
)
coord_names = click.option(
    "--coord_names",
    default=None,
    callback=values_from_csv_or_file,
    help=f"list of ref species chrom/coord names {csv_or_file_help}.",
)
cfgpath = click.option(
    "-c",
    "--configpath",
    type=pathlib.Path,
    help="Path to config file specifying databases, (only "
    "species or compara at present).",
)
download = click.option(
    "-d",
    "--download",
    type=pathlib.Path,
    help="Path to local download directory containing a cfg file.",
)
installed = click.option(
    "-i",
    "--installed",
    required=True,
    callback=get_installed_config_path,
    help="Path to root directory of an installation.",
)
outdir = click.option(
    "-od",
    "--outdir",
    required=True,
    type=pathlib.Path,
    help="Path to write files",
)
align_name = click.option(
    "--align_name",
    default=None,
    required=True,
    help="Ensembl alignment name or a glob pattern, e.g. '*primates*'.",
)
ref = click.option("--ref", default=None, help="Reference species.")
ref_genes = click.option(
    "--ref_genes",
    default=None,
    type=pathlib.Path,
    callback=stableids_from_tsv,
    help=".csv or .tsv file with a header containing a stableid column.",
)
ref_coords = click.option(
    "--ref_coords",
    type=pathlib.Path,
    default=None,
    callback=genome_coords_from_tsv,
    help="Path to tsv file with genomic coordinates in ref species.",
)
mask_ref = click.option(
    "--mask_ref",
    is_flag=True,
    help="Masking uses features from ref species only.",
)
limit = click.option(
    "--limit",
    type=int,
    default=None,
    help="Limit to this number of genes.",
    show_default=True,
)
verbose = click.option(
    "-v",
    "--verbose",
    is_flag=True,
)
force = click.option(
    "-f",
    "--force_overwrite",
    is_flag=True,
    help="Overwrite existing data.",
)
debug = click.option(
    "-d",
    "--debug",
    is_flag=True,
    help="Maximum verbosity, and reduces number of downloads, etc...",
)
dbrc_out = click.option(
    "-o",
    "--outpath",
    type=pathlib.Path,
    help="Path to directory to export all rc contents.",
)
nprocs = click.option(
    "-np",
    "--num_procs",
    type=int,
    default=1,
    help="Number of procs to use.",
    show_default=True,
)
