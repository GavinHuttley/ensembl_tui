import io
import pathlib
import re
import shutil
import typing

import cogent3
from rich.progress import Progress

from ensembl_tui import _config as eti_config
from ensembl_tui import _ftp_download as eti_ftp
from ensembl_tui import _ingest_annotation as eti_db_ingest
from ensembl_tui import _mysql_core_attr as eti_db_attr
from ensembl_tui import _name as eti_name
from ensembl_tui import _site_map as eti_site_map
from ensembl_tui import _species as eti_species
from ensembl_tui import _util as eti_util

if typing.TYPE_CHECKING:  # pragma: no cover
    from cogent3.core.table import Table
    from cogent3.core.tree import PhyloNode

DEFAULT_CFG = eti_util.get_resource_path("sample.cfg")

_valid_seq = re.compile(r"dna[.](nonchromosomal|toplevel)\.fa\.gz")


def valid_seq_file(name: str) -> bool:
    """unmasked genomic DNA sequences"""
    return _valid_seq.search(name) is not None


def _remove_tmpdirs(path: eti_util.PathType) -> None:
    """delete any tmp dirs left over from unsuccessful runs"""
    tmpdirs = [p for p in path.glob("tmp*") if p.is_dir()]
    for tmpdir in tmpdirs:
        shutil.rmtree(tmpdir)


def get_core_db_dirnames(
    config: eti_config.Config, site_map: eti_site_map.SiteMap
) -> dict[str, str]:
    """maps species name to ftp path to mysql core dbs"""
    remote_release_path = site_map.get_remote_release_path(config.release)
    # get all the mysql db names
    all_db_names = list(
        eti_ftp.listdir(site_map.site, f"{remote_release_path}/mysql"),
    )
    selected_species = {}
    core_db_names = set(config.get_core_db_names())
    for db_name in all_db_names:
        if "_core_" not in db_name:
            continue
        db = eti_name.EnsemblDbName(db_name.rsplit("/", maxsplit=1)[1])
        if db.db_type == "core" and db.prefix in core_db_names:
            selected_species[db.prefix] = db_name
    return selected_species


def get_remote_mysql_paths(db_name: str) -> list[str]:
    return [f"{db_name}/{name}" for name in eti_db_attr.make_mysqldump_names()]


def make_core_db_templates(
    *,
    config: eti_config.Config,
    site_map: eti_site_map.SiteMap,
    sp_db_map: dict[str, str],
    progress: Progress | None = None,
) -> None:
    """creates duckdb db files for importing Ensembl mysql data

    Parameters
    ----------
    config
        eti configuration
    sp_db_map
        mapping of species to mysql db names
    progress
        rich.progress context manager for tracking progress
    site_map
        site map stores attributes key attributes for the Ensembl site

    Notes
    -----
    Communicates with the Ensembl MySQL server to infer the table schema's.
    """
    table_names = eti_db_attr.get_all_tables()
    if progress is not None:
        msg = "Making db templates"
        make_templates = progress.add_task(
            total=len(table_names),
            description=msg,
        )

    template_dest = config.staging_template_path
    # get one species db name which wqe use to infer the db schema
    db_name = next(iter(sp_db_map.values())).split("/")[-1]
    template_dest.mkdir(parents=True, exist_ok=True)
    for table_name in table_names:
        eti_db_ingest.make_table_template(
            dest_dir=template_dest,
            db_name=db_name,
            table_name=table_name,
            db_host=site_map.db_host,
            db_port=site_map.db_port,
        )
        if progress is not None:
            progress.update(make_templates, description=msg, advance=1)


def download_species(
    *,
    site_map: eti_site_map.SiteMap,
    config: eti_config.Config,
    debug: bool,
    verbose: bool,
    progress: Progress | None = None,
) -> None:
    """download seq and annotation data"""
    remote_release_path = site_map.get_remote_release_path(config.release)
    remote_template = f"{remote_release_path}/" + "{}"
    if verbose:
        eti_util.print_colour(
            text=f"DOWNLOADING\n  ensembl release={config.release}",
            colour="green",
        )
        eti_util.print_colour(
            text="\n".join(f"  {d}" for d in config.species_dbs),
            colour="green",
        )
        eti_util.print_colour(
            text=f"\nWRITING to output path={config.staging_genomes}\n",
            colour="green",
        )

    sp_db_map = get_core_db_dirnames(config, site_map=site_map)

    # create the duckdb templates for the tables, if they don't exist
    make_core_db_templates(
        config=config,
        sp_db_map=sp_db_map,
        site_map=site_map,
        progress=progress,
    )

    msg = "Downloading genomes"
    if progress is not None:
        species_download = progress.add_task(
            total=len(config.species_dbs),
            description=msg,
        )

    for genome_name in config.species_dbs:
        abbrev = config.species_map.get_abbreviation(genome_name)
        db_prefix = config.species_map.get_ensembl_db_prefix(genome_name)
        if genome_name != db_prefix:
            # genome in a collection
            collection_name = db_prefix
            local_root = config.staging_genomes / genome_name
        else:
            collection_name = None
            local_root = config.staging_genomes / db_prefix

        local_root.mkdir(parents=True, exist_ok=True)

        # getting genome sequences
        remote = site_map.get_seqs_path(genome_name, collection_name=collection_name)
        remote_dir = remote_template.format(remote)
        remote_paths = list(
            eti_ftp.listdir(site_map.site, path=remote_dir, pattern=valid_seq_file),
        )
        if verbose:
            eti_util.print_colour(text=f"{remote_paths=}", colour="yellow")

        if debug:
            # we need the checksum files
            paths = [p for p in remote_paths if eti_util.is_signature(p)]
            # but fewer data files, to reduce time for debugging
            remote_paths = [p for p in remote_paths if not eti_util.dont_checksum(p)]
            remote_paths = remote_paths[:4] + paths

        dest_path = local_root / "fasta"
        dest_path.mkdir(parents=True, exist_ok=True)
        # cleanup previous download attempts
        _remove_tmpdirs(dest_path)
        icon = "🧬🧬"
        eti_ftp.download_data(
            host=site_map.site,
            local_dest=dest_path,
            remote_paths=remote_paths,
            description=f"{abbrev} {icon}",
            do_checksum=True,
            progress=progress,
        )

        # getting the annotations from mysql tables
        remote_dir = sp_db_map[db_prefix]
        remote_paths = get_remote_mysql_paths(remote_dir)
        # the mysql data will always be under the db_prefix,
        # even if it's a collection
        dest_path = config.staging_genomes / db_prefix / "mysql"
        dest_path.mkdir(parents=True, exist_ok=True)
        # cleanup previous download attempts
        _remove_tmpdirs(dest_path)
        icon = "📚"
        eti_ftp.download_data(
            host=site_map.site,
            local_dest=dest_path,
            remote_paths=remote_paths,
            description=f"{abbrev} {icon}",
            do_checksum=True,
            progress=progress,
        )

        if progress is not None:
            progress.update(species_download, description=msg, advance=1)


class valid_compara_align:  # noqa: N801
    """whole genome alignment data"""

    def __init__(self) -> None:
        self._valid = re.compile("([.](emf|maf)[.]gz|README|MD5SUM)")

    def __call__(self, name: str) -> bool:
        return self._valid.search(name) is not None


def download_aligns(
    *,
    config: eti_config.Config,
    site_map: eti_site_map.SiteMap,
    debug: bool,
    verbose: bool,
    progress: Progress | None = None,
) -> None:
    """download whole genome alignments"""
    if not config.align_names:
        return

    remote_template = f"{site_map.remote_path}/release-{config.release}/{site_map.alignments_path}/{{}}"

    msg = "Downloading alignments"
    if progress is not None:
        align_download = progress.add_task(
            total=len(config.species_dbs),
            description=msg,
        )

    valid_compara = valid_compara_align()
    for align_name in config.align_names:
        remote_path = remote_template.format(align_name)
        remote_paths = list(eti_ftp.listdir(site_map.site, remote_path, valid_compara))
        if verbose:
            print(remote_paths)

        if debug:
            # we need the checksum files
            paths = [p for p in remote_paths if eti_util.is_signature(p)]
            remote_paths = [p for p in remote_paths if not eti_util.is_signature(p)]
            remote_paths = remote_paths[:4] + paths

        local_dir = config.staging_aligns / align_name
        local_dir.mkdir(parents=True, exist_ok=True)
        _remove_tmpdirs(local_dir)
        eti_ftp.download_data(
            host=site_map.site,
            local_dest=local_dir,
            remote_paths=remote_paths,
            description=f"{align_name[:10]}...",
            do_checksum=True,
            progress=progress,
        )

        if progress is not None:
            progress.update(align_download, description=msg, advance=1)

    return


class valid_compara_homology:  # noqa: N801
    """homology tsv files"""

    def __init__(self) -> None:
        self._valid = re.compile("([.]tsv|[.]tsv[.]gz|README|MD5SUM)$")

    def __call__(self, name: str) -> bool:
        return self._valid.search(name) is not None


def download_homology(
    *,
    config: eti_config.Config,
    debug: bool,
    verbose: bool,
    site_map: eti_site_map.SiteMap,
    progress: Progress | None = None,
) -> None:
    """downloads tsv homology files for each genome"""
    if not config.homologies:
        return

    # change homologies path to take an argument, which modifies order of path/genome
    remote_template = f"{site_map.remote_path}/release-{config.release}/{site_map.homologies_path}/{{}}"

    local = config.staging_homologies

    msg = "Downloading homology"
    if progress is not None:
        species_download = progress.add_task(
            total=len(config.species_dbs),
            description=msg,
        )

    for genome_name in config.species_dbs:
        abbrev = config.species_map.get_abbreviation(genome_name)
        db_name = config.species_map.get_ensembl_db_prefix(genome_name)
        if db_name == genome_name:
            remote_path = remote_template.format(genome_name)
        else:
            # genome is in a collection
            remote_path = remote_template.format(f"{db_name}/{genome_name}")

        remote_paths = list(
            eti_ftp.listdir(site_map.site, remote_path, valid_compara_homology()),
        )
        if verbose:
            print(f"{remote_path=}", f"{remote_paths=}", sep="\n")

        if debug:
            # we need the checksum files
            remote_paths = [p for p in remote_paths if not eti_util.is_signature(p)]
            remote_paths = remote_paths[:4]

        local_dir = local / genome_name
        local_dir.mkdir(parents=True, exist_ok=True)
        _remove_tmpdirs(local_dir)
        eti_ftp.download_data(
            host=site_map.site,
            local_dest=local_dir,
            remote_paths=remote_paths,
            description=f"{abbrev}",
            do_checksum=False,  # no checksums for species homology files
            progress=progress,
        )

        if progress is not None:
            progress.update(species_download, description=msg, advance=1)

    return


def download_ensembl_tree(
    *,
    host: str,
    site_map: eti_site_map.SiteMap,
    release: str,
    tree_fname: str,
) -> typing.Optional["PhyloNode"]:
    """loads a tree from Ensembl"""
    if site_map.trees_path is None:
        return None
    url = f"https://{host}/{site_map.remote_path}/release-{release}/{site_map.trees_path}/{tree_fname}"
    return cogent3.load_tree(url, underscore_unmunge=False)


def get_ensembl_trees(
    *,
    host: str,
    release: str,
    site_map: eti_site_map.SiteMap,
) -> list[str]:
    """returns trees from ensembl compara"""
    if site_map.trees_path is None:
        return []

    path = f"{site_map.remote_path}/release-{release}/{site_map.trees_path}"
    return list(
        eti_ftp.listdir(host=host, path=path, pattern=lambda x: x.endswith(".nh")),
    )


def get_species_for_alignments(
    *,
    host: str,
    release: str,
    align_names: typing.Iterable[str],
    site_map: eti_site_map.SiteMap,
    species_map: eti_species.SpeciesNameMap,
) -> dict[str, list[str]]:
    """return the species for the indicated alignments"""
    ensembl_trees = get_ensembl_trees(
        site_map=site_map,
        host=host,
        release=release,
    )
    if not ensembl_trees:
        return {}

    aligns_trees = eti_util.trees_for_aligns(align_names, ensembl_trees)
    species = {}
    for tree_path in aligns_trees.values():
        tree = download_ensembl_tree(
            site_map=site_map,
            host=host,
            release=release,
            tree_fname=pathlib.Path(tree_path).name,
        )
        if tree is None:
            continue
        # dict structure is {common name: db prefix}, just use common name
        species |= {
            n: ["core"]
            for n in eti_species.species_from_ensembl_tree(
                tree, species_map=species_map
            )
        }
    return species


def download_species_table(
    *,
    site_map: eti_site_map.SiteMap,
) -> "Table":
    """downloads the species file for the given Ensembl division"""
    remote = f"{site_map.remote_path}/current/{site_map.species_file_name}"
    ftp = eti_ftp.configured_ftp(host=site_map.site)
    buff = io.BytesIO()
    ftp.retrbinary(f"RETR {remote}", buff.write)
    buff.seek(0)
    data = buff.getvalue().decode("utf-8").splitlines()
    header = data.pop(0).split("\t")
    header[0] = header[0].lstrip("#")
    # the file is tab delimited but does not have a consistent number of columns
    num_col = len(header)
    rows = [row.split("\t")[:num_col] for row in data]
    table = cogent3.make_table(header=header, data=rows)
    abbrevs = eti_species.make_unique_abbrevs(table.columns["species"])
    table = table.with_new_column("abbrev", lambda x: abbrevs[x], columns=["species"])
    old = ["name", "species", "core_db"]
    table = table.with_new_header(old, ["common_name", "genome_name", "db_prefix"])
    old = ["abbrev", *old]
    return table.get_columns(["abbrev"] + [c for c in table.header if c not in old])
