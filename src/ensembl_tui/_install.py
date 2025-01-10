import shutil

from rich.progress import Progress

from ensembl_tui import _align as eti_align
from ensembl_tui import _config as eti_config


def _make_src_dest_annotation_paths(
    src_dir: pathlib.Path,
    dest_dir: pathlib.Path,
) -> list[tuple[pathlib.Path, pathlib.Path]]:
    src_dir = src_dir / "gff3"
    dest = dest_dir / elt_genome.ANNOT_STORE_NAME
    paths = list(src_dir.glob("*.gff3.gz"))
    return [(path, dest) for path in paths]


from ensembl_tui import _genome as eti_genome
from ensembl_tui import _homology as eti_homology
from ensembl_tui import _maf as eti_maf
from ensembl_tui import _species as eti_species
from ensembl_tui import _util as eti_util


def local_install_genomes(
    config: eti_config.Config,
    force_overwrite: bool,
    max_workers: int | None,
    verbose: bool = False,
    progress: Progress | None = None,
) -> None:
    if force_overwrite:
        shutil.rmtree(config.install_genomes, ignore_errors=True)
    # we create the local installation
    config.install_genomes.mkdir(parents=True, exist_ok=True)
    # we create subdirectories for each species
    for db_name in list(config.db_names):
        sp_dir = config.install_genomes / db_name
        sp_dir.mkdir(parents=True, exist_ok=True)

    # for each species, we identify the download and dest paths for annotations
    db_names = list(config.db_names)
    if max_workers:
        max_workers = min(len(db_names) + 1, max_workers)

    if verbose:
        print(f"genomes {max_workers=}")

    # we load the individual gff3 files and write to annotation db's
    src_dest_paths = []
    for db_name in config.db_names:
        src_dir = config.staging_genomes / db_name
        dest_dir = config.install_genomes / db_name
        src_dest_paths.extend(_make_src_dest_annotation_paths(src_dir, dest_dir))

    msg = "Installing features 📚"
    if progress is not None:
        writing = progress.add_task(total=len(src_dest_paths), description=msg)

    tasks = elt_util.get_iterable_tasks(
        func=elt_genome.make_annotation_db,
        series=src_dest_paths,
        max_workers=max_workers,
    )
    for db_name, prefixes in tasks:
        if verbose:
            print(f"{db_name=} {prefixes=}")

        if prefixes:
            for prefix in prefixes:
                elt_species.Species.add_stableid_prefix(db_name, prefix)

        if progress is not None:
            progress.update(writing, description=msg, advance=1)

    species_table = elt_species.Species.to_table()
    species_table.write(config.install_genomes / elt_species.SPECIES_NAME)
    species_table = eti_species.Species.to_table()
    species_table.write(config.install_genomes / eti_species.SPECIES_NAME)
    if verbose:
        print("Finished installing features ")

    msg = "Installing  🧬🧬"
    if progress is not None:
        writing = progress.add_task(total=len(db_names), description=msg, advance=0)
    # we parallelise across databases
    writer = eti_genome.fasta_to_hdf5(config=config)
    tasks = eti_util.get_iterable_tasks(
        func=writer,
        series=db_names,
        max_workers=max_workers,
    )
    for result in tasks:
        if not result:
            print(result)
            raise RuntimeError(f"{result=}")

        if progress is not None:
            progress.update(writing, description=msg, advance=1)

    if verbose:
        print("Finished installing sequences ")


def local_install_alignments(
    config: eti_config.Config,
    force_overwrite: bool,
    max_workers: int | None,
    verbose: bool = False,
    progress: Progress | None = None,
):
    if force_overwrite:
        shutil.rmtree(config.install_aligns, ignore_errors=True)

    aln_loader = eti_maf.load_align_records(set(config.db_names))

    for align_name in config.align_names:
        src_dir = config.staging_aligns / align_name
        dest_dir = config.install_aligns
        dest_dir.mkdir(parents=True, exist_ok=True)
        # write out to a db with align_name
        output_path = dest_dir / f"{align_name}.{eti_align.ALIGN_STORE_SUFFIX}"
        db = eti_align.AlignDb(source=output_path)
        records = []
        paths = list(src_dir.glob(f"{align_name}*maf*"))

        if max_workers and max_workers > 1:
            # we adjust the maximum workers to the number of paths
            max_workers = min(len(paths) + 1, max_workers or 0)

        if verbose:
            print(f"{max_workers=}")

        series = eti_util.get_iterable_tasks(
            func=aln_loader,
            series=paths,
            max_workers=max_workers,
        )

        if progress is not None:
            msg = "Installing alignments"
            writing = progress.add_task(total=len(paths), description=msg, advance=0)

        for result in series:
            if not result:
                print(result)
                raise RuntimeError

            records.extend(result)

            if progress is not None:
                progress.update(writing, description=msg, advance=1)

        db.add_records(records=records)
        db.make_indexes()
        db.close()

    if verbose:
        print("Finished installing alignments")


def local_install_homology(
    config: eti_config.Config,
    force_overwrite: bool,
    max_workers: int | None,
    verbose: bool = False,
    progress: Progress | None = None,
):
    if force_overwrite:
        shutil.rmtree(config.install_homologies, ignore_errors=True)

    config.install_homologies.mkdir(parents=True, exist_ok=True)

    outpath = config.install_homologies / eti_homology.HOMOLOGY_STORE_NAME
    db = eti_homology.HomologyDb(source=outpath)

    dirnames = []
    for sp in config.db_names:
        path = config.staging_homologies / sp
        dirnames.extend(list(path.glob("*.tsv*")))

    if max_workers:
        max_workers = min(len(dirnames) + 1, max_workers)
    else:
        max_workers = 1

    if verbose:
        print(f"homologies {max_workers=}")

    loader = eti_homology.load_homologies(allowed_species=set(config.db_names))
    if max_workers > 1:
        loader = loader + eti_homology.pickler + eti_homology.compressor

    msg = "Installing homologies"
    if progress is not None:
        writing = progress.add_task(total=len(dirnames), description=msg, advance=0)

    tasks = eti_util.get_iterable_tasks(
        func=loader,
        series=dirnames,
        max_workers=max_workers,
    )
    for result in tasks:
        if max_workers > 1:
            # reconstitute the blosc compressed data
            result = eti_homology.inflate(result)

        for rel_type, records in result.items():
            db.add_records(records=records, relationship_type=rel_type)

        if progress is not None:
            progress.update(writing, description=msg, advance=1)

    no_records = len(db) == 0
    db.make_indexes()
    db.close()
    if no_records:
        outpath.unlink()

    if verbose:
        print("Finished installing homologies")
