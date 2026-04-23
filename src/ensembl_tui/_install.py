import shutil

from scinexus.progress import Progress

from ensembl_tui import _config as eti_config
from ensembl_tui import _genome as eti_genome
from ensembl_tui import _ingest_align as ingest_aln
from ensembl_tui import _ingest_annotation as eti_db_ingest
from ensembl_tui import _ingest_homology as homology_ingest
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
    db_names = [config.species_map.get_genome_name(sp) for sp in config.species_dbs]
    for db_name in db_names:
        sp_dir = config.install_genomes / db_name
        sp_dir.mkdir(parents=True, exist_ok=True)

    # for each species, we identify the download and dest paths for annotations
    if max_workers:
        max_workers = min(len(db_names) + 1, max_workers)

    if verbose:
        eti_util.print_colour(f"\nInstalling genomes {max_workers=}", "yellow")

    # we do this the installation of features in serial for now
    writer = eti_db_ingest.mysql_dump_to_parquet(config=config)
    tasks = eti_util.get_iterable_tasks(
        func=writer,
        series=db_names,
        max_workers=max_workers,
    )
    pbar = progress.child(leave=True) if progress is not None else progress
    task_iter = (
        pbar(tasks, total=len(db_names), msg="Installing features 📚")
        if pbar is not None
        else tasks
    )
    for result in task_iter:
        if not result:
            msg = f"{result=}"
            raise RuntimeError(msg)

    if verbose:
        eti_util.print_colour("\nFinished installing features", "yellow")

    # we parallelise across databases
    writer = eti_genome.fasta_to_hdf5(config=config)
    tasks = eti_util.get_iterable_tasks(
        func=writer,
        series=db_names,
        max_workers=max_workers,
    )
    pbar = progress.child(leave=True) if progress is not None else progress
    task_iter = (
        pbar(tasks, total=len(db_names), msg="Installing 🧬🧬")
        if pbar is not None
        else tasks
    )
    for result in task_iter:
        if not result:
            msg = f"{result=}"
            raise RuntimeError(msg)

    if verbose:
        eti_util.print_colour("\nFinished installing sequences", "yellow")


def local_install_alignments(
    config: eti_config.Config,
    force_overwrite: bool,
    max_workers: int | None,
    verbose: bool = False,
    progress: Progress | None = None,
) -> None:
    # check if alignments are specified in the config
    if not config.align_names:
        if verbose:
            eti_util.print_colour(
                "No alignments specified in the config. Skipping alignment installation.",
                "yellow",
            )
        return

    if force_overwrite:
        shutil.rmtree(config.install_aligns, ignore_errors=True)

    for align_name in config.align_names:
        ingest_aln.install_alignment(
            config=config,
            align_name=align_name,
            progress=progress,
            max_workers=max_workers,
        )

    if verbose:
        eti_util.print_colour("\nFinished installing alignments", "yellow")


def local_install_homology(
    config: eti_config.Config,
    force_overwrite: bool,
    max_workers: int | None,
    verbose: bool = False,
    progress: Progress | None = None,
) -> None:
    # check if homologies are specified in the config
    if not config.homologies:
        if verbose:
            eti_util.print_colour(
                "No homologies specified in the config. Skipping homology installation.",
                "yellow",
            )
        return

    if force_overwrite:
        shutil.rmtree(config.install_homologies, ignore_errors=True)

    config.install_homologies.mkdir(parents=True, exist_ok=True)

    dirnames = []
    for sp in config.species_dbs:
        path = config.staging_homologies / sp
        dirnames.extend(list(path.glob("*.tsv*")))

    max_workers = min(len(dirnames) + 1, max_workers) if max_workers else 1

    if verbose:
        eti_util.print_colour(f"homologies {max_workers=}", "yellow")

    loader = homology_ingest.load_homologies(
        allowed_species=set(config.species_dbs),
    )

    tasks = eti_util.get_iterable_tasks(
        func=loader,
        series=dirnames,
        max_workers=max_workers,
    )
    pbar = progress.child(leave=False) if progress is not None else progress
    task_iter = (
        pbar(tasks, total=len(dirnames), msg="Loading homologies")
        if pbar is not None
        else tasks
    )
    results = {}
    for result in task_iter:
        for rel_type, records in result.items():
            if rel_type not in results:
                results[rel_type] = []
            results[rel_type].extend(records)

    # we merge the homology groups
    items = results.items()
    pbar = progress.child(leave=False) if progress is not None else progress
    agg_iter = pbar(items, msg="Aggregating homologies") if pbar is not None else items
    for rel_type, records in agg_iter:
        results[rel_type] = homology_ingest.merge_grouped(records)

    # write the homology groups to in-memory db
    items = results.items()
    pbar = progress.child(leave=True) if progress is not None else progress
    write_iter = pbar(items, msg="Installing homologies") if pbar is not None else items
    db = homology_ingest.make_homology_aggregator_db()
    for rel_type, records in write_iter:
        db.add_records(records=records, relationship_type=rel_type)

    homology_ingest.write_homology_views(agg=db, outdir=config.install_homologies)
    if verbose:
        eti_util.print_colour("\nFinished installing homologies", "yellow")
