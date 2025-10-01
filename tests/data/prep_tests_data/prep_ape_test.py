import pathlib
import shutil

import click
import cogent3

from ensembl_tui import _align as eti_align
from ensembl_tui import _annotation as eti_ann
from ensembl_tui import _config as eti_config
from ensembl_tui import _genome as eti_genome


def drop_aligns(align_path: pathlib.Path, seqid: str = "22"):
    align_db = eti_align.AlignDb(source=align_path)
    align_db.conn.sql(f"DELETE FROM align_blocks WHERE seqid != {seqid!r}")
    dest = next(iter(align_path.glob(f"*{eti_align.ALIGN_STORE_SUFFIX}")))
    align_db.conn.sql(f"COPY align_blocks TO '{dest}' (FORMAT 'parquet')")


def drop_annotations(genome_dir: pathlib.Path, seqid: str = "22"):
    anno_db = eti_ann.Annotations(source=genome_dir)
    rv = anno_db.repeats

    # first figure out the seq_region_id corresponding to chromosome 22
    query = f"""
        SELECT sr.seq_region_id
        FROM seq_region sr
        JOIN coord_system cs ON sr.coord_system_id = cs.coord_system_id
        WHERE sr.name = '{seqid}' AND cs.attrib = 'default_version'
        """
    (chrom_id,) = rv.conn.sql(query).fetchone()
    # then "delete" those rows
    rv.conn.sql(f"DELETE FROM repeat_feature WHERE seq_region_id != {chrom_id}")

    repeat_feature_path = rv.source / "repeat_feature.parquet"
    # now overwrite the repeat_feature file
    rv.conn.sql(
        f"COPY repeat_feature TO '{repeat_feature_path}' (FORMAT 'parquet')",
    )

    # finally, remove any repeat_consensus rows that are non-22 chromosome
    rv.conn.sql(
        "DELETE FROM repeat_consensus WHERE repeat_consensus_id NOT IN (SELECT DISTINCT repeat_consensus_id FROM repeat_feature)",
    )
    repeat_con_path = rv.source / "repeat_consensus.parquet"
    # now overwrite the repeat_consensus file
    rv.conn.sql(
        f"COPY repeat_consensus TO '{repeat_con_path}' (FORMAT 'parquet')",
    )

    # remove any transcripts that are not on seqid
    gv = anno_db.genes
    transcript_attr_path = gv.source / "transcript_attr.parquet"
    sql = f"DELETE FROM transcript_attr WHERE seqid != {seqid!r}"
    gv.conn.sql(sql)
    gv.conn.sql(
        f"COPY transcript_attr TO '{transcript_attr_path}' (FORMAT 'parquet')",
    )

    # remove any genes that are not on seqid
    gene_attr_path = gv.source / "gene_attr.parquet"
    gv.conn.sql(f"DELETE FROM gene_attr WHERE seqid != {seqid!r}")
    gv.conn.sql(
        f"COPY gene_attr TO '{gene_attr_path}' (FORMAT 'parquet')",
    )


def drop_chrom(genome_dir: pathlib.Path, seqid: str = "22", check: bool = True):
    src = genome_dir / eti_genome.SEQ_STORE_NAME
    seqs = cogent3.load_unaligned_seqs(src, moltype="dna")
    if check:
        print(f"{genome_dir.name=}  {seqs.names}")
        return

    seqs = seqs.take_seqs([seqid])
    new_path = genome_dir / f"{seqid}.c3h5u"
    if new_path.exists():
        print(f"'{new_path}' exists, removing it")
        new_path.unlink()
    seqs.write(new_path)
    seqs.storage.close()
    del seqs
    new_path.rename(src)


def copy_maf(cfg_path: pathlib.Path, dest_dir: pathlib.Path, seqid: str = "22"):
    cfg = eti_config.read_config(config_path=cfg_path)
    align_path = cfg.staging_aligns / cfg.align_names[0]
    maf_file = min(align_path.glob(f"*.{seqid}*.maf.*"), key=lambda p: p.stat().st_size)
    dest_dir.mkdir(parents=True, exist_ok=True)
    dest_path = dest_dir / maf_file.name
    shutil.copy(maf_file, dest_path)


_click_command_opts = {
    "no_args_is_help": True,
    "context_settings": {"show_default": True},
}


@click.command(**_click_command_opts)
@click.argument("install_dir", type=pathlib.Path)
@click.argument("download_dir", type=pathlib.Path)
@click.option("--check", is_flag=True)
@click.option("--release", required=True, type=str)
def main(install_dir, download_dir, check, release):
    seqid = "22"
    # copy smallest maf file for chrom 22
    copy_maf(
        download_dir / "downloaded.cfg",
        pathlib.Path(f"apes-{release}-maf"),
        seqid=seqid,
    )
    cfg = eti_config.read_installed_cfg(install_dir)

    for db in cfg.list_genomes():
        genome_dir = cfg.installed_genome(db)
        drop_chrom(genome_dir, seqid=seqid, check=check)
        if check:
            continue
        drop_annotations(genome_dir=genome_dir, seqid=seqid)

    # now trim the alignment to just seqid
    align_path = cfg.path_to_alignment("*primate*", eti_align.ALIGN_STORE_SUFFIX)
    drop_aligns(align_path, seqid=seqid)


if __name__ == "__main__":
    main()
