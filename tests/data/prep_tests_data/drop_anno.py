import pathlib

import click

from ensembl_tui import _annotation as eti_ann
from ensembl_tui import _config as eti_config


@click.command()
@click.argument("install_dir", type=pathlib.Path)
def main(install_dir):
    seqid = "22"
    cfg = eti_config.read_installed_cfg(install_dir)
    for db in cfg.list_genomes():
        # drop non-22 chromosome repeat features and consensus
        anno_db = eti_ann.Annotations(source=cfg.installed_genome(db))
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

        old_repeat_feature_path = rv.source / "repeat_feature.parquet"
        new_repeat_feature_path = old_repeat_feature_path.with_suffix(".new-parquet")
        rv.conn.sql(
            f"COPY repeat_feature TO '{new_repeat_feature_path}' (FORMAT 'parquet')",
        )

        # finally, remove any repeat_consensus rows that are non-22 chromosome
        rv.conn.sql(
            "DELETE FROM repeat_consensus WHERE repeat_consensus_id NOT IN (SELECT DISTINCT repeat_consensus_id FROM repeat_feature)",
        )
        old_repeat_con_path = rv.source / "repeat_consensus.parquet"
        new_repeat_con_path = old_repeat_con_path.with_suffix(".new-parquet")
        rv.conn.sql(
            f"COPY repeat_consensus TO '{new_repeat_con_path}' (FORMAT 'parquet')",
        )


if __name__ == "__main__":
    main()
