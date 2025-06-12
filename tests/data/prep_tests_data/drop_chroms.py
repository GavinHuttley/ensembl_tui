import pathlib

import click
import cogent3

from ensembl_tui import _config as eti_config
from ensembl_tui import _genome as eti_genome


@click.command()
@click.argument("install_dir", type=pathlib.Path)
@click.option("--check", is_flag=True)
def main(install_dir, check):
    seqid = "22"
    cfg = eti_config.read_installed_cfg(install_dir)

    for db in cfg.list_genomes():
        genome_dir = cfg.installed_genome(db)

        src = genome_dir / eti_genome.SEQ_STORE_NAME
        seqs = cogent3.load_unaligned_seqs(src, moltype="dna", new_type=True)
        if check:
            print(f"{db=}  {seqs.names}")
            continue
        dest = src.with_suffix(f"{src.suffix}-temp")
        seqs = seqs.take_seqs([seqid])
        seqs.write(dest)
        dest.rename(str(dest).replace("-temp", ""))


if __name__ == "__main__":
    main()
