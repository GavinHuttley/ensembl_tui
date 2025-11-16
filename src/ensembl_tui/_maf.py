# parser for MAF, defined at
# https://genome.ucsc.edu/FAQ/FAQformat.html#format5

import typing

from cogent3 import open_

from ensembl_tui import _name as eti_name
from ensembl_tui import _util as eti_util


def _get_alignment_block_indices(data: list[str]) -> list[tuple[int, int]]:
    blocks = []
    start = None
    for i, line in enumerate(data):
        if line.startswith("a"):
            if start is not None:
                blocks.append((start, i))
            start = i

    if start is None:
        return []

    blocks.append((start, i))
    return blocks


def process_maf_line(line: str) -> tuple[eti_name.MafName, str]:
    # after the s token we have src.seqid, start, size, strand, src_size, seq
    _, src_coord, start, size, strand, coord_length, seq = line.strip().split()
    species, coord = src_coord.split(".", maxsplit=1)
    start, size, coord_length = int(start), int(size), int(coord_length)
    if strand == "-":
        start = coord_length - start - size

    stop = start + size
    n = eti_name.MafName(
        species=species,
        seqid=coord,
        start=start,
        stop=stop,
        strand=strand,
        coord_length=coord_length,
    )
    return n, seq


def _get_seqs(lines: list[str]) -> dict[eti_name.MafName, str]:
    alignment = {}
    for line in lines:
        if not line.startswith("s") or "ancestral" in line[:100]:
            continue
        n, seq = process_maf_line(line)
        alignment[n] = seq
    return alignment


def parse(
    path: eti_util.PathType,
) -> typing.Iterator[tuple[int, dict[eti_name.MafName, str]]]:
    with open_(path, mode="rb") as infile:
        data = infile.read()

    # the block ID's are made unique for each alignment by using
    # the str(sorted(str(MafNames)))
    data = data.decode("utf-8").splitlines()
    for block_start, block_end in _get_alignment_block_indices(data):
        alignment = _get_seqs(data[block_start + 1 : block_end])
        names = "".join(sorted(str(n) for n in alignment))
        block_id = eti_util.hash64(names.encode("utf-8"))
        yield block_id, alignment
