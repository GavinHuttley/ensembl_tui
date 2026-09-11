# parser for MAF, defined at
# https://genome.ucsc.edu/FAQ/FAQformat.html#format5

import typing

from scinexus.io_util import iter_splitlines

from ensembl_tui import _name as eti_name
from ensembl_tui import _util as eti_util


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


def _is_seq_line(line: str) -> bool:
    return line.startswith("s") and "ancestral" not in line[:100]


def _block_id(alignment: dict[eti_name.MafName, str]) -> int:
    # the block ID's are made unique for each alignment by using
    # the str(sorted(str(MafNames)))
    names = "".join(sorted(str(n) for n in alignment))
    return eti_util.hash64(names.encode("utf-8"))


def parse(
    path: eti_util.PathType,
) -> typing.Iterator[tuple[int, dict[eti_name.MafName, str]]]:
    """yields the block id and aligned sequences of each alignment block

    Notes
    -----
    A line beginning with "a" opens a block and the lines up to the next such
    line belong to it. The file is read a chunk at a time and only one block is
    held at once, which is what keeps the large Ensembl files off the heap.

    Two caveats come from iter_splitlines. It compares the size of the file on
    disk against its chunk size, so a compressed file smaller than that is read
    whole however large it expands to. It also opens the file in text mode with
    the encoding sniffed from the first 100 bytes, where this used to decode as
    utf-8 throughout.
    """
    alignment: dict[eti_name.MafName, str] = {}
    in_block = False
    for line in iter_splitlines(path):
        if line.startswith("a"):
            if in_block:
                yield _block_id(alignment), alignment
            alignment = {}
            in_block = True
        elif in_block and _is_seq_line(line):
            n, seq = process_maf_line(line)
            alignment[n] = seq

    if in_block:
        yield _block_id(alignment), alignment
