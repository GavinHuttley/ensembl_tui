import pytest

from ensembl_tui import _config as eti_config
from ensembl_tui import _genome as eti_genome
from ensembl_tui import _maf as eti_maf
from ensembl_tui import _util as eti_util


def get_expected_block_ids(alignments):
    expect = set()
    for alignment in alignments:
        concat = "".join(sorted(str(n) for n in alignment))
        expect.add(eti_util.hash64(concat.encode("utf-8")))
    return expect


def test_read(DATA_DIR):
    path = DATA_DIR / "sample.maf"
    blocks = list(eti_maf.parse(path))
    assert len(blocks) == 4
    block_ids, alignments = zip(*blocks, strict=False)
    expect = get_expected_block_ids(alignments)
    assert set(block_ids) == expect


_ONE_BLOCK = """\
##maf version=1
# id: 20060000040557
a
s homo_sapiens.1 100 7 + 1000 AC--GTA---CC
s mouse.2        200 12 + 2000 ACTGGTAGGTCC
s rat.3          300 8 + 3000 -CTGG--AGTC-
"""


@pytest.mark.parametrize("terminator", ("", "\n", "\n\n"))
def test_read_last_record_of_final_block(tmp_path, terminator):
    # real Ensembl maf files end with a blank line, but the final record must
    # survive whether or not the last block is terminated by one
    path = tmp_path / "sample.maf"
    path.write_text(_ONE_BLOCK + terminator)
    (block_id, alignment), *rest = list(eti_maf.parse(path))
    assert not rest
    assert sorted(n.species for n in alignment) == ["homo_sapiens", "mouse", "rat"]


_TWO_BLOCKS = """\
##maf version=1
# id: 20060000040557

a score=1
s homo_sapiens.1 100 7 + 1000 AC--GTA---CC
s mouse.2        200 12 + 2000 ACTGGTAGGTCC

a score=2
s homo_sapiens.1 300 8 + 1000 -CTGG--AGTC-
s rat.3          400 12 + 3000 ACTGGTAGGTCC
i rat.3 N 0 C 0
"""


def test_read_block_boundaries(tmp_path):
    # header lines belong to no block, an "a" line opens one, and everything
    # up to the next "a" line belongs to it. non "s" lines are ignored
    path = tmp_path / "two_blocks.maf"
    path.write_text(_TWO_BLOCKS)
    # unpacking is the assertion that there are exactly two blocks
    first, second = eti_maf.parse(path)
    assert sorted(n.species for n in first[1]) == ["homo_sapiens", "mouse"]
    assert sorted(n.species for n in second[1]) == ["homo_sapiens", "rat"]
    assert first[0] != second[0]


def test_read_decodes_as_it_goes(tmp_path):
    # a byte that cannot be decoded, placed past the reader's chunk size. the
    # blocks before it still come back, which they would not if the file were
    # read and decoded up front
    path = tmp_path / "bad_byte_late.maf"
    padding = "s pad.1 0 10 + 100 " + "A" * 400_000 + "\n"
    good = _TWO_BLOCKS + "".join(f"a score={i}\n{padding}" for i in range(4))
    path.write_bytes(good.encode("utf-8") + b"a score=9\ns bad.1 0 1 + 10 \xff\n")
    gen = eti_maf.parse(path)
    _, alignment = next(gen)
    assert sorted(n.species for n in alignment) == ["homo_sapiens", "mouse"]
    with pytest.raises(UnicodeDecodeError):
        list(gen)


def test_read_yields_blocks_incrementally(tmp_path):
    # blocks come out one at a time, so a malformed block does not stop the
    # ones before it being returned
    path = tmp_path / "trailing_junk.maf"
    path.write_text(_TWO_BLOCKS + "\na score=3\ns not enough fields\n")
    gen = eti_maf.parse(path)
    block_id, alignment = next(gen)
    assert sorted(n.species for n in alignment) == ["homo_sapiens", "mouse"]
    assert block_id
    with pytest.raises(ValueError, match="unpack"):
        list(gen)


@pytest.mark.parametrize(
    "line",
    (
        "s pan_paniscus.11 2 7 + 13 ACTCTCCAGATGA",
        "s pan_paniscus.11 4 7 - 13 ACTCTCCAGATGA",
    ),
)
def test_process_maf_line_plus(line):
    n, s = eti_maf.process_maf_line(line)
    assert s == "ACTCTCCAGATGA"
    # maf is zero based
    assert n.start == 2
    assert n.stop == 2 + 7


def test_process_maf_line_minus():
    from cogent3 import DNA

    seq_plus = "ACCTTTTGGGGGG"
    seq_minus = "CCCCCCAAAAGGT"
    # just check our rc is correct
    assert seq_plus == DNA.rc(seq_minus)
    #            01234567890123
    maf_start = 6
    maf_size = 4
    maf_line = f"s hg38.chr22 {maf_start} {maf_size} - {len(seq_minus)} {seq_minus}"
    n, _ = eti_maf.process_maf_line(maf_line)
    expected_start = seq_plus.find("T")
    expected_stop = expected_start + maf_size
    assert n.start == expected_start
    assert n.stop == expected_stop
    assert seq_plus[n.start : n.stop] == DNA.rc(
        seq_minus[maf_start : maf_start + maf_size],
    )


def get_human_record(blocks, species, strand):
    for _, seqs in blocks:
        for name, seq in seqs.items():
            if name.species == species and name.strand == strand:
                return name, seq
    msg = f"No human record found for {species} with strand {strand}"
    raise ValueError(msg)


def test_compare_maf_with_genome(apes_install_path, apes_maf_install_path):
    config = eti_config.read_installed_cfg(apes_install_path)
    species = "homo_sapiens"
    hsap = eti_genome.load_genome(config=config, species=species)
    chr22 = hsap.seqs["22"]

    # check that the inferred coordinates from an alignment block
    # and the extracted sequence match the genome for those coordinates
    # need to add this file to cached data
    blocks = list(eti_maf.parse(apes_maf_install_path))
    # plus strand
    name, seq = get_human_record(blocks, species, strand="+")
    got = seq.replace("-", "").upper()
    assert got == str(chr22[name.start : name.stop])

    # minus strand
    name, seq = get_human_record(blocks, species, strand="-")
    got = seq.replace("-", "").upper()
    assert got == str(chr22[name.start : name.stop].rc())
