import pytest

from ensembl_tui import _config as eti_config
from ensembl_tui import _genome as eti_genome
from ensembl_tui import _maf as eti_maf


def test_read(DATA_DIR):
    path = DATA_DIR / "sample.maf"
    blocks = list(eti_maf.parse(path))
    assert len(blocks) == 4
    block_ids = {b for b, *_ in blocks}
    assert block_ids == {20060000040557, 20060000042317, 20060000132559, 20060000102888}


def test_process_id_line():
    got = eti_maf.process_id_line("# id: 20060000042317 \n")
    assert got == 20060000042317


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
    n, s = eti_maf.process_maf_line(maf_line)
    expected_start = seq_plus.find("T")
    expected_stop = expected_start + maf_size
    assert n.start == expected_start
    assert n.stop == expected_stop
    assert seq_plus[n.start : n.stop] == DNA.rc(
        seq_minus[maf_start : maf_start + maf_size],
    )


