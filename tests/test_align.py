import pickle

import duckdb
import numpy
import pytest

from ensembl_tui import _align as eti_align
from ensembl_tui import _annotation as eti_annots
from ensembl_tui import _config as eti_config
from ensembl_tui import _genome as eti_genome
from ensembl_tui import _ingest_align as eti_ingest_align


def make_gene_attr(records: list[dict]) -> eti_annots.GeneView:
    schema = (
        "stable_id TEXT",
        "biotype TEXT",
        "seqid TEXT",
        "start INTEGER",
        "stop INTEGER",
        "strand TINYINT",
        "canonical_transcript_id INTEGER",
        "symbol TEXT",
        "gene_id INTEGER",
        "description TEXT",
    )
    columns = [c.split()[0] for c in schema]
    sql = f"""CREATE TABLE IF NOT EXISTS gene_attr ({",".join(schema)})"""
    conn = duckdb.connect(":memory:")
    conn.sql(sql)
    sql = "INSERT INTO gene_attr VALUES (?, ?, ?, ?, ?, ?, ?, ?, ?, ?)"
    rows = [[r.get(c) for c in columns] for r in records]
    conn.executemany(sql, parameters=rows)
    return eti_annots.GeneView(source=":memory:", db=conn)


def get_annotation_db() -> dict[str, eti_annots.Annotations]:
    gene_attr = {
        "s1": make_gene_attr(
            [
                {
                    "seqid": "s1",
                    "biotype": "protein_coding",
                    "stable_id": "not-on-s2",
                    "start": 4,
                    "stop": 7,
                },
            ],
        ),
        "s2": make_gene_attr(
            [
                {
                    "seqid": "s2",
                    "biotype": "protein_coding",
                    "stable_id": "includes-s2-gap",
                    "start": 2,
                    "stop": 6,
                },
            ],
        ),
        "s3": make_gene_attr(
            [
                {
                    "seqid": "s3",
                    "biotype": "protein_coding",
                    "stable_id": "includes-s3-gap",
                    "start": 22,
                    "stop": 27,
                },
            ],
        ),
    }
    for sp, gv in gene_attr.items():
        gene_attr[sp] = eti_annots.Annotations(source=":memory:", genes=gv)
    return gene_attr


def small_seqs():
    from cogent3 import make_aligned_seqs

    seqs = {
        "s1": "GTTGAAGTAGTAGAAGTTCCAAATAATGAA",
        "s2": "GTG------GTAGAAGTTCCAAATAATGAA",
        "s3": "GCTGAAGTAGTGGAAGTTGCAAAT---GAA",
    }
    return make_aligned_seqs(
        data=seqs,
        moltype="dna",
        array_align=False,
        info={"species": {"s1": "human", "s2": "mouse", "s3": "dog"}},
        new_type=True,
    )


def make_records(start, end, block_id):
    aln = small_seqs()[start:end]
    records = []
    species = aln.info.species
    for seq in aln.seqs:
        seqid, seq_start, seq_end, seq_strand = seq.parent_coordinates(seq_coords=True)
        gs = seq.gapped_seq
        imap, s = gs.parse_out_gaps()
        if imap.num_gaps:
            gap_spans = numpy.array(
                [imap.gap_pos, imap.get_gap_lengths()],
                dtype=numpy.int32,
            ).T
        else:
            gap_spans = numpy.array([], dtype=numpy.int32)
        record = eti_align.AlignRecord(
            source="blah",
            species=species[seq.name],
            block_id=block_id,
            seqid=seqid,
            start=seq_start,
            stop=seq_end,
            strand=seq_strand,
            gap_spans=gap_spans,
        )
        records.append(record)
    return records


@pytest.fixture
def small_records():
    return make_records(1, 5, 0)


def empty_align_agg_gap_store():
    return eti_ingest_align.make_alignment_aggregator_db()


def test_aligndb_records_match_input(small_records):
    import copy

    orig_records = copy.deepcopy(small_records)
    agg = empty_align_agg_gap_store()
    eti_ingest_align.add_records(records=small_records, conn=agg)
    db = eti_align.AlignDb(db=agg, source=":memory:")
    got = next(iter(db.get_records_matching(species="human", seqid="s1")))
    assert got == set(orig_records)


def test_aligndb_records_skip_duplicated_block_ids(small_records):
    agg = empty_align_agg_gap_store()
    eti_ingest_align.add_records(conn=agg, records=small_records)
    sql = "SELECT COUNT(*) FROM align_blocks"
    count = agg.sql(sql).fetchone()[0]
    eti_ingest_align.add_records(conn=agg, records=small_records)
    assert agg.sql(sql).fetchone()[0] == count


# fixture to make synthetic GenomeSeqsDb and alignment db
# based on a given alignment
@pytest.fixture
def genomedbs_aligndb(small_records):
    agg = empty_align_agg_gap_store()
    eti_ingest_align.add_records(records=small_records, conn=agg)
    align_db = eti_align.AlignDb(source=":memory:", db=agg)
    seqs = small_seqs().degap()
    species = seqs.info.species
    data = seqs.to_dict()
    genomes = {}
    for name, seq in data.items():
        genome = eti_genome.SeqsDataHdf5(
            source=f"{name}",
            species=species[name],
            mode="w",
            in_memory=True,
        )
        genome.add_records(records=[(name, seq)])
        genomes[species[name]] = eti_genome.Genome(
            seqs=genome,
            annots=None,
            species=species[name],
        )

    return genomes, align_db


def test_building_alignment(genomedbs_aligndb, namer):
    genomes, align_db = genomedbs_aligndb
    got = next(
        iter(
            eti_align.get_alignment(
                align_db,
                genomes,
                ref_species="mouse",
                seqid="s2",
                namer=namer,
            ),
        ),
    )
    orig = small_seqs()[1:5]
    assert got.to_dict() == orig.to_dict()


@pytest.mark.parametrize(
    "kwargs",
    ({"ref_species": "dodo", "seqid": "s2"},),
)
def test_building_alignment_invalid_details(genomedbs_aligndb, kwargs):
    genomes, align_db = genomedbs_aligndb
    with pytest.raises(ValueError):
        list(eti_align.get_alignment(align_db, genomes, **kwargs))


def make_sample(two_aligns=False):
    aln = small_seqs()
    species = aln.info.species
    # make annotation db's
    annot_dbs = get_annotation_db()

    # we will reverse complement the s2 genome compared to the original
    # this means our coordinates for alignment records from that genome
    # also need to be rc'ed
    genomes = {}
    for seq in aln.seqs:
        name = seq.name
        seq = seq.seq
        if seq.name == "s2":
            seq = seq.rc()
            s2_genome = str(seq)
        genome = eti_genome.SeqsDataHdf5(
            source=f"{name}",
            mode="w",
            in_memory=True,
            species=species[seq.name],
        )
        genome.add_records(records=[(name, str(seq))])
        genomes[species[name]] = eti_genome.Genome(
            seqs=genome,
            annots=annot_dbs[name],
            species=species[name],
        )

    # define two alignment blocks that incorporate features
    align_records = _update_records(s2_genome, aln, "0", 1, 12)
    if two_aligns:
        align_records += _update_records(s2_genome, aln, "1", 22, 30)

    maker = eti_ingest_align.make_alignment_aggregator_db()
    eti_ingest_align.add_records(conn=maker, records=align_records)
    align_db = eti_align.AlignDb(db=maker, source=":memory:")
    return genomes, align_db


def _update_records(s2_genome, aln, block_id, start, end):
    # start, stop are the coordinates used to slice the alignment
    align_records = make_records(start, end, block_id)
    # in the alignment, s2 is in reverse complement relative to its genome
    # In order to be sure what "genome" coordinates are for s2, we first slice
    # the alignment
    aln = aln[start:end]
    # then get the ungapped sequence
    s2 = aln.get_seq("s2")
    # and reverse complement it ...
    selected = s2.rc()
    # so we can get the genome coordinates for this segment on the s2 genome
    start = s2_genome.find(str(selected))
    end = start + len(selected)
    for record in align_records:
        if record.seqid == "s2":
            record.start = start
            record.stop = end
            record.strand = -1
            break
    return align_records


@pytest.mark.parametrize(
    "start_end",
    [
        (None, None),
        (None, 11),
        (3, None),
        (3, 13),
    ],
)
@pytest.mark.parametrize(
    "species_coord",
    [
        ("human", "s1"),
        ("dog", "s3"),
    ],
)
def test_select_alignment_plus_strand(species_coord, start_end, namer):
    species, seqid = species_coord
    start, end = start_end
    aln = small_seqs()
    # the sample alignment db has an alignment block that
    # starts at 1 and ends at 12. The following slice is to
    # get the expected answer
    expect = aln[max(1, start or 1) : min(end or 12, 12)]
    # one sequence is stored in reverse complement
    genomes, align_db = make_sample()
    got = list(
        eti_align.get_alignment(
            align_db=align_db,
            genomes=genomes,
            ref_species=species,
            seqid=seqid,
            ref_start=start,
            ref_end=end,
            namer=namer,
        ),
    )
    assert len(got) == 1
    assert got[0].to_dict() == expect.to_dict()


@pytest.mark.parametrize(
    "start_end",
    [
        (None, None),
        (None, 5),
        (2, None),
        (2, 7),
    ],
)
def test_select_alignment_minus_strand(start_end, namer):
    species, seqid = "mouse", "s2"
    start, end = start_end
    aln = small_seqs()
    ft = aln.add_feature(
        biotype="custom",
        name="selected",
        seqid="s2",
        spans=[(max(1, start or 0), min(end or 12, 12))],
    )
    expect = aln[ft.map.start : min(ft.map.end, 12)]
    # mouse sequence is on minus strand, so need to adjust
    # coordinates for query
    s2 = aln.get_seq("s2")
    s2_ft = next(iter(s2.get_features(name="selected")))
    if not any([start is None, end is None]):
        start = len(s2) - s2_ft.map.end
        end = len(s2) - s2_ft.map.start
    elif start == None != end:  # noqa E711
        start = len(s2) - s2_ft.map.end
        end = None
    elif start != None == end:  # noqa E711
        end = len(s2) - s2_ft.map.start
        start = None

    # mouse sequence is on minus strand, so need to adjust
    # coordinates for query

    genomes, align_db = make_sample(two_aligns=False)
    got = list(
        eti_align.get_alignment(
            align_db=align_db,
            genomes=genomes,
            ref_species=species,
            seqid=seqid,
            ref_start=start,
            ref_end=end,
            namer=namer,
        ),
    )
    # drop the strand info
    assert len(got) == 1, f"{s2_ft=}"
    assert got[0].to_dict() == expect.to_dict()


@pytest.mark.parametrize(
    "coord",
    [
        ("human", "s1", None, 11),  # finish within
        ("human", "s1", 3, None),  # start within
        ("human", "s1", 3, 9),  # within
        ("human", "s1", 3, 13),  # extends past
    ],
)
def test_get_alignment_features(coord):
    kwargs = dict(
        zip(("ref_species", "seqid", "ref_start", "ref_end"), coord, strict=False),
    )
    genomes, align_db = make_sample(two_aligns=False)
    got = next(
        iter(eti_align.get_alignment(align_db=align_db, genomes=genomes, **kwargs)),
    )
    assert len(got.annotation_db) == 1


@pytest.mark.parametrize(
    "coord",
    [
        ("human", "s1", None, 11),  # finish within
        ("human", "s1", 3, None),  # start within
        ("human", "s1", 3, 9),  # within
        ("human", "s1", 3, 13),  # extends past
    ],
)
def test_get_alignment_masked_features(coord):
    kwargs = dict(
        zip(("ref_species", "seqid", "ref_start", "ref_end"), coord, strict=False),
    )
    kwargs["mask_features"] = ["gene"]
    genomes, align_db = make_sample(two_aligns=False)
    got = next(
        iter(eti_align.get_alignment(align_db=align_db, genomes=genomes, **kwargs)),
    )
    assert len(got.annotation_db) == 1


@pytest.mark.parametrize(
    "coord",
    [
        ("human", "s1", None, 11),  # finish within
        ("human", "s1", 3, None),  # start within
        ("human", "s1", 3, 9),  # within
        ("human", "s1", 3, 13),  # extends past
    ],
)
def test_align_db_get_records(coord):
    kwargs = dict(zip(("species", "seqid", "start", "stop"), coord, strict=False))
    # records are, we should get a single hit from each query
    # [('blah', 0, 'human', 's1', 1, 12, '+', array([], dtype=int32)),
    _, align_db = make_sample(two_aligns=True)
    got = list(align_db.get_records_matching(**kwargs))
    assert len(got) == 1


@pytest.mark.parametrize(
    "coord",
    [
        ("human", "s1"),
        ("mouse", "s2"),
        ("dog", "s3"),
    ],
)
def test_align_db_get_records_required_only(coord):
    kwargs = dict(zip(("species", "seqid"), coord, strict=False))
    # two hits for each species
    _, align_db = make_sample(two_aligns=True)
    got = list(align_db.get_records_matching(**kwargs))
    assert len(got) == 2


@pytest.mark.parametrize(
    "coord",
    [
        ("human", "s2"),
        ("mouse", "xx"),
        ("blah", "s3"),
    ],
)
def test_align_db_get_records_no_matches(coord):
    kwargs = dict(zip(("species", "seqid"), coord, strict=False))
    # no hits at all
    _, align_db = make_sample()
    got = list(align_db.get_records_matching(**kwargs))
    assert not len(got)


def test_get_species():
    _, align_db = make_sample()
    assert set(align_db.get_species_names()) == {"dog", "human", "mouse"}


def test_write_alignments():
    genomes, align_db = make_sample(two_aligns=True)
    locations = eti_genome.get_gene_segments(
        annot_db=genomes["human"].annotation_db,
        species="human",
        stableids=["not-on-s2"],
    )
    app = eti_align.construct_alignment(
        align_db=align_db,
        genomes=genomes,
    )
    aln = app(locations[0])  # pylint: disable=not-callable
    assert len(aln[0]) == 3


@pytest.fixture
def db_align(DATA_DIR, tmp_dir):
    staging_path = tmp_dir / "staging"
    staging_path.mkdir(parents=True, exist_ok=True)
    align_name = "align_name"
    install_path = tmp_dir / "install"

    cfg = eti_config.Config(
        host="localhost",
        remote_path="",
        staging_path=staging_path,
        install_path=install_path,
        species_dbs={},
        release="113",
        align_names=[align_name],
        tree_names=[],
        homologies=True,
    )
    align_dir = cfg.staging_aligns / align_name
    align_dir.mkdir(parents=True, exist_ok=True)
    maf = align_dir / f"{align_name}.maf"
    maf.write_text((DATA_DIR / "tiny.maf").read_text())

    parquet_path = eti_ingest_align.install_alignment(cfg, align_name)
    return eti_align.AlignDb(source=parquet_path.parent)


def test_db_align(db_align):
    orig = len(db_align)
    source = db_align.source
    db_align.close()
    got = eti_align.AlignDb(source=source)
    assert len(got) == orig


@pytest.mark.parametrize("func", [str, repr])
def test_db_align_repr(db_align, func):
    got = func(db_align)
    assert "AlignDb" in got


def test_pickling_db(db_align):
    # should not fail
    pkl = pickle.dumps(db_align)  # nosec B301
    upkl = pickle.loads(pkl)  # nosec B301  # noqa: S301
    assert db_align.source == upkl.source


def test_aligndb_post_init_failure(tmp_path):
    with pytest.raises(FileNotFoundError):
        eti_align.AlignDb(source="/non/existent/directory")

    outfile = tmp_path / "textfile.txt"
    outfile.write_text("blah")
    with pytest.raises(OSError):
        eti_align.AlignDb(source="/non/existent/directory")


def test_aligndb_close(db_align):
    db_align.close()
    with pytest.raises(duckdb.duckdb.ConnectionException):
        db_align.num_records()


def test_load_align_records():
    maf_record = {
        "species": "chlorocebus_sabaeus",
        "seqid": "28",
        "start": 2617243,
        "stop": 2645355,
        "strand": "-",
        "block_id": 20060000081647,
        "source": "10_primates.epo.other_6.maf.gz",
        # mixed case seq because this must be upper cased in seq2gap for cogent3
        # to deem it a valid seq
        "seq": "-acCAGGAAG",
    }
    got = eti_ingest_align.seq2gaps(maf_record)
    assert (got.gap_spans == numpy.array([[0, 1]], dtype=numpy.int32)).all()
