from ensembl_tui import _genome as eti_genome
from ensembl_tui import _species as eti_species


def test_get_gene_segments_limit(yeast_db):
    segs = eti_genome.get_gene_segments(
        annot_db=yeast_db,
        species="saccharomyces_cerevisiae",
        limit=10,
    )
    assert len(segs) == 10


def test_get_seq_feature_seq_correct(yeast):
    seq = yeast.seqs["III"][309069:310155]
    raw_seq = str(seq)
    assert raw_seq.startswith("ATGCTTTACCCAGAAAAATTTCA")  # expected from ensembl.org
    assert raw_seq.endswith("ATAAGAAATTCCATAAATAG")


def test_get_seq_feature_seq_correct_name(yeast):
    # need to modify cogent3 so it applies the feature name
    # to the new sequence
    seq = yeast.seqs["III"][309069:310155]
    got = next(iter(seq.get_features()))
    feat_seq = got.get_slice()
    assert feat_seq.name == "YCR105W"


def test_get_gene_table_for_species(yeast_db):
    from cogent3.core.table import Table

    # we do not check values here, only the Type and that we have > 0 records
    got = eti_genome.get_gene_table_for_species(annot_db=yeast_db)
    assert isinstance(got, Table)
    assert len(got) > 0


def test_get_species_gene_summary(yeast_db):
    from cogent3.core.table import Table

    got = eti_genome.get_species_gene_summary(
        annot_db=yeast_db, species_map=eti_species.make_species_map(None)
    )
    # we do not check values here, only the Type and that we have > 0 records
    assert isinstance(got, Table)
    assert len(got) > 0
    assert "biotype" in got.header


def test_get_species_repeat_summary(yeast_db):
    from cogent3.core.table import Table

    got = eti_genome.get_species_repeat_summary(
        annot_db=yeast_db, species_map=eti_species.make_species_map(None)
    )
    # we do not check values here, only the Type and that we have > 0 records
    assert isinstance(got, Table)
    assert len(got) > 0
    assert "repeat_type" in got.header


def test_genome_coord_names(yeast_db):
    counts = yeast_db.count_distinct(seqid=True)
    assert counts.shape[0] == 17


def test_get_seq(yeast):
    mt = yeast.seqs["Mito"]
    assert len(mt) == 85779  # number for ensembl.org


def test_species_setting(yeast):
    # note that species are converted into the Ensembl db prefix
    assert yeast.info.species == "saccharomyces_cerevisiae"


def test_gene_description(worm_db):
    gene = next(iter(worm_db.genes.get_by_stable_id(stable_id="WBGene00010790")))
    assert "alcohol dehydrogenase" in gene.description.lower()


def test_featuredb_num_records(worm_db):
    db = worm_db
    num_repeats = db.repeats.num_records()
    num_genes = db.genes.num_records()
    assert db.num_records() == num_repeats + num_genes
    assert num_repeats > 1000
    assert num_genes > 1000
    num_biotypes = db.biotypes.num_records()
    assert num_biotypes > 10


def test_get_features_matching(yeast_db):
    got = list(yeast_db.get_features_matching(biotype="rRNA"))
    assert len(got) > 20  # number from ensembl.org


def test_get_feature_child_parent(worm_db):
    db = worm_db.genes
    gene = next(iter(db.get_by_stable_id(stable_id="WBGene00010790")))
    transcript = next(iter(db.get_feature_children(gene)))
    cds = next(iter(db.get_feature_children(transcript)))
    assert gene.canonical_transcript_id == transcript.transcript_id == cds.transcript_id
    got = db.get_feature_parent(cds)
    assert isinstance(got, type(transcript))
    assert got.stable_id == transcript.stable_id
    got = db.get_feature_parent(transcript)
    assert isinstance(got, type(gene))
    assert got.stable_id == gene.stable_id


def test_genome_segment():
    segment = eti_genome.genome_segment(
        species="abcd_efg",
        seqid="1",
        start=20,
        stop=40,
        strand="+",
    )
    assert segment.unique_id == "abcd_efg-1-20-40"
    segment = eti_genome.genome_segment(
        species="abcd_efg",
        seqid="1",
        start=20,
        stop=40,
        strand="+",
        unique_id="gene:NICE",
    )
    assert segment.unique_id == "NICE"


def test_get_gene_segments_value(worm_db):
    segments = eti_genome.get_gene_segments(
        annot_db=worm_db,
        species="caenorhabditis_elegans",
    )
    assert len(segments) == len(
        list(worm_db.get_features_matching(biotype="protein_coding")),
    )
    assert "WBGene00010790" in {s.unique_id for s in segments}


def test_get_gene_segments_stableids(worm_db):
    segments = eti_genome.get_gene_segments(
        annot_db=worm_db,
        species="caenorhabditis_elegans",
        stableids=["WBGene00010790"],
    )
    assert len(segments) == 1
    segment = segments[0]
    assert segment.unique_id == "WBGene00010790"
    assert segment.source == "WBGene00010790"


def test_get_features(yeast):
    features = list(yeast.get_features(biotype="rRNA", limit=10))
    # because of how feature querying works, one seqid at a time,
    # the limit argument will be applied to each seqid
    assert len(features) >= 10


def test_get_celegans_cds(worm):
    cds = next(iter(worm.get_features(name="WBGene00021347", biotype="cds")))
    seq = cds.get_slice()
    aa = seq.get_translation()
    assert aa == "MIIPIRCFTCGKVIGDKWETYLGFLQSEYSEGDALDALGLRRYCCRRMLLAHVDLIEKLLNYHPLEK"


def test_yeast_repeat(yeast):
    repeat = next(
        iter(
            yeast.get_features(
                biotype="dust",
                name="dust",
                seqid="VI",
                limit=20,
                start=8680,
                stop=9000,
            ),
        ),
    )
    seq = repeat.get_slice()
    assert str(seq) == "AAAAAAAAAA"
