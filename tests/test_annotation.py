import pytest

import ensembl_tui._annotation as eti_anno


@pytest.fixture
def gene_view(genome_dir):
    return eti_anno.GeneView(source=genome_dir)


def test_select_protein_coding(worm_genes):
    prot = list(
        worm_genes.get_features_matching(biotype="protein_coding", seqid="I"),
    )
    assert len(prot) > 2000, len(prot)
    prot = list(
        worm_genes.get_features_matching(biotype="protein_coding", seqid="MtDNA"),
    )
    assert len(prot) == 12


def test_gene_view_count_distinct(worm_genes):
    seqid_counts = worm_genes.count_distinct(seqid=True)
    assert seqid_counts.shape[0] == 7  # c. elegans has 7 chromosomes
    assert seqid_counts.shape[1] == 2

    biotype_counts = worm_genes.count_distinct(biotype=True)
    assert biotype_counts.shape[0] >= 10  # from direct inspection of db
    assert biotype_counts.shape[1] == 2

    biotype_by_seqid = worm_genes.count_distinct(seqid=True, biotype=True)
    assert biotype_by_seqid.shape[0] > seqid_counts.shape[0]
    assert biotype_by_seqid.shape[1] == 3

    # constrained by string value
    seqid_biotype = worm_genes.count_distinct(biotype=True, seqid="I")
    assert seqid_biotype.shape[0] > 1
    assert seqid_biotype.shape[1] == 3


def test_biotype_view(worm_biotypes):
    bt = worm_biotypes
    distinct = bt.distinct
    assert "protein_coding" in distinct
    assert "miRNA" in distinct
    counts = bt.count_distinct()
    assert counts["protein_coding", "count"] > 10_000


def test_repeat_view_count_distinct(worm_repeats):
    seqid_counts = worm_repeats.count_distinct(seqid=True)
    assert seqid_counts.shape[0] == 7
    assert seqid_counts.shape[1] == 2

    class_counts = worm_repeats.count_distinct(repeat_class=True)
    assert class_counts.shape[0] > 10
    assert class_counts.shape[1] == 2

    type_counts = worm_repeats.count_distinct(repeat_type=True)
    assert type_counts.shape[0] == 10
    assert type_counts.shape[1] == 2

    class_by_seqid = worm_repeats.count_distinct(seqid=True, repeat_class=True)
    assert class_by_seqid.shape[0] > seqid_counts.shape[0]
    assert class_by_seqid.shape[1] == 3

    type_by_seqid = worm_repeats.count_distinct(seqid=True, repeat_type=True)
    assert type_by_seqid.shape[0] > seqid_counts.shape[0]
    assert type_by_seqid.shape[1] == 3

    type_by_class = worm_repeats.count_distinct(repeat_class=True, repeat_type=True)
    assert type_by_class.shape[0] > 10
    assert type_by_class.shape[1] == 3

    seqid_by_class_by_type = worm_repeats.count_distinct(
        seqid=True,
        repeat_class=True,
        repeat_type=True,
    )
    assert seqid_by_class_by_type.shape[0] > 10
    assert seqid_by_class_by_type.shape[1] == 4

    # constrained by string value
    seqid_class = worm_repeats.count_distinct(repeat_class=True, seqid="I")
    assert seqid_class.shape[0] > 1
    assert seqid_class.shape[1] == 3


def test_create_genome(genome_dir):
    g = eti_anno.Annotations(source=genome_dir)
    prot = list(g.get_features_matching(biotype="protein_coding"))
    assert prot


def test_get_feature_by_symbol(genome_dir):
    g = eti_anno.GeneView(source=genome_dir)
    features = list(g.get_by_symbol(symbol="sms-2"))
    assert features
    # validate that all genes have a span with single start and stop
    assert all(ft["spans"].shape[1] == 2 for ft in features)
    assert all(ft["symbol"] == "sms-2" for ft in features)
    assert all(ft["stable_id"] == "WBGene00004893" for ft in features)
    assert all(ft["name"] == "WBGene00004893" for ft in features)


def test_get_feature_by_description(genome_dir):
    g = eti_anno.GeneView(source=genome_dir)
    features = list(g.get_by_description(description="Alcohol dehydrogenase"))
    assert features
    assert all("alcohol dehydrogenase" in ft["description"].lower() for ft in features)


def test_canonical_cds(genome_dir):
    g = eti_anno.GeneView(source=genome_dir)
    cds = next(
        iter(
            g.get_features_matching(
                stable_id="WBGene00004893",
                biotype="cds",
                canonical=True,
            ),
        ),
    )
    assert cds.stable_id == "F53H8.4.1"  # ID from ensembl.org
    assert cds.spans.shape[1] == 2
    assert cds.spans.shape[0] > 1


@pytest.mark.parametrize("biotype", ["mrna", "transcript"])
def test_canonical_transcript(genome_dir, biotype):
    g = eti_anno.GeneView(source=genome_dir)
    feature = next(
        iter(g.get_features_matching(stable_id="WBGene00004893", biotype=biotype)),
    )
    assert feature.stable_id == "F53H8.4.1"  # ID from ensembl.org
    assert feature.spans.shape[1] == 2
    assert feature.spans.shape[0] > 1


def test_canonical_transcript_cds_differ(genome_dir):
    g = eti_anno.GeneView(source=genome_dir)
    cds = next(
        iter(g.get_features_matching(stable_id="WBGene00004893", biotype="cds")),
    )
    mrna = next(
        iter(g.get_features_matching(stable_id="WBGene00004893", biotype="mrna")),
    )
    assert (cds.spans != mrna.spans).any()


@pytest.mark.parametrize("biotype", ["mrna", "cds"])
def test_get_features_matching_not_canonical(genome_dir, biotype):
    g = eti_anno.GeneView(source=genome_dir)
    features = list(
        g.get_features_matching(
            stable_id="WBGene00004893",
            biotype=biotype,
            canonical=False,
        ),
    )
    assert features


def test_featuredb(genome_dir):
    db = eti_anno.GeneView(source=genome_dir)
    gene = next(iter(db.get_by_stable_id(stable_id="WBGene00000138")))
    cds = next(iter(db.get_feature_children(gene)))
    assert cds.stable_id.startswith("B0019.1")


def test_convert_to_dict():
    import numpy

    raw = {
        "seqid": "I",
        "coord_system_name": "chromosome",
        "start": 11701629,
        "stop": 11703698,
        "spans": numpy.array(
            [
                [11701629, 11701909],
                [11702155, 11702573],
                [11703021, 11703268],
                [11703457, 11703564],
            ],
            dtype=numpy.int32,
        ),
        "strand": 1,
        "name": "T22H2.6a.1",
        "biotype": None,
        "stable_id": "T22H2.6a.1",
        "gene_stable_id": "WBGene00011936",
        "transcript_id": 16967,
    }
    cds = eti_anno.CdsData(**raw)
    got = dict(cds)
    assert isinstance(got, dict)


def test_get_ids_for_biotype(small_install_cfg):
    config = small_install_cfg
    genome = eti_anno.GeneView(
        source=config.installed_genome(species="caenorhabditis_elegans"),
    )
    stable_ids = genome.get_ids_for_biotype(biotype="protein_coding", limit=10)
    assert len(stable_ids) == 10
    assert all(stable_id.startswith("WBGene") for stable_id in stable_ids)
    stable_ids = genome.get_ids_for_biotype(biotype="protein_coding", seqid="I")
    assert len(stable_ids) > 2000


def test_repeat_query(worm_repeats):
    import re

    limit = 10
    repeats = list(
        worm_repeats.get_features_matching(repeat_class="Simple_repeat", limit=limit),
    )
    assert len(repeats) == limit
    rpt = repeats[0]
    got = dict(rpt)
    assert got["xattr"]["repeat_class"] == "Simple_repeat"
    # the order of returned repeats is not guaranteed, so we
    # only check the basic structure
    assert re.search(r"\([A-Z]+\)n", got["xattr"]["repeat_name"])


def test_view_species(worm_db):
    assert worm_db.species == "caenorhabditis_elegans"
    assert worm_db.genes.species == "caenorhabditis_elegans"
    assert worm_db.repeats.species == "caenorhabditis_elegans"
    assert worm_db.biotypes.species == "caenorhabditis_elegans"


@pytest.fixture
def multi_species_db(yeast_db, worm_db):
    dbs = {"caenorhabditis_elegans": worm_db, "saccharomyces_cerevisiae": yeast_db}
    name_map = {
        "worm-WBGene00011936": eti_anno.get_species_seqid(
            species="caenorhabditis_elegans",
            seqid="I",
        ),
        "yeast-YCR105W": eti_anno.get_species_seqid(
            species="saccharomyces_cerevisiae",
            seqid="III",
        ),
    }
    return eti_anno.MultispeciesAnnotations(species_annotations=dbs, name_map=name_map)


def test_multi_species(multi_species_db):
    stable_id = "WBGene00011936"
    kwargs = {"name": stable_id, "biotype": "protein_coding"}
    gene = next(
        iter(
            multi_species_db.get_features_matching(seqid=f"worm-{stable_id}", **kwargs),
        ),
    )
    expect = next(
        iter(
            multi_species_db.species_annotations[
                "caenorhabditis_elegans"
            ].get_features_matching(**kwargs),
        ),
    )
    got = dict(gene)
    got["spans"] = got["spans"].tolist()
    expect = dict(expect)
    expect["spans"] = expect["spans"].tolist()
    assert got == expect


def test_get_ids_for_biotype2(yeast_db):
    features = list(yeast_db.get_ids_for_biotype(biotype="rRNA", limit=10))
    assert len(features) == 10


@pytest.mark.slow
def test_get_ids_for_biotype_seqid(yeast_db, yeast):
    stable_ids = list(
        yeast_db.get_ids_for_biotype(biotype="protein_coding", seqid="III"),
    )
    assert len(stable_ids) == 184  # from direct inspection of sql count distinct
    stable_ids = list(
        yeast_db.get_ids_for_biotype(biotype="protein_coding", seqid=["III", "XVI"]),
    )
    assert len(stable_ids) == 184 + 511  # from direct inspection of sql count distinct
    # make sure the seqid match the input
    seqids = {"III", "XVI"}
    got = {
        r.seqid for stable_id in stable_ids for r in yeast.get_features(name=stable_id)
    }
    assert got == seqids
