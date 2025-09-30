import pytest
from cogent3 import make_tree
from cogent3.core.table import Table

from ensembl_tui import _download as eti_download
from ensembl_tui import _site_map as eti_smap
from ensembl_tui import _species as eti_species


@pytest.fixture(scope="module")
def species():
    return eti_species.make_species_map(species_path=None)


def test_get_name_type(species):
    """should return the (latin|common) name given a latin, common or ensembl
    db prefix names"""
    assert species.get_species_name("human") == "Homo sapiens"
    assert species.get_species_name("homo_sapiens") == "Homo sapiens"
    assert (
        species.get_species_name("canis_lupus_familiaris") == "Canis lupus familiaris"
    )
    assert species.get_common_name("Mus musculus") == "mouse"
    assert species.get_common_name("mus_musculus") == "mouse"


def test_get_species_from_common(species):
    common = "Microbat"
    assert common in species
    got = species.get_species_name(common)
    assert got == "Myotis lucifugus"


def test_get_ensembl_format(species):
    """should take common or latin names and return the corresponding
    ensembl db prefix"""
    assert species.get_ensembl_db_prefix("human") == "homo_sapiens"
    assert species.get_ensembl_db_prefix("mouse") == "mus_musculus"
    assert species.get_ensembl_db_prefix("Mus musculus") == "mus_musculus"
    assert (
        species.get_ensembl_db_prefix("Canis lupus familiaris")
        == "canis_lupus_familiaris"
    )


def test_get_ensembl_prefix_invalid(species):
    assert species.get_ensembl_db_prefix(1) is None
    assert species.get_ensembl_db_prefix("not present") is None


def test_get_genome_name(species):
    got = species.get_genome_name("Sheep - Polled Dorset")
    assert got.startswith("ovis_aries")


def test_get_genome_name_invalid(species):
    assert not species.get_genome_name(1)
    assert not species.get_genome_name("not present")


@pytest.mark.parametrize("arg", ["Human", "human", "homo_sapiens", "Homo sapiens"])
def test_get_abreviation(arg, species):
    table = species.to_table()
    human_identifiers = table.filtered(
        lambda x: x == "human", columns="common_name"
    ).columns.to_dict()
    expect = human_identifiers.pop("abbrev")[0]
    got = species.get_abbreviation(arg)
    assert got == expect
    got = species.get_abbreviation(expect)
    assert got == expect


def test_lookup_raises(species):
    """setting level to raise should create exceptions"""
    with pytest.raises(ValueError):  # noqa: PT011
        species.get_species_name("failme", level="raise")
    with pytest.raises(ValueError):  # noqa: PT011
        species.get_common_name("failme", level="raise")
    with pytest.raises(ValueError):  # noqa: PT011
        species.get_ensembl_db_prefix("failme", level="raise")


def test_lookup_warns(species, capsys):
    """setting level to warn should create warnings"""
    species.get_species_name("failme", level="warn")
    captured = capsys.readouterr()
    assert captured.out.startswith("WARN:")


def test_lookup_latin(species):
    # lonchura_striata_domestica
    query = "lonchura striata"
    got = species.get_abbreviation(query, level="raise")
    assert got == "lon-stri-dome"


def test_to_table(species):
    """returns a table object"""
    table = species.to_table()
    assert isinstance(table, Table)
    assert table.shape[0] > 20
    assert table.shape[1] == len(eti_species.TABLE_COLUMNS)


@pytest.mark.parametrize(
    "name",
    ["Human", "human", "Anas platyrhynchos", "Dog - Basenji"],
)
def test_contains(name, species):
    assert name in species


def test_make_unique_abbrevs():
    names = ["danaus_plexippus", "danaus_plexippus_gca018135715v1"]
    got = eti_species.make_unique_abbrevs(names)
    assert got == dict(zip(names, ["dan-plex", "dan-plex-2"], strict=False))


def test_for_storage(species):
    d = species.for_storage()
    # first value in header should be the genome name
    delim = "\t"
    header = d["header"].split(delim)
    assert header[0] == "genome_name"
    assert set(header) == set(eti_species.TABLE_COLUMNS)
    human_data = dict(
        zip(
            header,
            ["homo_sapiens", *d["homo_sapiens"].split(delim)],
            strict=True,
        )
    )
    expect = {
        "abbrev": "hom-sapi",
        "genome_name": "homo_sapiens",
        "common_name": "human",
        "db_prefix": "homo_sapiens",
    }
    assert human_data == expect


def test_from_storage(species):
    table = species.to_table().sorted(columns="abbrev")
    d = species.for_storage()
    inflated = species.from_storage(d)
    got = inflated.to_table().sorted(columns="abbrev")
    assert got.to_list() == table.to_list()


def test_get_subset(species):
    """should take common or latin names and return the corresponding
    ensembl db prefix"""
    subset = species.get_subset(["human", "Mus musculus", "canis_lupus_familiaris"])
    assert subset.to_table().shape[0] == 3
    got = subset.get_abbreviation("saccharomyces_cerevisiae", level="ignore")
    assert got is None


def test_get_subset_invalid(species):
    """should take common or latin names and return the corresponding
    ensembl db prefix"""
    with pytest.raises(ValueError):  # noqa: PT011
        species.get_subset(["does not exist"])


def test_spec_map_dunder(species):
    assert repr(species)
    assert str(species)
    assert species._repr_html_()


@pytest.mark.internet
@pytest.mark.timeout(10)
def test_species_from_tree(species):
    smap = eti_smap.get_site_map("main")
    ens_tree = eti_download.download_ensembl_tree(
        host=smap.site,
        site_map=smap,
        release="115",
        tree_fname="10_primates_EPO_default.nh",
    )

    genome_tip_map = eti_species.species_from_ensembl_tree(ens_tree, species)
    assert len(genome_tip_map) == len(ens_tree.get_tip_names())


def test_species_from_invalid_tree(species):
    tree = make_tree(tip_names=["a_b_c_d", "e_f_g_h", "i_j_k_l"])
    with pytest.raises(ValueError):  # noqa: PT011
        _ = eti_species.species_from_ensembl_tree(tree, species)
