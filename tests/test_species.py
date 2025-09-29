import pytest
from cogent3.core.table import Table

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


def test_get_genome_name(species):
    got = species.get_genome_name("Sheep - Polled Dorset")
    assert got.startswith("ovis_aries")


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
