import pytest
from cogent3.core.table import Table

from ensembl_tui._species import Species


def test_get_name_type():
    """should return the (latin|common) name given a latin, common or ensembl
    db prefix names"""
    assert Species.get_species_name("human") == "Homo sapiens"
    assert Species.get_species_name("homo_sapiens") == "Homo sapiens"
    assert (
        Species.get_species_name("canis_lupus_familiaris") == "Canis lupus familiaris"
    )
    assert Species.get_common_name("Mus musculus") == "Mouse"
    assert Species.get_common_name("mus_musculus") == "Mouse"


def test_get_ensembl_format():
    """should take common or latin names and return the corresponding
    ensembl db prefix"""
    assert Species.get_ensembl_db_prefix("human") == "homo_sapiens"
    assert Species.get_ensembl_db_prefix("mouse") == "mus_musculus"
    assert Species.get_ensembl_db_prefix("Mus musculus") == "mus_musculus"
    assert (
        Species.get_ensembl_db_prefix("Canis lupus familiaris")
        == "canis_lupus_familiaris"
    )


def test_add_new_species():
    """should correctly add a new species/common combination and infer the
    correct ensembl prefix"""
    species_name, common_name = "Otolemur garnettii", "Bushbaby"
    Species.amend_species(species_name, common_name)
    assert Species.get_species_name(species_name) == species_name
    assert Species.get_species_name("Bushbaby") == species_name
    assert Species.get_species_name(common_name) == species_name
    assert Species.get_common_name(species_name) == common_name
    assert Species.get_common_name("Bushbaby") == common_name
    assert Species.get_ensembl_db_prefix("Bushbaby") == "otolemur_garnettii"
    assert Species.get_ensembl_db_prefix(species_name) == "otolemur_garnettii"
    assert Species.get_ensembl_db_prefix(common_name) == "otolemur_garnettii"


def test_amend_existing():
    """should correctly amend an existing species"""
    species_name = "Ochotona princeps"
    common_name1 = "american pika"
    common_name2 = "pika"
    ensembl_pref = "ochotona_princeps"
    Species.amend_species(species_name, common_name1)
    assert Species.get_common_name(species_name) == common_name1
    Species.amend_species(species_name, common_name2)
    assert Species.get_species_name(common_name2) == species_name
    assert Species.get_species_name(ensembl_pref) == species_name
    assert Species.get_common_name(species_name) == common_name2
    assert Species.get_common_name(ensembl_pref) == common_name2
    assert Species.get_ensembl_db_prefix(species_name) == ensembl_pref
    assert Species.get_ensembl_db_prefix(common_name2) == ensembl_pref


def test_lookup_raises():
    """setting level  to raise should create exceptions"""
    with pytest.raises(ValueError):
        Species.get_species_name("failme", level="raise")
    with pytest.raises(ValueError):
        Species.get_common_name("failme", level="raise")
    with pytest.raises(ValueError):
        Species.get_ensembl_db_prefix("failme")


def test_to_table():
    """returns a table object"""
    table = Species.to_table()
    assert isinstance(table, Table)
    assert table.shape[0] > 20
    assert table.shape[1] == 3


@pytest.mark.parametrize(
    "name",
    ("Human", "human", "Anas platyrhynchos", "Dog - Basenji"),
)
def test_contains(name):
    assert name in Species


@pytest.mark.parametrize("feature_prefix", ("E", "FM", "G", "GT", "P", "R", "T"))
def test_db_prefix_from_stableid(feature_prefix):
    stableid = f"ENSPTI{feature_prefix}00000025624"
    expect = "panthera_tigris_altaica"
    got = Species.get_db_prefix_from_stableid(stableid)
    assert got == expect


@pytest.mark.parametrize("invalid", ("blah", "ENSPTIX012345678911"))
def test_db_prefix_from_stableid_fail(invalid):
    with pytest.raises(ValueError):
        Species.get_db_prefix_from_stableid(invalid)


@pytest.mark.parametrize(
    "stableid,expect",
    (
        ("ENSCSAG00000025624", "chlorocebus_sabaeus"),
        ("ENSGGOG00000041075", "gorilla_gorilla"),
        ("ENSPPAG00000009915", "pan_paniscus"),
        ("ENSG00000278672", "homo_sapiens"),
    ),
)
def test_db_prefixes_from_stablesids(stableid, expect):
    got = Species.get_db_prefix_from_stableid(stableid)
    assert got == expect


def test_make_unique_abbrevs():
    names = ["danaus_plexippus", "danaus_plexippus_gca018135715v1"]
    got = eti_species.make_unique_abbrevs(names)
    assert got == dict(zip(names, ["dan-plex", "dan-plex-2"], strict=False))
