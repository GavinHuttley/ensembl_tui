import pytest

from ensembl_tui import _site_map as eti_smap


@pytest.mark.parametrize("site", ["ftp.ensembl.org"])
def test_correct_site(site):
    smap = eti_smap.get_site_map(site)
    assert smap.site == site


def test_standard_smp():
    sm = eti_smap.get_site_map("ftp.ensembl.org")
    assert sm.get_seqs_path("abcd", collection_name=None) == "fasta/abcd/dna"
    assert sm.get_annotations_path("abcd") == "mysql/abcd"


def test_get_site_map_names():
    names = eti_smap.get_site_map_names()
    assert "main" in names
    assert "metazoa" in names
    assert "protists" in names


def test_get_default_site_map_species():
    smap = eti_smap.get_site_map("main")
    assert smap.species_file_name == "species_EnsemblVertebrates.txt"


@pytest.mark.parametrize("site", ["main", "metazoa"])
def test_get_remote_path(site):
    smap = eti_smap.get_site_map(site)
    release = "62"
    prefix = "pub" if site == "main" else f"pub/{site}"
    expect = f"{prefix}/release-{release}"
    assert smap.get_remote_release_path(release) == expect
