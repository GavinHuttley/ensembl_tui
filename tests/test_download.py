import pytest

import ensembl_tui._config as eti_config
import ensembl_tui._download as eti_download
from ensembl_tui import _mysql_core_attr as eti_db_attr
from ensembl_tui import _site_map as eti_smap


@pytest.mark.internet
@pytest.mark.timeout(10)
def test_get_db_names(tmp_config_just_yeast, ENSEMBL_RELEASE_VERSION):
    cfg = eti_config.read_config(config_path=tmp_config_just_yeast)
    db_names = eti_download.get_core_db_dirnames(cfg, eti_smap.get_site_map("main"))
    assert db_names == {
        "saccharomyces_cerevisiae": f"pub/release-{ENSEMBL_RELEASE_VERSION}/mysql/saccharomyces_cerevisiae_core_{ENSEMBL_RELEASE_VERSION}_4",
    }


def test_make_dumpfiles():
    table_names = set(eti_db_attr.make_mysqldump_names())
    assert "CHECKSUMS" in table_names
    assert {
        n.split(".")[0] for n in table_names - {"CHECKSUMS"}
    } == eti_db_attr.get_all_tables()


@pytest.mark.internet
@pytest.mark.timeout(10)
def test_download_tree():
    from cogent3.core.tree import PhyloNode

    smap = eti_smap.get_site_map("main")
    got = eti_download.download_ensembl_tree(
        host=smap.site,
        site_map=smap,
        release="115",
        tree_fname="10_primates_EPO_default.nh",
    )

    assert isinstance(got, PhyloNode)
