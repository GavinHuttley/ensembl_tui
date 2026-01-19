import configparser
import pathlib

import pytest

from ensembl_tui import _align as eti_align
from ensembl_tui import _config as eti_config
from ensembl_tui import _util as eti_util


def test_installed_genome(default_species_map):
    cfg = eti_config.InstalledConfig(
        release="110",
        install_path="abcd",
        software_versions={},
        species_map=default_species_map,
    )
    assert cfg.installed_genome("human") == pathlib.Path("abcd/genomes/homo_sapiens")


def test_installed_aligns(default_species_map):
    cfg = eti_config.InstalledConfig(
        release="110",
        install_path="abcd",
        software_versions={},
        species_map=default_species_map,
    )
    assert cfg.aligns_path == pathlib.Path("abcd/compara/aligns")


def test_installed_homologies(default_species_map):
    cfg = eti_config.InstalledConfig(
        release="110",
        install_path="abcd",
        software_versions={},
        species_map=default_species_map,
    )
    assert cfg.homologies_path == pathlib.Path("abcd/compara/homologies")


@pytest.fixture
def installed_cfg_path(tmp_config_just_yeast):
    config = eti_config.read_config(config_path=tmp_config_just_yeast)
    return eti_config.write_installed_cfg(config)


def test_read_installed(installed_cfg_path):
    got = eti_config.read_installed_cfg(installed_cfg_path)
    assert str(got.installed_genome("sac-cere")) == str(
        got.install_path / "genomes/saccharomyces_cerevisiae",
    )


def test_read_installed_software_versions(installed_cfg_path):
    import cogent3
    import cogent3_h5seqs

    import ensembl_tui

    config = eti_config.read_installed_cfg(installed_cfg_path)
    assert config.software_versions["ensembl_tui"] == ensembl_tui.__version__
    assert config.software_versions["cogent3"] == cogent3.__version__
    assert config.software_versions["cogent3_h5seqs"] == cogent3_h5seqs.__version__


def test_read_installed_get_version_table(installed_cfg_path):
    import ensembl_tui

    config = eti_config.read_installed_cfg(installed_cfg_path)
    table = config.get_version_table()
    assert table["ensembl_tui", "version"] == ensembl_tui.__version__


def test_installed_config_hash(default_species_map):
    ic = eti_config.InstalledConfig(
        release="11",
        install_path="abcd",
        software_versions={},
        species_map=default_species_map,
    )
    assert hash(ic) == id(ic)
    v = {ic}
    assert len(v) == 1


@pytest.fixture
def installed_aligns(tmp_path, default_species_map):
    align_dir = tmp_path / eti_config._COMPARA_NAME / eti_config._ALIGNS_NAME
    # make two alignment paths with similar names
    names = "10_primates.epo", "24_primates.epo_extended"
    for name in names:
        dirname = align_dir / name
        dirname.mkdir(parents=True, exist_ok=True)
        db = (dirname / f"align_blocks.{eti_align.ALIGN_STORE_SUFFIX}").open(mode="w")
        db.close()
    return eti_config.InstalledConfig(
        release="11",
        install_path=tmp_path,
        software_versions={},
        species_map=default_species_map,
    )


@pytest.fixture
def incomplete_installed(installed_aligns):
    align_path = installed_aligns.aligns_path
    for path in align_path.glob(
        f"*/*.{eti_align.ALIGN_STORE_SUFFIX}",
    ):
        path.unlink()
    return installed_aligns


@pytest.mark.parametrize("pattern", ["10*", "1*prim*", "10_p*", "10_primates.epo"])
def test_get_alignment_path(installed_aligns, pattern):
    got = installed_aligns.path_to_alignment(pattern, eti_align.ALIGN_STORE_SUFFIX)
    assert got.name == "10_primates.epo"


def test_get_alignment_path_incomplete(incomplete_installed):
    with pytest.raises(FileNotFoundError):
        incomplete_installed.path_to_alignment("10*", eti_align.ALIGN_STORE_SUFFIX)


@pytest.mark.parametrize("pattern", ["10pri*", "blah-blah", ""])
def test_get_alignment_path_invalid(installed_aligns, pattern):
    assert (
        installed_aligns.path_to_alignment(pattern, eti_align.ALIGN_STORE_SUFFIX)
        is None
    )


@pytest.mark.parametrize("pattern", ["*pri*", "*epo*"])
def test_get_alignment_path_multiple(installed_aligns, pattern):
    with pytest.raises(ValueError):
        installed_aligns.path_to_alignment(pattern, eti_align.ALIGN_STORE_SUFFIX)


@pytest.fixture
def empty_cfg(tmp_dir, ENSEMBL_RELEASE_VERSION):
    parser = configparser.ConfigParser()
    parser.read(eti_util.get_resource_path("sample.cfg"))
    parser.remove_section("Caenorhabditis elegans")
    parser.remove_section("Saccharomyces cerevisiae")
    parser.remove_section("compara")
    parser.set("local path", "staging_path", value=str(tmp_dir / "staging"))
    parser.set("local path", "install_path", value=str(tmp_dir / "install"))
    parser.set("release", "release", value=ENSEMBL_RELEASE_VERSION)
    return tmp_dir, parser


@pytest.fixture
def cfg_just_aligns(empty_cfg):
    tmp_dir, parser = empty_cfg
    parser.add_section("compara")
    parser.set("compara", "align_names", value="10_primates.epo")
    download_cfg = tmp_dir / "download.cfg"
    with open(download_cfg, "w") as out:
        parser.write(out)

    return download_cfg


COMMON_NAMES = (
    "Crab-eating macaque",
    "Human",
    "Bonobo",
    "Chimpanzee",
    "Mouse Lemur",
    "Gorilla",
    "Gibbon",
    "Vervet-AGM",
    "Sumatran orangutan",
    "Macaque",
)


@pytest.fixture
def cfg_just_genomes(empty_cfg):
    tmp_dir, parser = empty_cfg
    parser.add_section("compara")
    parser.set("compara", "align_names", value="10_primates.epo")
    download_cfg = tmp_dir / "download.cfg"
    for name in COMMON_NAMES:
        parser.add_section(name)
        parser.set(name, "db", value="core")

    with open(download_cfg, "w") as out:
        parser.write(out)

    return download_cfg
    # we make a config using common names


@pytest.mark.internet
@pytest.mark.timeout(10)
def test_read_config_compara_genomes(cfg_just_aligns, default_species_map):
    config = eti_config.read_config(config_path=cfg_just_aligns)
    expected = {default_species_map.get_genome_name(n) for n in COMMON_NAMES}
    assert set(config.species_dbs.keys()) == expected


@pytest.mark.internet
@pytest.mark.timeout(10)
def test_read_config_genomes(cfg_just_genomes, default_species_map):
    config = eti_config.read_config(config_path=cfg_just_genomes)
    expected = {default_species_map.get_genome_name(n) for n in COMMON_NAMES}
    assert set(config.species_dbs.keys()) == expected


def test_read_config_with_domain(tmp_config_domain_format):
    """Test reading config with new 'domain' option"""
    config = eti_config.read_config(config_path=tmp_config_domain_format)
    assert config.domain == "main"  # Domain stores user's choice


def test_read_config_with_host_shows_deprecation(tmp_dir):
    """Test that using 'host' triggers deprecation warning"""
    parser = configparser.ConfigParser()
    parser.read(eti_util.get_resource_path("sample.cfg"))
    # Remove domain, add host (to simulate old config format)
    parser.remove_option("remote path", "domain")
    parser.set("remote path", "host", value="ftp.ensembl.org")
    parser.remove_section("Caenorhabditis elegans")
    parser.remove_section("compara")
    parser.set("local path", "staging_path", value=str(tmp_dir / "staging"))
    parser.set("local path", "install_path", value=str(tmp_dir / "install"))
    cfg_path = tmp_dir / "old_format.cfg"
    with open(cfg_path, "w") as out:
        parser.write(out)

    with pytest.warns(DeprecationWarning, match="'host' option.*deprecated"):
        config = eti_config.read_config(config_path=cfg_path)
    assert config.domain == "ftp.ensembl.org"  # domain set from host


def test_domain_takes_precedence_over_host(tmp_dir):
    """Test that 'domain' takes precedence when both present"""
    parser = configparser.ConfigParser()
    parser.read(eti_util.get_resource_path("sample.cfg"))
    parser.set("remote path", "domain", value="metazoa")
    parser.set("remote path", "host", value="ftp.ensembl.org")
    parser.set("local path", "staging_path", value=str(tmp_dir / "staging"))
    parser.set("local path", "install_path", value=str(tmp_dir / "install"))
    parser.remove_section("Caenorhabditis elegans")
    parser.remove_section("compara")
    cfg_path = tmp_dir / "both.cfg"
    with open(cfg_path, "w") as out:
        parser.write(out)

    with pytest.warns(DeprecationWarning):
        config = eti_config.read_config(config_path=cfg_path)
    assert config.domain == "metazoa"  # domain from config file


def test_invalid_domain_raises(tmp_dir):
    parser = configparser.ConfigParser()
    parser.read(eti_util.get_resource_path("sample.cfg"))
    parser.remove_option("remote path", "host")
    parser.set("remote path", "domain", value="invalid_domain")
    parser.set("local path", "staging_path", value=str(tmp_dir / "staging"))
    parser.set("local path", "install_path", value=str(tmp_dir / "install"))
    parser.remove_section("Caenorhabditis elegans")
    parser.remove_section("compara")
    cfg_path = tmp_dir / "invalid.cfg"
    with open(cfg_path, "w") as out:
        parser.write(out)

    with pytest.raises(ValueError):
        eti_config.read_config(config_path=cfg_path)


def test_missing_both_domain_and_host_raises(tmp_dir):
    """Test that missing both domain and host causes clear error"""
    parser = configparser.ConfigParser()
    parser.read(eti_util.get_resource_path("sample.cfg"))
    parser.remove_option("remote path", "domain")
    parser.set("local path", "staging_path", value=str(tmp_dir / "staging"))
    parser.set("local path", "install_path", value=str(tmp_dir / "install"))
    parser.remove_section("Caenorhabditis elegans")
    parser.remove_section("compara")
    cfg_path = tmp_dir / "missing.cfg"
    with open(cfg_path, "w") as out:
        parser.write(out)

    with pytest.raises(ValueError):
        eti_config.read_config(config_path=cfg_path)


def test_write_config_uses_domain_not_host(tmp_config_domain_format):
    """Test that Config.write() writes 'domain' not 'host'"""
    config = eti_config.read_config(config_path=tmp_config_domain_format)
    config.write()

    # Read back the written config
    parser = configparser.ConfigParser()
    written_path = config.staging_path / eti_config.DOWNLOADED_CONFIG_NAME
    parser.read(written_path)

    # Should have 'domain' and NOT 'host'
    assert parser.has_option("remote path", "domain")
    assert not parser.has_option("remote path", "host")


@pytest.mark.parametrize(
    "domain",
    [
        "main",
        "vertebrates",
        "ftp.ensembl.org",
        "metazoa",
        "ftp.ensemblgenomes.org",
        "protists",
    ],
)
def test_all_registered_domains_work(tmp_dir, domain):
    """Test that all registered domains can be used in config"""
    parser = configparser.ConfigParser()
    parser.add_section("remote path")
    parser.set("remote path", "domain", value=domain)
    parser.add_section("local path")
    parser.set("local path", "staging_path", value=str(tmp_dir / "staging"))
    parser.set("local path", "install_path", value=str(tmp_dir / "install"))
    parser.add_section("release")
    parser.set("release", "release", value="115")
    parser.add_section("Saccharomyces cerevisiae")
    parser.set("Saccharomyces cerevisiae", "db", value="core")

    cfg_path = tmp_dir / f"test_{domain.replace('.', '_')}.cfg"
    with open(cfg_path, "w") as out:
        parser.write(out)

    # Should not raise
    config = eti_config.read_config(config_path=cfg_path)
    assert config.domain == domain  # Domain stores user's choice
