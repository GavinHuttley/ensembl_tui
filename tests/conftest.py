import pathlib
import shutil
from configparser import ConfigParser

import pytest

from ensembl_tui import _align as eti_align
from ensembl_tui import _config as eti_config
from ensembl_tui import _genome as eti_genome
from ensembl_tui import _util as eti_util


@pytest.fixture(scope="session")
def DATA_DIR():
    return pathlib.Path(__file__).parent / "data"


@pytest.fixture(scope="session")
def ENSEMBL_RELEASE_VERSION() -> str:
    return "114"


@pytest.fixture
def tmp_dir(tmp_path_factory):
    return tmp_path_factory.mktemp("cli")


@pytest.fixture
def tmp_config(tmp_dir):
    # create a simpler download config
    # we want a very small test set
    parser = ConfigParser()
    parser.read(eti_util.get_resource_path("sample.cfg"))
    parser.remove_section("Caenorhabditis elegans")
    parser.remove_section("compara")
    parser.set("local path", "staging_path", value=str(tmp_dir / "staging"))
    parser.set("local path", "install_path", value=str(tmp_dir / "install"))
    download_cfg = tmp_dir / "download.cfg"
    with open(download_cfg, "w") as out:
        parser.write(out)

    return download_cfg


def name_as_seqid(species, seqid, start, end):  # noqa: ARG001
    return seqid


@pytest.fixture
def namer():
    return name_as_seqid


TEST_DATA_URL = "https://www.dropbox.com/scl/fi/a3dkt04z7d1t2p3io1pp6/small-114.zip?rlkey=di9ty6diu1kusjsopam891zyg&dl=1"
SMALL_DATA_DIRNAME = "small-114"


@pytest.fixture(scope="session")
def small_path(DATA_DIR):
    import urllib
    import zipfile

    small_data_path = DATA_DIR / SMALL_DATA_DIRNAME
    if not small_data_path.exists():
        dest = DATA_DIR / f"{SMALL_DATA_DIRNAME}.zip"
        urllib.request.urlretrieve(TEST_DATA_URL, dest)
        with zipfile.ZipFile(dest, "r") as zip_ref:
            zip_ref.extractall(DATA_DIR)
    return small_data_path


@pytest.fixture(scope="session")
def small_download_path(small_path):
    return small_path / "download"


@pytest.fixture(scope="session")
def small_install_path(small_path):
    return small_path / "install"


@pytest.fixture(scope="session")
def small_download_cfg(small_download_path):
    return eti_config.read_config(small_download_path)


@pytest.fixture(scope="session")
def small_install_cfg(small_install_path):
    return eti_config.read_installed_cfg(small_install_path)


@pytest.fixture(scope="module")
def tmp_downloaded(tmp_path_factory, small_download_path):
    tmp_path = tmp_path_factory.mktemp("downloaded")
    dest = tmp_path / small_download_path.name
    shutil.copytree(small_download_path, dest)
    return dest


@pytest.fixture(scope="session")
def yeast(small_install_cfg):
    return eti_genome.load_genome(
        config=small_install_cfg,
        species="saccharomyces_cerevisiae",
    )


@pytest.fixture
def yeast_db(yeast):
    return yeast.annotation_db


@pytest.fixture
def yeast_seqs(yeast):
    return yeast.seqs


@pytest.fixture
def yeast_genes(yeast_db):
    return yeast_db.genes


@pytest.fixture
def yeast_biotypes(yeast_db):
    return yeast_db.biotypes


@pytest.fixture
def yeast_repeats(yeast_db):
    return yeast_db.repeats


@pytest.fixture(scope="session")
def worm(small_install_cfg):
    return eti_genome.load_genome(
        config=small_install_cfg,
        species="caenorhabditis_elegans",
    )


@pytest.fixture
def worm_db(worm):
    return worm.annotation_db


@pytest.fixture
def worm_genes(worm_db):
    return worm_db.genes


@pytest.fixture
def worm_biotypes(worm_db):
    return worm_db.biotypes


@pytest.fixture
def worm_repeats(worm_db):
    return worm_db.repeats


@pytest.fixture
def genome_dir(small_install_cfg):
    return small_install_cfg.installed_genome(species="caenorhabditis_elegans")


@pytest.fixture(scope="module")
def tmp_config_no_compara(tmp_path_factory, small_download_path):
    # modify the small_download_path to remove the compara folder
    tmp_path = tmp_path_factory.mktemp("downloaded")
    dest = tmp_path / small_download_path.name
    shutil.copytree(small_download_path, dest)
    shutil.rmtree(dest / "genomes" / "saccharomyces_cerevisiae")
    shutil.rmtree(dest / "compara")

    # create a config with one specie but without compara section
    parser = ConfigParser()
    parser.read(next(iter(small_download_path.glob("*cfg"))))
    parser.remove_section("saccharomyces_cerevisiae")
    parser.remove_section("compara")

    download_cfg = dest / "downloaded.cfg"
    with open(download_cfg, "w") as out:
        parser.write(out)

    return dest


TEST_APES_DATA_URL = "https://www.dropbox.com/scl/fi/cyr1p5aqteffsggtlqjo7/apes-114.zip?rlkey=sbq1h0kx37fz7gsmlblherxr5&dl=1"
APES_DATA_DIRNAME = "apes-114"


def apes_install(data_url, data_dir, data_name):
    import urllib
    import zipfile

    data_path = data_dir / data_name
    if not data_path.exists():
        dest = data_dir / f"{data_name}.zip"
        urllib.request.urlretrieve(data_url, dest)  # noqa: S310
        with zipfile.ZipFile(dest, "r") as zip_ref:
            zip_ref.extractall(data_dir)
    return data_path


@pytest.fixture(scope="session")
def apes_install_path(DATA_DIR):
    return apes_install(TEST_DATA_URL, DATA_DIR, APES_DATA_DIRNAME)


TEST_APES_MAF_URL = "https://www.dropbox.com/scl/fi/9kc57hitzhwwifq35je8l/apes-114-maf.zip?rlkey=mxeytmuv672cpfar7iirh7emm&dl=1"
APES_MAF_DIRNAME = "apes-114-maf"


@pytest.fixture(scope="session")
def apes_maf_install_path(DATA_DIR):
    path = apes_install(TEST_APES_MAF_URL, DATA_DIR, APES_MAF_DIRNAME)
    return next(iter(path.glob("*.maf.gz")))


@pytest.fixture
def apes(apes_install_path):
    config = eti_config.read_installed_cfg(apes_install_path)
    return {
        sp: eti_genome.load_genome(config=config, species=sp)
        for sp in config.list_genomes()
    }


@pytest.fixture
def apes_aligndb(apes_install_path):
    config = eti_config.read_installed_cfg(apes_install_path)
    return eti_align.load_aligndb(config=config, align_name="primate")
