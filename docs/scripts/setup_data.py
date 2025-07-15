import pathlib
import shutil
import urllib.request
import zipfile

root_dir = pathlib.Path(__file__)
while not (root_dir / "data").exists():
    root_dir = root_dir.parent

DATA_DIR = root_dir / "data"

SMALL_DATA_URL = "https://www.dropbox.com/scl/fi/a3dkt04z7d1t2p3io1pp6/small-114.zip?rlkey=di9ty6diu1kusjsopam891zyg&dl=1"
SMALL_DATA_DIRNAME = "small-114"

APES_DATA_URL = "https://www.dropbox.com/scl/fi/cyr1p5aqteffsggtlqjo7/apes-114.zip?rlkey=sbq1h0kx37fz7gsmlblherxr5&dl=1"
APES_DATA_DIRNAME = "apes-114"


def setup_installed(url: str, dest: str) -> str:
    zip_dest = DATA_DIR / f"{dest}.zip"
    expected = DATA_DIR / dest
    if zip_dest.exists() and expected.exists():
        shutil.rmtree(expected)
    elif not zip_dest.exists():
        urllib.request.urlretrieve(url, filename=zip_dest)  # noqa: S310

    with zipfile.ZipFile(zip_dest, "r") as zip_ref:
        zip_ref.extractall(DATA_DIR)

    return dest


def on_pre_build(*args, **kwargs) -> None:
    for url, dirname in [
        (SMALL_DATA_URL, SMALL_DATA_DIRNAME),
        (APES_DATA_URL, APES_DATA_DIRNAME),
    ]:
        setup_installed(url, dirname)


if __name__ == "__main__":
    on_pre_build()
