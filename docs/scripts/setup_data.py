import pathlib
import shutil
import urllib.request
import zipfile

root_dir = pathlib.Path(__file__)
while not (root_dir / "docs").exists():
    root_dir = root_dir.parent

ROOT_DIR = root_dir / "docs"

DATA_URL = "https://www.dropbox.com/scl/fi/wgh6yrgh0hb4qrd442q13/ensembl_tui_doc_data.zip?rlkey=jc0n0bd1p92xbdg5s2tw5e5ci&dl=1"


def cleanup_data() -> None:
    for dirname in (
        "demo",
        "apes-115",
        "small-download",
        "human_data",
        "apes_homologs",
        "apes_aligns",
        "worm_yeast",
        "worm",
        "yeast",
    ):
        temp_dir = ROOT_DIR / dirname
        shutil.rmtree(temp_dir, ignore_errors=True)


def setup_installed(url: str, dest_zip: str, dest: str) -> str:
    zip_dest = ROOT_DIR / dest_zip
    unzipped_dest = ROOT_DIR / zip_dest.stem
    expected = ROOT_DIR / dest
    if zip_dest.exists() and expected.exists():
        # we will inflate zip archive each time
        shutil.rmtree(expected)
    elif not zip_dest.exists():
        urllib.request.urlretrieve(url, filename=zip_dest)  # noqa: S310

    with zipfile.ZipFile(zip_dest, "r") as zip_ref:
        zip_ref.extractall(ROOT_DIR)

    unzipped_dest.rename(expected)
    return dest


def setup_fonts() -> None:
    font_dir = ROOT_DIR / "overrides/assets/fonts"
    font_dir.mkdir(parents=True, exist_ok=True)
    # Download and install fonts here
    font_dest = font_dir / "roboto-fonts.zip"
    if not font_dest.exists():
        font_url = "https://fonts.google.com/download?family=Roboto"
        urllib.request.urlretrieve(font_url, filename=font_dest)  # noqa: S310


def on_pre_build(*args, **kwargs) -> None:
    cleanup_data()
    setup_fonts()
    demo = ROOT_DIR / "demo"
    shutil.rmtree(demo, ignore_errors=True)

    setup_installed(DATA_URL, "ensembl_tui_doc_data.zip", "data")


def on_post_build(*args, **kwargs) -> None:
    """Clean up temporary data files after build."""
    cleanup_data()


if __name__ == "__main__":
    on_pre_build()
