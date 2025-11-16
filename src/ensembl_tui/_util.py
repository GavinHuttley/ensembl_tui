import contextlib
import functools
import os
import pathlib
import re
import shutil
import subprocess
import sys
import tempfile
import typing
import uuid
from collections.abc import Callable, Hashable
from hashlib import md5
from typing import IO

import numba
import numpy
from cogent3.app.composable import define_app
from cogent3.util.parallel import as_completed
from rich import text as rich_text

if typing.TYPE_CHECKING:
    from cogent3.core.table import Table


PathType = str | pathlib.Path | os.PathLike

try:
    from wakepy.keep import running as keep_running

    # trap flaky behaviour on linux
    with keep_running():
        ...

except (NotImplementedError, ImportError):
    keep_running = contextlib.nullcontext

CWD = pathlib.Path.cwd()


def md5sum(data: bytes, *args) -> str:
    """computes MD5SUM

    Notes
    -----
    *args is for signature compatability with checksum
    """
    return md5(data).hexdigest()


# based on https://www.reddit.com/r/learnpython/comments/9bpgjl/implementing_bsd_16bit_checksum/
# and https://www.gnu.org/software/coreutils/manual/html_node/sum-invocation.html#sum-invocation
@numba.jit(nopython=True)
def checksum(data: bytes, size: int) -> tuple[int, int]:  # pragma: no cover
    """computes BSD style checksum"""
    # equivalent to command line BSD sum
    nb = numpy.ceil(size / 1024)
    cksum = 0
    for c in data:
        cksum = (cksum >> 1) + ((cksum & 1) << 15)
        cksum += c
        cksum &= 0xFFFF
    return cksum, int(nb)


def _get_resource_dir() -> PathType:
    """returns path to resource directory"""
    if "ENSEMBLDBRC" in os.environ:
        path = os.environ["ENSEMBLDBRC"]
    else:
        from ensembl_tui import data

        path = pathlib.Path(data.__file__).parent

    path = pathlib.Path(path).expanduser().absolute()
    if not path.exists():
        msg = f"ENSEMBLDBRC directory {str(path)!r} does not exist"
        raise ValueError(msg)

    return pathlib.Path(path)


def get_resource_path(resource: PathType) -> PathType:
    path = ENSEMBLDBRC / resource
    assert path.exists()
    return path


# the following is where essential files live, such as
# the species/common name map and sample download.cfg
ENSEMBLDBRC = _get_resource_dir()


def exec_command(
    cmnd: str,
    stdout: int = subprocess.PIPE,
    stderr: int = subprocess.PIPE,
) -> str | None:
    """executes shell command and returns stdout if completes exit code 0

    Parameters
    ----------

    cmnd : str
      shell command to be executed
    stdout, stderr : streams
      Default value (PIPE) intercepts process output, setting to None
      blocks this."""
    proc = subprocess.Popen(cmnd, shell=True, stdout=stdout, stderr=stderr)
    out, err = proc.communicate()
    if proc.returncode != 0:
        msg = err
        sys.stderr.writelines(f"FAILED: {cmnd}\n{msg}")
        sys.exit(proc.returncode)
    return out.decode("utf8") if out is not None else None


def load_ensembl_checksum(path: pathlib.Path) -> dict:
    """loads the BSD checksums from Ensembl CHECKSUMS file"""
    result = {}
    for line in path.read_text().splitlines():
        line = line.strip()
        if not line:
            continue
        s, b, p, *_ = line.split()
        result[p] = int(s), int(b)
    result.pop("README", None)
    return result


def load_ensembl_md5sum(path: pathlib.Path) -> dict:
    """loads the md5 sum from Ensembl MD5SUM file"""
    result = {}
    for line in path.read_text().splitlines():
        line = line.strip()
        if not line:
            continue
        s, p = line.split()
        result[p] = s
    result.pop("README", None)
    return result


def hash64(data: bytes) -> int:
    """returns 63-bit hash of numpy array as a signed integer"""
    h = md5(data, usedforsecurity=False)
    # Take first 8 bytes and mask to 63 bits (keep it positive)
    hash_val = int.from_bytes(h.digest()[:8], byteorder="big")
    return hash_val & 0x7FFFFFFFFFFFFFFF  # Mask to 63 bits


class atomic_write:  # noqa: N801
    """performs atomic write operations, cleans up if fails"""

    def __init__(
        self,
        path: pathlib.Path,
        tmpdir: str | None = None,
        mode: str = "wb",
        encoding: str | None = None,
    ) -> None:
        """

        Parameters
        ----------
        path
            path to file
        tmpdir
            directory where temporary file will be created
        mode
            file writing mode
        encoding
            text encoding
        """
        path = pathlib.Path(path).expanduser()

        self._path = path
        self._mode = mode
        self._file = None
        self._encoding = encoding
        self._tmppath = self._make_tmppath(tmpdir)

        self.succeeded = None
        self._close_func = self._close_rename_standard

    def _make_tmppath(self, tmpdir: str | None) -> pathlib.Path:
        """returns path of temporary file

        Parameters
        ----------
        tmpdir: Path
            to directory

        Returns
        -------
        full path to a temporary file

        Notes
        -----
        Uses a random uuid as the file name, adds suffixes from path
        """
        suffixes = "".join(self._path.suffixes)
        parent = self._path.parent
        name = f"{uuid.uuid4()}{suffixes}"
        tmpdir: pathlib.Path = (
            pathlib.Path(tempfile.mkdtemp(dir=parent))
            if tmpdir is None
            else pathlib.Path(tmpdir)
        )

        if not tmpdir.exists():
            msg = f"{tmpdir} directory does not exist"
            raise FileNotFoundError(msg)

        return tmpdir / name

    def _get_fileobj(self) -> IO:
        """returns file to be written to"""
        if self._file is None:
            self._file = open(self._tmppath, self._mode)  # noqa: SIM115

        return self._file

    def __enter__(self) -> IO:
        return self._get_fileobj()

    def _close_rename_standard(self, src):
        dest = pathlib.Path(self._path)
        try:
            dest.unlink()
        except FileNotFoundError:
            pass
        finally:
            src.rename(dest)

        shutil.rmtree(src.parent)

    def __exit__(self, exc_type, exc_val, exc_tb):  # noqa: ANN204
        if self._file is not None:
            self._file.close()

        if exc_type is None:
            self._close_func(self._tmppath)
            self.succeeded = True
        else:
            self.succeeded = False

        shutil.rmtree(self._tmppath.parent, ignore_errors=True)

    def write(self, text: str) -> None:
        """writes text to file"""
        fileobj = self._get_fileobj()
        fileobj.write(text)

    def close(self) -> None:
        """closes file"""
        self.__exit__(None, None, None)


_sig_load_funcs = {"CHECKSUMS": load_ensembl_checksum, "MD5SUM": load_ensembl_md5sum}
_sig_calc_funcs = {"CHECKSUMS": checksum, "MD5SUM": md5sum}
_dont_checksum = re.compile("(CHECKSUMS|MD5SUM|README)")
_sig_file = re.compile("(CHECKSUMS|MD5SUM)")


def dont_checksum(path: PathType) -> bool:
    return _dont_checksum.search(str(path)) is not None


@functools.singledispatch
def is_signature(path: pathlib.Path) -> bool:
    return _sig_file.search(path.name) is not None


@is_signature.register
def _(path: str) -> bool:
    return _sig_file.search(path) is not None


@functools.singledispatch
def get_sig_calc_func(sig_path) -> Callable:  # noqa: ANN001
    """returns signature calculating function based on Ensembl path name"""
    msg = f"{type(sig_path)} not supported"
    raise NotImplementedError(msg)


@get_sig_calc_func.register
def _(sig_path: str) -> Callable:
    return _sig_calc_funcs[sig_path]


def get_signature_data(path: pathlib.Path) -> dict:
    return _sig_load_funcs[path.name](path)


def rich_display(c3t: "Table", title_justify: str = "left") -> None:
    """converts a cogent3 Table to a Rich Table and displays it"""
    from rich.console import Console
    from rich.table import Table

    cols = c3t.columns
    columns = []
    for c in c3t.header:
        if tmplt := c3t._column_templates.get(c, None):
            col = [tmplt(v) for v in cols[c]]
        else:
            col = cols[c]
        columns.append(col)

    rich_table = Table(
        title=c3t.title,
        highlight=True,
        title_justify=title_justify,
        title_style="bold blue",
    )
    for col in c3t.header:
        numeric_type = any(v in cols[col].dtype.name for v in ("int", "float"))
        j = "right" if numeric_type else "left"
        rich_table.add_column(col, justify=j, no_wrap=numeric_type)

    for row in zip(*columns, strict=False):
        rich_table.add_row(*row)

    console = Console()
    console.print(rich_table)


_seps = re.compile(r"[-._\s]")


def _name_parts(path: str) -> list[str]:
    return _seps.split(pathlib.Path(path).name.lower())


def _simple_check(align_parts: str, tree_parts: str) -> int:
    """evaluates whether the start of the two paths match"""
    matches = 0
    for a, b in zip(align_parts, tree_parts, strict=False):
        if a != b:
            break
        matches += 1

    return matches


def trees_for_aligns(aligns, trees) -> dict[str, str]:
    aligns = {p: _name_parts(p) for p in aligns}
    trees = {p: _name_parts(p) for p in trees}
    result = {}
    for align, align_parts in aligns.items():
        dists = [
            (_simple_check(align_parts, tree_parts), tree)
            for tree, tree_parts in trees.items()
        ]
        v, p = max(dists)
        if v == 0:
            msg = f"no tree for {align}"
            raise ValueError(msg)

        result[align] = p

    return result


@define_app
def _str_to_bytes(data: str) -> bytes:
    """converts string to bytes"""
    return data.encode("utf8")


@define_app
def _bytes_to_str(data: bytes) -> str:
    """converts bytes into string"""
    return data.decode("utf8")


_biotypes = re.compile(r"(gene|transcript|exon|mRNA|rRNA|protein):")


def sanitise_stableid(stableid: str) -> str:
    """remove <biotype>:E.. from Ensembl stable ID

    Notes
    -----
    The GFF3 files from Ensembl store identifiers as <biotype>:<identifier>,
    this function removes redundant biotype component.
    """
    return _biotypes.sub("", stableid)


_quotes = re.compile(r"^[\'\"]|[\'\"]$")


def strip_quotes(text: str) -> str:
    return _quotes.sub("", text)


def get_iterable_tasks(
    *,
    func: typing.Callable,
    series: typing.Sequence,
    max_workers: int | None,
    **kwargs: dict,
) -> typing.Iterator:
    max_workers = max_workers or 1
    if max_workers == 1:
        return map(func, series)
    return as_completed(func, series, max_workers=max_workers, **kwargs)


# From http://mart.ensembl.org/info/genome/stable_ids/prefixes.html
# The Ensembl stable id structure is
# [species prefix][feature type prefix][a unique eleven digit number]
# feature type prefixes are
# E exon
# FM Ensembl protein family
# G gene
# GT gene tree
# P protein
# R regulatory feature
# T transcript
_feature_type_1 = {"E", "G", "P", "R", "T"}
_feature_type_2 = {"FM", "GT"}


def get_stableid_prefix(stableid: str) -> str:
    """returns the prefix component of a stableid"""
    if len(stableid) < 15:
        msg = f"{stableid!r} too short"
        raise ValueError(msg)

    if stableid[-13:-11] in _feature_type_2:
        return stableid[:-13]
    if stableid[-12] not in _feature_type_1:
        msg = f"{stableid!r} has unknown feature type {stableid[-13]!r}"
        raise ValueError(msg)
    return stableid[:-12]


class _printer:  # noqa: N801
    from rich.console import Console

    def __init__(self) -> None:
        self._console = self.Console()

    def __call__(self, text: str, colour: str, style: str = "") -> None:
        """print text in colour"""
        msg = rich_text.Text.from_markup(text, style=style)
        msg.stylize(colour)
        self._console.print(msg)


print_colour = _printer()


class unique_value_indexer:  # noqa: N801
    """creates indexes for unique values

    Notes
    -----
    Instance is callable and will return the unique index for a value.
    Indexes are determined by the order of first appearance.
    """

    __slots__ = ("_values",)

    def __init__(self) -> None:
        self._values = {}

    def __call__(self, value: Hashable) -> int:
        if not (index := self._values.get(value, 0)):
            index = len(self._values) + 1
            self._values[value] = index
        return index

    def __iter__(self) -> typing.Iterator[tuple[int, Hashable]]:
        for value, index in self._values.items():
            yield index, value


class category_indexer:  # noqa: N801
    """maps multiple keys to the same index value for distinct categories"""

    __slots__ = ("_index", "_values")

    def __init__(self) -> None:
        self._values = {}
        self._index = 0

    def __call__(self, category: str | int, vals: set[str] | set[int]) -> int:
        if not vals:
            message = f"no values provided for {category}"
            raise ValueError(message)
        if (entries := self._values.get(category)) and (
            shared := entries.keys() & vals
        ):
            key = next(iter(shared))
            index = entries[key]
            new_vals = vals - entries.keys()
        else:
            entries = self._values.get(category, {})
            self._index += 1
            index = self._index
            new_vals = vals
        entries |= dict.fromkeys(new_vals, index)
        self._values[category] = entries
        return index

    def __iter__(self) -> typing.Iterator[tuple[int, str | int, str | int]]:
        for cat_name, cat in self._values.items():
            for value, index in cat.items():
                yield index, cat_name, value


@contextlib.contextmanager
def tempdir(working_dir: pathlib.Path | str | None = None) -> pathlib.Path:
    """context manager returns a temporary directory in working_dir"""
    with tempfile.TemporaryDirectory(dir=working_dir) as temp_dir:
        yield pathlib.Path(temp_dir)


def make_column_constant(schema: tuple[str, ...]) -> tuple[str, ...]:
    return tuple(c.split()[0] for c in schema)


_has_wildcard = re.compile(r"[*?\[\]]")


def contains_glob_pattern(s: str) -> bool:
    return _has_wildcard.search(s) is not None
