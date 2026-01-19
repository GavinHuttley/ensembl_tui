import dataclasses
import functools
import io
import pathlib

import duckdb
import numpy

ReturnType = tuple[str, tuple]  # the sql statement and corresponding values


@functools.singledispatch
def array_to_blob(data: numpy.ndarray) -> bytes:
    with io.BytesIO() as out:
        numpy.save(out, data)
        out.seek(0)
        return out.read()


@array_to_blob.register
def _(data: bytes) -> bytes:
    # already a blob
    return data


@functools.singledispatch
def blob_to_array(data: bytes | numpy.ndarray) -> numpy.ndarray:
    with io.BytesIO(data) as out:
        out.seek(0)
        return numpy.load(out)


@blob_to_array.register
def _(data: numpy.ndarray) -> numpy.ndarray:
    return data


class ViewMixin:
    _source: pathlib.Path  # override in subclass

    @property
    def species(self) -> str:
        return self._source.name


@dataclasses.dataclass(slots=True)
class DuckdbParquetBase:
    source: dataclasses.InitVar[pathlib.Path]
    # db is for testing purposes
    db: dataclasses.InitVar[duckdb.DuckDBPyConnection | None] = None
    _source: pathlib.Path = dataclasses.field(init=False)
    _conn: duckdb.DuckDBPyConnection = dataclasses.field(init=False, default=None)  # type: ignore
    _tables: tuple[str, ...] | tuple = ()

    def __post_init__(
        self,
        source: pathlib.Path,
        db: duckdb.DuckDBPyConnection | None,
    ) -> None:
        source = pathlib.Path(source)
        self._source = source
        if db:
            self._conn = db
            return

        self._conn = None

        if not source.exists():
            msg = f"{self._source} does not exist"
            raise FileNotFoundError(msg)
        if not source.is_dir():
            msg = f"{self._source} is not a directory"
            raise OSError(msg)

        if hasattr(self, "_post_init"):
            self._post_init()

    @property
    def conn(self) -> duckdb.DuckDBPyConnection:
        if self._conn is None:
            self._conn = duckdb.connect(":memory:")
            for table in self._tables:
                parquet_file = self._source / f"{table}.parquet"
                if not parquet_file.exists():
                    msg = f"{parquet_file} does not exist"
                    raise FileNotFoundError(msg)

                sql = f"CREATE TABLE {table} AS SELECT * FROM read_parquet('{parquet_file}')"
                self._conn.sql(sql)

        return self._conn

    def __len__(self) -> int:
        return self.num_records()

    def __bool__(self) -> bool:
        # run an efficient check to see if the db is non-empty
        for table in self._tables:
            parquet_file = self._source / f"{table}.parquet"
            if not parquet_file.exists():
                # possibly in memory, so we run a query
                sql = f"SELECT EXISTS(SELECT 1 FROM {table} LIMIT 1)"
                r = self.conn.sql(sql).fetchone()
                if r and r[0]:
                    return True
            elif parquet_file.stat().st_size > 1024:
                # ad hoc threshold of 1kb
                return True

        return False

    def __eq__(self, other: object) -> bool:
        if not isinstance(other, DuckdbParquetBase):
            return False
        return other.conn is self.conn

    @property
    def source(self) -> pathlib.Path:
        return self._source

    def close(self) -> None:
        self.conn.close()

    def num_records(self) -> int:  # pragma: no cover
        # override in subclass
        raise NotImplementedError
