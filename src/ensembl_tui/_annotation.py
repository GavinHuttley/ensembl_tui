import dataclasses
import functools
import pathlib
import types
import typing

import cogent3
import duckdb
import numpy
from cogent3.core.annotation_db import (
    AnnotationDbABC,
    FeatureDataType,
    SqliteAnnotationDbMixin,
)

import ensembl_tui._mysql_core_attr as core_tables
from ensembl_tui import _storage_mixin as eti_storage

if typing.TYPE_CHECKING:
    from cogent3.core.table import Table

OptInt = int | None
OptStr = str | None
OptBool = bool | None
StrOrBool = str | bool
FeatureDictVals = str | int | numpy.ndarray

GENE_ATTR_COLUMNS = core_tables.GENE_ATTR_COLUMNS


class FeatureDataMixin:  # supports getitem as a dict on properties
    def __getitem__(self, key: str) -> FeatureDictVals:
        return getattr(self, key)

    def __setitem__(self, key: str, value: FeatureDictVals) -> None:
        setattr(self, key, value)

    def pop(self, key: str) -> FeatureDictVals:
        value = getattr(self, key)
        delattr(self, key)
        return value

    def get(self, key: str, default: typing.Any = None) -> FeatureDictVals:
        return getattr(self, key, default)


@dataclasses.dataclass(slots=True)
class FeatureDataBase(FeatureDataMixin):
    seqid: str = dataclasses.field(kw_only=True)
    coord_system_name: str = dataclasses.field(kw_only=True)
    start: int = dataclasses.field(kw_only=True)
    stop: int = dataclasses.field(kw_only=True)
    spans: numpy.ndarray[numpy.int32] = dataclasses.field(kw_only=True)
    strand: int = dataclasses.field(kw_only=True)
    name: str = dataclasses.field(kw_only=True)
    biotype: str | None = dataclasses.field(kw_only=True, default=None)

    def __dict__(self) -> dict:
        return dataclasses.asdict(self)

    def __iter__(
        self,
    ) -> typing.Iterator[tuple[str, str | int | numpy.ndarray | dict | None]]:
        feature_fields = {"seqid", "biotype", "name", "spans", "strand"}
        xattr = {}
        seen_name = False
        for field in dataclasses.fields(self):
            if field.name in {"start", "stop"}:
                # start and stop are not to be included
                # as they can be derived from spans
                continue

            if field.name not in feature_fields:
                # other fields are considered xattr
                xattr[field.name] = getattr(self, field.name)
                continue

            yield field.name, getattr(self, field.name)
            if field.name == "name":
                seen_name = True

        if not seen_name:
            yield "name", self.name

        if xattr:
            yield "xattr", xattr


@dataclasses.dataclass(slots=True)
class GeneData(FeatureDataBase):
    canonical_transcript_id: int
    stable_id: str
    gene_id: str
    symbol: str
    name: str | None = dataclasses.field(kw_only=True, default=None)
    description: str | None = dataclasses.field(kw_only=True, default=None)

    def __post_init__(self) -> None:
        self.name = self.stable_id


@dataclasses.dataclass(slots=True)
class TranscriptData(FeatureDataBase):
    transcript_id: int
    stable_id: str
    gene_stable_id: str
    name: str | None = dataclasses.field(kw_only=True, default=None)
    gene_id: int | None = dataclasses.field(kw_only=True, default=None)
    symbol: str | None = dataclasses.field(kw_only=True, default=None)

    def __post_init__(self) -> None:
        self.name = self.stable_id


@dataclasses.dataclass(slots=True)
class CdsData(FeatureDataBase):
    stable_id: str
    gene_stable_id: str
    name: str | None = dataclasses.field(kw_only=True, default=None)
    transcript_id: int | None = dataclasses.field(kw_only=True, default=None)

    def __post_init__(self) -> None:
        self.name = self.stable_id


@dataclasses.dataclass(slots=True)
class RepeatData(FeatureDataBase):
    repeat_type: str
    repeat_class: str
    repeat_name: str


def _matching_conditions(
    equals_conds: dict[str, str | int] | None = None,
    like_conds: dict[str, str] | None = None,
    allow_partial: bool = True,
) -> str:
    """creates WHERE clause

    Parameters
    ----------
    equals_conds
        column name and values to be matched by equals
    like_conds
        column name and values to be matched by ILIKE (case-insensitive)
    allow_partial
        if False, only records within start, stop are included. If True,
        all records that overlap the segment defined by start, stop are included.

    Returns
    -------
    str, tuple
        the SQL statement and the tuple of values
    """
    equals_conds = equals_conds or {}
    start = equals_conds.pop("start", None)
    stop = equals_conds.pop("stop", None)

    sql = []
    if equals_conds:
        conds = []
        for col, val in equals_conds.items():
            # conditions are filtered for None before here, so we should add
            # an else where the op is assigned !=
            if isinstance(val, tuple | set | list):
                vals = ",".join(f"{v!r}" for v in val)
                conds.append(f"{col} IN ({vals})")
            elif val is not None:
                conds.append(f"{col} = {val!r}")

        sql.append(" AND ".join(conds))
    if like_conds:
        sql.extend(f"{col} ILIKE '%{val}%'" for col, val in like_conds.items())
    if start is not None and stop is not None:
        if allow_partial:
            # allow matches that overlap the segment
            cond = [
                f"(start >= {start} AND stop <= {stop})",  # lies within the segment
                f"(start <= {start} AND stop > {start})",  # straddles beginning of segment
                f"(start < {stop} AND stop >= {stop})",  # straddles stop of segment
                f"(start <= {start} AND stop >= {stop})",  # includes segment
            ]
            cond = " OR ".join(cond)
        else:
            # only matches within bounds
            cond = f"start >= {start} AND stop <= {stop}"
        sql.append(f"({cond})")
    elif start is not None:
        # if query has no stop, then any feature containing start
        cond = f"(start <= {start} AND {start} < stop)"
        sql.append(f"({cond})")
    elif stop is not None:
        # if query has no start, then any feature containing stop
        cond = f"(start <= {stop} AND {stop} < stop)"
        sql.append(f"({cond})")

    return " AND ".join(sql)


def _select_records_sql(
    *,
    table_name: str,
    equals_conds: dict[str, str | int] | None = None,
    like_conds: dict[str, str] | None = None,
    columns: typing.Sequence[str] | None = None,
    allow_partial: bool = True,
) -> str:
    """create SQL select statement and values

    Parameters
    ----------
    table_name
        containing the data to be selected from
    columns
        values to select
    equals_conds
        the WHERE condition = value
    like_conds
        the WHERE condition ILIKE '%value%'
    start, stop
        select records whose (start, stop) values lie between start and stop,
        or overlap them if (allow_partial is True)
    allow_partial
        if False, only records within start, stop are included. If True,
        all records that overlap the segment defined by start, stop are included.

    Returns
    -------
    str, tuple
        the SQL statement and the tuple of values
    """

    conditions = _matching_conditions(
        equals_conds=equals_conds,
        like_conds=like_conds,
        allow_partial=allow_partial,
    )
    cols = "*" if columns is None else f"{', '.join(columns)}"
    sql = f"SELECT {cols} FROM {table_name}"
    return f"{sql} WHERE {conditions}" if conditions else sql


def is_derived_biotype(biotype: str | None) -> bool:
    """returns True if the biotype is derived from a transcript_attr record"""
    biotype = biotype or ""
    return biotype.lower() in _derived_biotypes


def gene_from_gene_record(record: dict) -> GeneData:
    """returns a GeneData record from a gene record"""
    start, stop = (
        record.get("start"),
        record.get("stop"),
    )
    record["spans"] = numpy.array([sorted([start, stop])], dtype=numpy.int32)  # type: ignore
    return GeneData(**record)


def cds_from_gene_record(transcript: dict) -> CdsData:
    """returns a cds record from a transcript_attr record"""
    if not (spans := transcript.pop("cds_spans", None)):
        msg = f"No CDS spans found for {transcript['cds_stable_id']=!r}"
        raise ValueError(msg)

    spans = eti_storage.blob_to_array(spans)
    stable_id = transcript.pop("cds_stable_id")
    gene_stable_id = transcript.pop("gene_stable_id")
    transcript.pop("transcript_stable_id", None)
    transcript.pop("transcript_spans", None)
    return CdsData(
        **{
            **transcript,
            "spans": spans,
            "stable_id": stable_id,
            "gene_stable_id": gene_stable_id,
        },
    )


def transcript_from_gene_record(transcript: dict) -> TranscriptData:
    """returns a transcript record from a transcript_attr record"""
    if not (spans := transcript.pop("transcript_spans", None)):
        msg = f"No transcript spans found for {transcript['transcript_stable_id']=!r}"
        raise ValueError(msg)

    spans = eti_storage.blob_to_array(spans)
    stable_id = transcript.pop("transcript_stable_id")
    gene_stable_id = transcript.pop("gene_stable_id")
    transcript.pop("cds_stable_id", None)
    transcript.pop("cds_spans", None)
    return TranscriptData(
        **{
            **transcript,
            "spans": spans,
            "stable_id": stable_id,
            "gene_stable_id": gene_stable_id,
        },
    )


# make the following module level dicts immutable by using the mapping proxy
## maps derived biotype to the ensembl biotype
_derived_biotypes = types.MappingProxyType(
    {
        "gene": ("protein_coding",),
        "cds": ("protein_coding",),
        "mrna": ("protein_coding",),
        "transcript": ("protein_coding",),
    },
)

# maps derived biotype to the function for creating the derived data instance
_derived_biotype_funcs = types.MappingProxyType(
    {
        "gene": gene_from_gene_record,
        "cds": cds_from_gene_record,
        "mrna": transcript_from_gene_record,
        "transcript": transcript_from_gene_record,
    },
)


@dataclasses.dataclass
class BiotypeView(eti_storage.DuckdbParquetBase, eti_storage.ViewMixin):
    _tables: tuple[str] = ("gene_attr",)

    def num_records(self) -> int:
        """returns the number of distinct biotypes"""
        return len(self.distinct)

    @functools.cached_property
    def distinct(self) -> tuple[str, ...]:
        sql = f"SELECT DISTINCT biotype FROM {self._tables[0]}"
        return tuple(r[0] for r in self.conn.sql(sql).fetchall())

    def count_distinct(self) -> "Table":
        sql = (
            f"SELECT biotype, COUNT(*) AS freq FROM {self._tables[0]} GROUP BY biotype"
        )
        got = self.conn.sql(sql).fetchall()
        return cogent3.make_table(
            header=["biotype", "count"],
            data=got,
            index_name="biotype",
        )


@dataclasses.dataclass
class GeneView(eti_storage.DuckdbParquetBase, eti_storage.ViewMixin):
    _tables: tuple[str, str] = ("gene_attr", "transcript_attr")

    def num_records(self) -> int:
        """returns the number of distinct genes as identified by stable_id"""
        sql = "SELECT COUNT(DISTINCT stable_id) FROM gene_attr"
        return self.conn.sql(sql).fetchone()[0]

    def _get_features_matching(
        self,
        *,
        seqid: OptStr = None,
        biotype: OptStr = None,
        stable_id: OptStr = None,
        start: OptInt = None,
        stop: OptInt = None,
        strand: OptStr = None,
        symbol: OptStr = None,
        description: OptStr = None,
        **kwargs,  # noqa: ANN003
    ) -> typing.Iterator[dict]:
        # add supoport for querying by symbol and description
        stable_id = stable_id or kwargs.pop("name", None)
        limit = kwargs.pop("limit", None)
        local_vars = locals()
        if kwargs := {
            k: v
            for k, v in local_vars.items()
            if k not in ("self", "kwargs", "limit") and v is not None
        }:
            like_conds = (
                {"description": kwargs.pop("description")} if description else None
            )
            sql = _select_records_sql(
                equals_conds=kwargs,
                like_conds=like_conds,
                table_name="gene_attr",
                columns=GENE_ATTR_COLUMNS,
            )
        else:
            sql = f"SELECT {','.join(GENE_ATTR_COLUMNS)} FROM gene_attr"

        sql += f" LIMIT {limit}" if limit else ""
        for record in self.conn.sql(sql).fetchall():
            yield dict(zip(GENE_ATTR_COLUMNS, record, strict=True))

    def _transcript_from_gene(
        self,
        *,
        gene_field_name: str,
        transcript_field_name: str,
        **kwargs,
    ) -> dict:
        columns = (
            "transcript_id",
            "seqid",
            "coord_system_name",
            "start",
            "stop",
            "strand",
            "transcript_spans",
            "transcript_stable_id",
            "cds_spans",
            "cds_stable_id",
        )
        sql = f"SELECT {','.join(columns)} FROM transcript_attr WHERE {transcript_field_name} = ?"
        for record in self._get_features_matching(**kwargs):
            field_value = record.get(gene_field_name)
            gene_stable_id = {"gene_stable_id": record.get("stable_id")}
            for tr_record in self.conn.sql(sql, params=(field_value,)).fetchall():
                yield dict(zip(columns, tr_record, strict=True)) | gene_stable_id

    def get_features_matching(
        self,
        *,
        seqid: OptStr = None,
        biotype: OptStr = None,
        stable_id: OptStr = None,
        start: OptInt = None,
        stop: OptInt = None,
        strand: OptStr = None,
        symbol: OptStr = None,
        description: OptStr = None,
        canonical: bool = True,
        **kwargs,  # noqa: ANN003
    ) -> typing.Iterator[GeneData | TranscriptData | CdsData]:
        local_vars = locals()
        local_vars = {
            k: v
            for k, v in local_vars.items()
            if k not in ("self", "kwargs", "local_vars", "canonical") and v is not None
        }
        kwargs |= local_vars
        if not is_derived_biotype(kwargs.get("biotype")):
            for record in self._get_features_matching(**kwargs):
                yield gene_from_gene_record(record)
            return

        if biotype == "gene":
            for ensembl_biotype in _derived_biotypes[biotype]:
                kwargs["biotype"] = ensembl_biotype
                yield from self.get_features_matching(**kwargs)
            return

        if canonical:
            gene_field_name = "canonical_transcript_id"
            transcript_field_name = "transcript_id"
        else:
            gene_field_name = "gene_id"
            transcript_field_name = "gene_id"

        derived_biotype = kwargs.pop("biotype").lower()
        func = _derived_biotype_funcs[derived_biotype]
        for ensembl_biotype in _derived_biotypes[derived_biotype]:
            kwargs["biotype"] = ensembl_biotype
            for transcript in self._transcript_from_gene(
                gene_field_name=gene_field_name,
                transcript_field_name=transcript_field_name,
                **kwargs,
            ):
                yield func(transcript)

    def get_by_stable_id(self, stable_id: str) -> typing.Iterator[GeneData]:
        yield from self.get_features_matching(stable_id=stable_id, biotype="gene")

    def get_by_symbol(self, symbol: str) -> typing.Iterator[GeneData]:
        yield from self.get_features_matching(symbol=symbol, biotype="gene")

    def get_by_description(self, description: str) -> typing.Iterator[GeneData]:
        yield from self.get_features_matching(description=description, biotype="gene")

    @functools.singledispatchmethod
    def get_feature_children(
        self,
        feature: FeatureDataBase,
    ) -> typing.Iterator[FeatureDataBase]:
        msg = f"{type(feature)=} not supported"
        raise NotImplementedError(msg)

    @get_feature_children.register(GeneData)
    def _(self, gene: GeneData) -> typing.Iterator[TranscriptData]:
        columns = (
            "transcript_id",
            "seqid",
            "coord_system_name",
            "start",
            "stop",
            "strand",
            "transcript_spans",
            "transcript_stable_id",
            "gene_id",
        )
        sql = f"SELECT {','.join(columns)} FROM transcript_attr WHERE transcript_id = {gene['canonical_transcript_id']}"
        for record in self.conn.sql(sql).fetchall():
            transcript = dict(zip(columns, record, strict=True))
            spans = transcript.pop("transcript_spans")
            spans = eti_storage.blob_to_array(spans)
            stable_id = transcript.pop("transcript_stable_id")
            yield TranscriptData(
                **{
                    **transcript,
                    "spans": spans,
                    "biotype": gene.biotype,
                    "stable_id": stable_id,
                    "gene_stable_id": gene.stable_id,
                    "symbol": gene.symbol,
                },
            )

    @get_feature_children.register
    def _(self, transcript: TranscriptData) -> typing.Iterator[CdsData]:
        columns = (
            "transcript_id",
            "seqid",
            "coord_system_name",
            "start",
            "stop",
            "strand",
            "cds_spans",
            "cds_stable_id",
        )
        sql = f"SELECT {','.join(columns)} FROM transcript_attr WHERE transcript_id = {transcript['transcript_id']}"
        if not (record := self.conn.sql(sql).fetchone()):
            msg = f"No CDS spans found for {transcript=}"
            raise ValueError(msg)

        data = dict(zip(columns, record, strict=True))
        if not (spans := data.pop("cds_spans", None)):
            msg = f"No CDS spans found for {transcript=}"
            raise ValueError(msg)

        spans = eti_storage.blob_to_array(spans)
        stable_id = data.pop("cds_stable_id")
        yield CdsData(
            **{
                **data,
                "spans": spans,
                "stable_id": stable_id,
                "gene_stable_id": transcript.gene_stable_id,
            },
        )

    @functools.singledispatchmethod
    def get_feature_parent(self, feature: FeatureDataType) -> FeatureDataBase:
        msg = f"type {feature=} has no parents"
        raise ValueError(msg)

    @get_feature_parent.register
    def _(self, transcript: TranscriptData) -> GeneData:
        sql = f"SELECT {','.join(GENE_ATTR_COLUMNS)} FROM gene_attr WHERE gene_id = {transcript['gene_id']}"
        if not (record := self.conn.sql(sql).fetchone()):
            msg = f"No gene spans found for {transcript=}"
            raise ValueError(msg)

        gene = dict(zip(GENE_ATTR_COLUMNS, record, strict=True))
        spans = numpy.array([sorted([gene["start"], gene["stop"]])], dtype=numpy.int32)
        stable_id = gene.pop("stable_id")
        return GeneData(**{**gene, "spans": spans, "stable_id": stable_id})

    @get_feature_parent.register
    def _(self, cds: CdsData) -> TranscriptData:
        columns = (
            "transcript_id",
            "seqid",
            "coord_system_name",
            "start",
            "stop",
            "strand",
            "transcript_spans",
            "transcript_stable_id",
        )
        sql = f"SELECT {','.join(columns)} FROM transcript_attr WHERE transcript_id = {cds['transcript_id']}"
        if not (record := self.conn.sql(sql).fetchone()):
            msg = f"No transcript spans found for {cds=}"
            raise ValueError(msg)

        transcript = dict(zip(columns, record, strict=True))
        spans = transcript.pop("transcript_spans")
        spans = eti_storage.blob_to_array(spans)
        stable_id = transcript.pop("transcript_stable_id")
        return TranscriptData(
            **{
                **transcript,
                "spans": spans,
                "biotype": cds["biotype"],
                "stable_id": stable_id,
                "gene_stable_id": cds.gene_stable_id,
            },
        )

    def count_distinct(
        self,
        *,
        seqid: StrOrBool = False,
        biotype: OptBool = False,
    ) -> "Table | None":
        if not any((seqid, biotype)):
            return None

        local_vars = locals()
        if constraints := {k: v for k, v in local_vars.items() if isinstance(v, str)}:
            where_clause = f"WHERE {_matching_conditions(equals_conds=constraints)}"
        else:
            where_clause = ""

        header = [c for c in ("biotype", "seqid") if local_vars[c]]
        sql = (
            f"SELECT {', '.join(header)}, COUNT(*) as count FROM gene_attr"
            f" {where_clause} GROUP BY {', '.join(header)};"
        )
        return cogent3.make_table(
            header=[*header, "count"],
            data=self.conn.sql(sql).fetchall(),
            column_templates={"count": lambda x: f"{x:,}"},
        )

    def get_ids_for_biotype(
        self,
        biotype: str,
        seqid: str | None = None,
        limit: int | None = None,
    ) -> list[str]:
        sql = "SELECT stable_id from gene_attr WHERE biotype=?"
        val = (biotype,)
        if seqid:
            sql += " AND seqid=?"
            val = (*val, seqid)
        if limit:
            sql += " LIMIT ?"
            val = (*val, limit)
        return [r[0] for r in self.conn.sql(sql, params=val).fetchall()]

    @functools.cached_property
    def gene_table(self) -> "Table":
        """return a Table with all gene data"""
        columns = (
            "species",
            "name",
            "seqid",
            "source",
            "biotype",
            "start",
            "stop",
            "strand",
            "symbol",
            "description",
        )
        rows = []
        rows.extend(
            [self.species] + [record.get(c, None) for c in columns[1:]]
            for record in self.get_features_matching()
        )
        header = ["stableid" if c == "name" else c for c in columns]
        table = cogent3.make_table(header=header, data=rows)
        # get the numbers of transcripts per gene
        sql = """SELECT ga.stable_id, COUNT(DISTINCT ta.transcript_id) AS distinct_transcript_count
                 FROM transcript_attr ta
                 JOIN gene_attr ga ON ta.gene_id = ga.gene_id
                 GROUP BY ga.stable_id
                 """
        transcript_counts = dict(self.conn.sql(sql).fetchall())
        table = table.with_new_column(
            "num_transcripts",
            lambda x: transcript_counts.get(x, 0),
            columns="stableid",
        )
        # get the transcript biotypes
        sql = """SELECT ga.stable_id, STRING_AGG(DISTINCT ta.transcript_biotype, ',') AS transcript_biotypes
                 FROM transcript_attr ta
                 JOIN gene_attr ga ON ta.gene_id = ga.gene_id
                 GROUP BY ga.stable_id
                 """
        transcript_biotypes = dict(self.conn.sql(sql).fetchall())
        table = table.with_new_column(
            "transcript_biotypes",
            lambda x: transcript_biotypes.get(x, ""),
            columns="stableid",
        )
        columns = (
            "species",
            "seqid",
            "seqid",
            "source",
            "biotype",
            "transcript_biotypes",
            "num_transcripts",
            "start",
            "stop",
            "strand",
            "symbol",
            "description",
        )
        return table.get_columns(columns=columns)


@dataclasses.dataclass
class RepeatView(eti_storage.DuckdbParquetBase, eti_storage.ViewMixin):
    _tables: tuple[str, ...] = tuple(
        core_tables.collect_table_names(*core_tables.REPEAT_ATTR),
    )

    @property
    def conn(self) -> duckdb.DuckDBPyConnection:
        if self._conn is None:
            # trigger the creation of the view using property on super
            super().conn  # noqa: B018
            sql = """CREATE VIEW IF NOT EXISTS repeat_view AS
                    SELECT 
                        rc.repeat_type AS repeat_type,
                        sr.name AS seqid,
                        cs.name AS coord_system_name,
                        rc.repeat_class AS repeat_class,
                        rc.repeat_name AS repeat_name,
                        rf.seq_region_start AS start,
                        rf.seq_region_end AS stop,
                        rf.seq_region_strand AS strand
                    FROM repeat_consensus rc
                    JOIN repeat_feature rf ON rc.repeat_consensus_id = rf.repeat_consensus_id
                    JOIN seq_region sr ON rf.seq_region_id = sr.seq_region_id
                    JOIN coord_system cs ON sr.coord_system_id = cs.coord_system_id
                    """
            self._conn.sql(sql)
        return self._conn

    def num_records(self) -> int:
        """returns the number of rows in repeat_feature"""
        sql = "SELECT COUNT(*) FROM repeat_feature"
        return self.conn.sql(sql).fetchone()[0]

    def get_features_matching(
        self,
        *,
        seqid: OptStr = None,
        biotype: OptStr = None,
        name: OptStr = None,
        start: OptInt = None,
        stop: OptInt = None,
        strand: OptStr = None,
        repeat_type: OptStr = None,
        repeat_class: OptStr = None,
        **kwargs,  # noqa: ANN003
    ) -> typing.Iterator[RepeatData]:
        limit = kwargs.pop("limit", None)
        repeat_type = repeat_type or biotype
        biotype = "repeat"
        repeat_class = repeat_class or name
        name = repeat_class
        local_vars = locals()
        local_vars = {
            k: v
            for k, v in local_vars.items()
            if k not in ("self", "kwargs", "limit", "local_vars", "name", "biotype")
            and v is not None
        }
        core_cols = "seqid", "start", "stop", "strand", "coord_system_name"
        repeat_cols = "repeat_type", "repeat_class", "repeat_name"
        if kwargs := {k: v for k, v in local_vars.items() if v is not None}:
            like_conds = {k: v for k, v in kwargs.items() if k in repeat_cols}
            equals_conds = {k: v for k, v in kwargs.items() if k not in repeat_cols}
        else:
            like_conds = None
            equals_conds = None

        sql = _select_records_sql(
            table_name="repeat_view",
            equals_conds=equals_conds,
            like_conds=like_conds,
            columns=core_cols + repeat_cols,
        )
        sql += f" LIMIT {limit}" if limit else ""
        columns = core_cols + repeat_cols
        for record in self.conn.sql(sql).fetchall():
            data = dict(zip(columns, record, strict=True))
            spans = numpy.array(
                [(data.pop("start"), data.pop("stop"))],
                dtype=numpy.int32,
            )
            data["spans"] = spans
            data["biotype"] = "repeat"
            data["name"] = data["repeat_name"]
            data["start"] = spans.min()
            data["stop"] = spans.max()
            yield RepeatData(**data)

    def get_children_matching(self, **kwargs):
        return ()

    def count_distinct(
        self,
        seqid: StrOrBool = False,
        repeat_type: OptBool = False,
        repeat_class: OptBool = False,
    ) -> "Table | None":
        if not any((seqid, repeat_type, repeat_class)):
            return None

        local_vars = locals()
        if constraints := {k: v for k, v in local_vars.items() if isinstance(v, str)}:
            where_clause = f"WHERE {_matching_conditions(equals_conds=constraints)}"
        else:
            where_clause = ""

        header = [c for c in ("seqid", "repeat_type", "repeat_class") if local_vars[c]]
        sql = (
            f"SELECT {', '.join(header)}, COUNT(*) as count FROM repeat_view"
            f" {where_clause} GROUP BY {', '.join(header)};"
        )
        return cogent3.make_table(
            header=[*header, "count"],
            data=self.conn.sql(sql).fetchall(),
            column_templates={"count": lambda x: f"{x:,}"},
        )


@dataclasses.dataclass
class Annotations(AnnotationDbABC, eti_storage.ViewMixin):
    """virtual genome annotation database that provides access to gene and repeat features"""

    source: dataclasses.InitVar[pathlib.Path | str]
    _source: pathlib.Path = dataclasses.field(init=False)
    biotypes: BiotypeView | None = dataclasses.field(init=True, default=None)
    genes: GeneView | None = dataclasses.field(init=True, default=None)
    repeats: RepeatView | None = dataclasses.field(init=True, default=None)

    def __post_init__(
        self,
        source: pathlib.Path,
    ) -> None:
        source = pathlib.Path(source)
        self._source = source
        if source.is_dir():
            self.biotypes = BiotypeView(source=source)
            self.genes = GeneView(source=source)
            self.repeats = RepeatView(source=source)

    @property
    def source(self) -> pathlib.Path:
        return self._source

    def get_features_matching(
        self,
        *,
        biotype: str,
        **kwargs,
    ) -> typing.Iterator[FeatureDataBase]:
        biotype = biotype or "protein_coding"
        gene_biotypes = set(self.biotypes.distinct) if self.biotypes else set()
        kwargs["biotype"] = biotype
        if biotype in gene_biotypes or is_derived_biotype(biotype):
            view = self.genes
        else:
            kwargs.pop("canonical", None)
            view = self.repeats
        if not view:
            return

        yield from view.get_features_matching(**kwargs)

    def __len__(self) -> int:
        return self.num_records()

    def get_feature_children(self, **kwargs):
        raise NotImplementedError

    def get_feature_parent(self, **kwargs):
        raise NotImplementedError

    def num_matches(self, **kwargs):
        raise NotImplementedError

    def get_ids_for_biotype(self, biotype: str, limit: int | None = None) -> list[str]:
        return self.genes.get_ids_for_biotype(biotype=biotype, limit=limit)

    def count_distinct(self, **kwargs) -> "Table | None":
        return None if self.genes is None else self.genes.count_distinct(**kwargs)

    def num_records(self) -> int:
        """returns the total number of genes and repeat features"""
        num_genes = 0 if self.genes is None else self.genes.num_records()
        num_repeats = 0 if self.repeats is None else self.repeats.num_records()
        return num_genes + num_repeats

    def close(self) -> None:
        if self.biotypes:
            self.biotypes.close()
        if self.genes:
            self.genes.close()
        if self.repeats:
            self.repeats.close()

    def get_ids_for_biotype(
        self,
        *,
        biotype: str,
        seqid: str | list[str] | None = None,
        limit: int | None = None,
    ) -> typing.Iterable[str]:
        if self.genes is None:
            msg = f"no gene data for {self.species}"
            raise ValueError(msg)
        seqids = [seqid] if isinstance(seqid, str | type(None)) else seqid
        for seqid in seqids:
            yield from self.genes.get_ids_for_biotype(
                biotype=biotype,
                seqid=seqid,
                limit=limit,
            )

    def get_records_matching(self, *, seqid: str, **kwargs):
        yield from self.get_features_matching(seqid=seqid, **kwargs)

    def compatible(
        self, other_db: SqliteAnnotationDbMixin, symmetric: bool = True
    ) -> bool:
        return super().compatible(other_db, symmetric)


@dataclasses.dataclass(frozen=True)
class species_seqid:
    species: str
    seqid: str


@functools.cache
def get_species_seqid(*, species: str, seqid: str) -> species_seqid:
    return species_seqid(species, seqid)


@dataclasses.dataclass
class MultispeciesAnnotations(AnnotationDbABC):
    name_map: dict[str, species_seqid]
    species_annotations: dict[str, Annotations]

    def __len__(self) -> int:
        return sum(len(ann) for ann in self.species_annotations.values())

    def get_features_matching(self, *, seqid: str | None = None, **kwargs):
        if seqid not in self.name_map:
            return ()
        sp_sid = self.name_map[seqid]
        db = self.species_annotations[sp_sid.species]
        return db.get_features_matching(seqid=sp_sid.seqid, **kwargs)

    def get_feature_children(self, seqid: str, **kwargs):
        if seqid not in self.name_map:
            return ()
        sp_sid = self.name_map[seqid]
        db = self.species_annotations[sp_sid.species]
        return db.get_feature_children(seqid=sp_sid.seqid, **kwargs)

    def get_feature_parent(self, seqid: str, **kwargs):
        if seqid not in self.name_map:
            return ()
        sp_sid = self.name_map[seqid]
        db = self.species_annotations[sp_sid.species]
        return db.get_feature_parent(seqid=sp_sid.seqid, **kwargs)

    def num_matches(self, seqid: str, **kwargs) -> int:
        """number of records matching arguments in the specified seqid"""
        if seqid not in self.name_map:
            return 0
        sp_sid = self.name_map[seqid]
        db = self.species_annotations[sp_sid.species]
        return db.num_matches(seqid=sp_sid.seqid, **kwargs)

    def get_records_matching(self, *, seqid: str, **kwargs):
        if seqid not in self.name_map:
            return ()
        sp_sid = self.name_map[seqid]
        db = self.species_annotations[sp_sid.species]
        return db.get_features_matching(seqid=sp_sid.seqid, **kwargs)

    def compatible(
        self, other_db: SqliteAnnotationDbMixin, symmetric: bool = True
    ) -> bool:
        return super().compatible(other_db, symmetric)

    def close(self) -> None:
        for db in self.species_annotations.values():
            db.close()
