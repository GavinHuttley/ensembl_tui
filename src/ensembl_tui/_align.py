import dataclasses
import sys
import typing
from collections import defaultdict

import cogent3
import numpy
from cogent3.core import alignment as c3_align
from cogent3.core.location import _DEFAULT_GAP_DTYPE, IndelMap
from scinexus.composable import define_app

from ensembl_tui import _annotation as eti_ann
from ensembl_tui import _config as eti_config
from ensembl_tui import _genome as eti_genome
from ensembl_tui import _storage_mixin as eti_storage
from ensembl_tui import _util as eti_util

DNA = cogent3.get_moltype("dna")

_no_gaps = numpy.array([], dtype=_DEFAULT_GAP_DTYPE)

ALIGN_STORE_SUFFIX = "parquet"

ALIGN_ATTR_SCHEMA = (
    "align_id INTEGER PRIMARY KEY DEFAULT nextval('align_id_seq')",
    "source TEXT",
    "block_id BIGINT",
    "species TEXT",
    "seqid TEXT",
    "start INTEGER",
    "stop INTEGER",
    "strand TINYINT",
    "gap_spans BLOB",
)
ALIGN_ATTR_COLS = eti_util.make_column_constant(ALIGN_ATTR_SCHEMA)

VT = str | int | numpy.ndarray


@dataclasses.dataclass(slots=True)
class AlignRecord:
    """a record from an AlignDb

    Notes
    -----
    Can return fields as attributes or like a dict using the field name as
    a string.
    """

    source: str
    block_id: int
    species: str
    seqid: str
    start: int
    stop: int
    strand: int
    gap_spans: numpy.ndarray

    def __post_init__(self) -> None:
        if isinstance(self.strand, str):
            self.strand = -1 if self.strand.startswith("-") else 1

    def __getitem__(self, item: str) -> VT:
        return getattr(self, item)

    def __setitem__(self, item: str, value: VT) -> None:
        setattr(self, item, value)

    @property
    def identity(self) -> tuple[int, str, str, int, int, int]:
        """the fields that identify this record within an alignment store"""
        return (
            self.block_id,
            self.species,
            self.seqid,
            self.start,
            self.stop,
            self.strand,
        )

    def __eq__(self, other: object) -> bool:
        if not isinstance(other, AlignRecord):
            return False

        # array_equal rather than comparing the arrays directly, which raises
        # when their shapes differ. equal identities with different gaps hash
        # alike, so they do meet in the set get_records_matching builds
        return self.identity == other.identity and numpy.array_equal(
            self.gap_spans,
            other.gap_spans,
        )

    def __hash__(self) -> int:
        return hash(self.identity)

    @property
    def cum_gap_data(self) -> tuple[numpy.ndarray, numpy.ndarray]:
        """gap positions in sequence coordinates and cumulative gap lengths

        Notes
        -----
        The lengths are cumulative, which is how IndelMap stores them, so
        they are passed to it as cum_gap_lengths and not as gap_lengths.
        """
        if len(self.gap_spans):
            gap_pos, cum_gap_lengths = self.gap_spans.T
        else:
            gap_pos, cum_gap_lengths = _no_gaps.copy(), _no_gaps.copy()

        return gap_pos, cum_gap_lengths

    def to_dict(self) -> dict:
        return dataclasses.asdict(self)

    def to_record(self, columns: tuple[str, ...]) -> tuple:
        data = self.to_dict()
        data["gap_spans"] = eti_storage.array_to_blob(self.gap_spans)
        return tuple(data[c] for c in columns)


ReturnType = tuple[str, tuple]  # the sql statement and corresponding values


# TODO add a table and methods to support storing the species tree used
#  for the alignment and for getting the species tree
@dataclasses.dataclass(slots=True)
class AlignDb(eti_storage.DuckdbParquetBase):
    _tables: tuple[str, ...] = ("align_blocks",)

    def _get_block_id(
        self,
        *,
        species: str,
        seqid: str,
        start: int | None,
        stop: int | None,
    ) -> list[tuple[int]]:
        sql = (
            f"SELECT DISTINCT block_id from {self._tables[0]} "
            "WHERE species = ? AND seqid = ?"
        )
        values = species, seqid
        # a block matches if it overlaps the query interval at all. coordinates
        # are half open, so a block must finish strictly after the query start
        # and begin strictly before the query stop. an absent bound means the
        # query is open ended on that side. testing the two bounds separately
        # keeps blocks that lie wholly inside the query, which a test of only
        # the query end points would discard
        if start is not None:
            sql = f"{sql} AND ? < stop"
            values += (start,)
        if stop is not None:
            sql = f"{sql} AND start < ?"
            values += (stop,)

        return self.conn.sql(sql, params=values).fetchall()

    def get_records_matching(
        self,
        *,
        species: str,
        seqid: str,
        start: int | None = None,
        stop: int | None = None,
    ) -> typing.Iterable[set[AlignRecord]]:
        # make sure python, not numpy, integers
        start = None if start is None else int(start)
        stop = None if stop is None else int(stop)

        # We need the block IDs for all records for a species whose coordinates
        # lie in the range (start, stop). We then search for all records with
        # each block id. We return full records.
        # Client code is responsible for creating Aligned sequence instances
        # and the Alignment.

        # TODO: there's an issue here with records being duplicated, solved
        #   for now by making AlignRecord hashable and using a set for block_ids
        block_ids = {
            r[0]
            for r in self._get_block_id(
                species=species,
                seqid=seqid,
                start=start,
                stop=stop,
            )
        }
        if not block_ids:
            return []
        columns = tuple(c for c in ALIGN_ATTR_COLS if c != "align_id")
        col_order = ", ".join(columns)
        values = ", ".join("?" * len(block_ids))
        sql = f"SELECT {col_order} from {self._tables[0]} WHERE block_id IN ({values})"
        results = defaultdict(set)
        for record in self.conn.sql(sql, params=tuple(block_ids)).fetchall():
            data = dict(zip(columns, record, strict=True))
            data["gap_spans"] = eti_storage.blob_to_array(data["gap_spans"])
            results[data["block_id"]].add(AlignRecord(**data))

        return results.values()

    def get_species_names(self) -> list[str]:
        """return the list of species names"""
        return list(self.get_distinct("species"))

    def get_distinct(self, field: str) -> list[str]:
        return [
            r[0]
            for r in self.conn.sql(
                f"SELECT DISTINCT {field} from {self._tables[0]}",
            ).fetchall()
        ]

    def num_records(self) -> int:
        return typing.cast(
            "int",
            self.conn.sql(f"SELECT COUNT(*) from {self._tables[0]}").fetchone()[0],
        )

    def close(self) -> None:
        """closes duckdb storage"""
        self.conn.close()


def _overlapping_ref_records(
    block: set[AlignRecord],
    ref_species: str,
    seqid: str,
    ref_start: int | None,
    ref_end: int | None,
) -> list[AlignRecord]:
    """records for the reference sequence that overlap the query interval

    Notes
    -----
    A single alignment block can contain more than one segment of the reference
    sequence, for example where a region has been duplicated. Only the segments
    overlapping the query are wanted, and they are returned in genomic order so
    the result does not depend on the iteration order of the block set.
    """
    records = [r for r in block if r.species == ref_species and r.seqid == seqid]
    if ref_start is not None:
        records = [r for r in records if ref_start < r.stop]
    if ref_end is not None:
        records = [r for r in records if r.start < ref_end]
    return sorted(records, key=lambda r: (r.start, r.stop))


def get_alignment(
    align_db: AlignDb,
    genomes: dict[str, c3_align.SequenceCollection],
    ref_species: str,
    seqid: str,
    ref_start: int | None = None,
    ref_end: int | None = None,
    namer: typing.Callable[[str, str, int, int], str] | None = None,
    mask_features: list[str] | None = None,
    shadow: bool = False,
    mask_ref: bool = False,
) -> typing.Iterable[c3_align.Alignment]:
    """yields cogent3 new type Alignments"""

    if ref_species not in genomes:
        msg = f"unknown species {ref_species!r}"
        raise ValueError(msg)

    align_records = align_db.get_records_matching(
        species=ref_species,
        seqid=seqid,
        start=ref_start,
        stop=ref_end,
    )
    # sample the sequences
    for block in align_records:
        # we get the gaps corresponding to the reference sequence
        # and convert them to a IndelMap instance. We then convert
        # the ref_start, ref_end into align_start, align_end. Those values are
        # used for all other species -- they are converted into sequence
        # coordinates for each species -- selecting their sequence,
        # building the Aligned instance, and selecting the annotation subset.
        ref_records = _overlapping_ref_records(
            block,
            ref_species,
            seqid,
            ref_start,
            ref_end,
        )
        if not ref_records:
            msg = f"no matching alignment record for {ref_species!r}"
            raise ValueError(msg)

        for ref_record in ref_records:
            yield _make_alignment(
                block=block,
                ref_record=ref_record,
                genomes=genomes,
                ref_species=ref_species,
                ref_start=ref_start,
                ref_end=ref_end,
                namer=namer,
                mask_features=mask_features,
                shadow=shadow,
                mask_ref=mask_ref,
            )


def _make_alignment(
    *,
    block: set[AlignRecord],
    ref_record: AlignRecord,
    genomes: dict[str, c3_align.SequenceCollection],
    ref_species: str,
    ref_start: int | None,
    ref_end: int | None,
    namer: typing.Callable[[str, str, int, int], str] | None,
    mask_features: list[str] | None,
    shadow: bool,
    mask_ref: bool,
) -> c3_align.Alignment:
    """builds the alignment for one segment of the reference sequence"""
    # ref_start, ref_end are genomic positions and the ref_record
    # start / stop are also genomic positions
    genome_start = ref_record.start
    genome_end = ref_record.stop
    gap_pos, cum_gap_lengths = ref_record.cum_gap_data
    imap = IndelMap(
        gap_pos=gap_pos,
        cum_gap_lengths=cum_gap_lengths,
        parent_length=genome_end - genome_start,
    )

    # We use the IndelMap object to identify the alignment
    # positions the ref_start / ref_end correspond to. The alignment
    # positions are used below for slicing each sequence in the
    # alignment.

    # make sure the sequence start and stop are within this
    # aligned block
    seq_start = max(ref_start or genome_start, genome_start)
    seq_end = min(ref_end or genome_end, genome_end)
    # make these coordinates relative to the aligned segment
    if ref_record.strand == -1:
        # if record is on minus strand, then genome stop is
        # the alignment start
        seq_start, seq_end = genome_end - seq_end, genome_end - seq_start
    else:
        seq_start = seq_start - genome_start
        seq_end = seq_end - genome_start

    align_start = imap.get_align_index(seq_start)
    align_end = imap.get_align_index(seq_end)

    seqs = {}
    gaps = {}
    offsets = {}
    reversed_seqs = set()
    seqid_species = {}
    ann_dbs = {}
    for align_record in block:
        record_species = align_record.species
        genome = genomes[record_species]
        # We need to convert the alignment coordinates into sequence
        # coordinates for this species.
        genome_start = align_record.start
        genome_end = align_record.stop
        gap_pos, cum_gap_lengths = align_record.cum_gap_data
        imap = IndelMap(
            gap_pos=gap_pos,
            cum_gap_lengths=cum_gap_lengths,
            parent_length=genome_end - genome_start,
        )

        # We use the alignment indices derived for the reference sequence
        # above
        seq_start = imap.get_seq_index(align_start)
        seq_end = imap.get_seq_index(align_end)
        seq_length = seq_end - seq_start
        if align_record.strand == -1:
            # if it's neg strand, the alignment start is the genome stop
            seq_start = imap.parent_length - seq_end

        start = genome_start + seq_start
        stop = genome_start + seq_start + seq_length
        s = genome.seqs[align_record.seqid][start:stop]

        if namer:
            name = namer(align_record.species, align_record.seqid, start, stop)
        else:
            name = f"{align_record.species}:{align_record.seqid}:{start}-{stop}"

        s.name = name
        s.replace_annotation_db(None)
        # we now trim the gaps for this sequence to the sub-alignment
        imap = imap[align_start:align_end]

        if not namer:
            s.name = f"{s.name}:{align_record.strand}"

        if s.name in seqs:
            eti_util.print_colour(f"duplicated {s.name}", colour="yellow")

        if align_record.strand == -1:
            s = s.rc()
            reversed_seqs.add(s.name)

        seqs[s.name] = numpy.array(s)
        gaps[s.name] = imap.array
        if mask_ref and record_species != ref_species:
            # limit features to only those from the reference genome
            continue

        offsets[s.name] = genome_start + seq_start
        seqid_species[s.name] = eti_ann.get_species_seqid(
            species=record_species,
            seqid=align_record.seqid,
        )
        ann_dbs[record_species] = genome.annotation_db

    aln_data = c3_align.AlignedSeqsData.from_seqs_and_gaps(
        seqs=seqs,
        gaps=gaps,
        alphabet=DNA.most_degen_alphabet(),
        offset=offsets,
        reversed_seqs=reversed_seqs,
    )
    aln = c3_align.Alignment(seqs_data=aln_data, moltype=DNA)

    ann_db = eti_ann.MultispeciesAnnotations(
        name_map=seqid_species,
        species_annotations=ann_dbs,
    )
    aln.annotation_db = ann_db

    if mask_features:
        aln = aln.with_masked_annotations(biotypes=mask_features, shadow=shadow)

    return aln


@define_app
class construct_alignment:  # noqa: N801
    """reassemble an alignment that maps to a given genomic segment

    If the segment spans multiple alignments these are joined using
    the sep character.
    """

    def __init__(
        self,
        align_db: AlignDb,
        genomes: dict[str, c3_align.SequenceCollection],
        mask_features: list[str] | None = None,
        shadow: bool = False,
        mask_ref: bool = False,
        sep: str = "?",
    ) -> None:
        self._align_db = align_db
        self._genomes = genomes
        self._mask_features = mask_features
        self._shadow = shadow
        self._sep = sep
        self._ref_only = mask_ref

    def main(self, segment: eti_genome.genome_segment) -> list[c3_align.Alignment]:
        results = []
        for aln in get_alignment(
            self._align_db,
            self._genomes,
            segment.species,
            segment.seqid,
            segment.start,
            segment.stop,
            mask_features=self._mask_features,
            shadow=self._shadow,
            mask_ref=self._ref_only,
        ):
            aln.source = segment.source
            results.append(aln)

        return results


def load_aligndb(config: eti_config.InstalledConfig, align_name: str) -> AlignDb:
    """returns an AlignDb instance for the given config"""
    align_name = eti_util.strip_quotes(align_name)
    align_path = config.path_to_alignment(align_name, ALIGN_STORE_SUFFIX)
    if align_path is None:
        eti_util.print_colour(
            text=f"{align_name!r} does not match any alignments under '{config.aligns_path}'",
            colour="red",
        )
        available = "\n".join(
            [
                fn.stem
                for fn in config.aligns_path.glob("*")
                if not fn.name.startswith(".") and fn.is_dir()
            ],
        )
        eti_util.print_colour(text=f"Available alignments:\n{available}", colour="red")
        sys.exit(1)

    return AlignDb(source=align_path)
