import dataclasses
import pathlib
import sys
import typing
from collections.abc import Callable

import cogent3
import cogent3_h5seqs as c3h5
import numpy
from cogent3.app.composable import define_app
from cogent3.core import alphabet as c3alpha
from cogent3.core.sequence import Sequence
from cogent3.core.table import Table
from cogent3.parse.fasta import iter_fasta_records

import ensembl_tui._annotation as eti_annots
from ensembl_tui import _config as eti_config
from ensembl_tui import _species as eti_species
from ensembl_tui import _util as eti_util

SEQ_STORE_NAME = f"genome-seqs.{c3h5.UNALIGNED_SUFFIX}"

DNA = cogent3.get_moltype("dna")
alphabet = DNA.most_degen_alphabet()  # type: ignore  # noqa: PGH003
bytes_to_array = c3alpha.bytes_to_array(
    chars=alphabet.as_bytes(),
    dtype=numpy.uint8,
    delete=b" \n\r\t",
)


def _rename(label: str) -> str:
    return label.split()[0]


@define_app
class fasta_to_hdf5:  # noqa: N801
    def __init__(
        self,
        config: eti_config.Config,
        label_to_name: Callable[[str], str] = _rename,
    ) -> None:
        self.config = config
        self.label_to_name = label_to_name

    def main(self, db_name: str) -> bool:
        src_dir = self.config.staging_genomes / db_name
        dest_dir = self.config.install_genomes / db_name
        dest_dir.mkdir(parents=True, exist_ok=True)
        out_path = dest_dir / SEQ_STORE_NAME
        # we directly use the cogent3_h5seqs library to create the
        # unaligned hdf5 file. This library is valid storage for
        # cogent3 collections
        alpha = cogent3.get_moltype("dna").most_degen_alphabet()
        seq_store = c3h5.make_unaligned(out_path, alphabet=alpha, mode="w")
        seq_store.set_attr(
            "species",
            # we have to coerce the species name from a case insensitive string
            # to a standard python string
            str(self.config.species_map.get_genome_name(db_name)),
            force=True,
        )
        src_dir = src_dir / "fasta"
        for path in src_dir.glob("*.fa.gz"):
            seqs = dict(
                iter_fasta_records(
                    path,
                    converter=bytes_to_array,
                    label_to_name=self.label_to_name,
                )
            )
            seq_store.add_seqs(seqs, force_unique_keys=False)
            del seqs

        seq_store.close()

        return True


@dataclasses.dataclass(slots=True)
class genome_segment:  # noqa: N801
    species: str
    seqid: str
    start: int
    stop: int
    strand: str
    unique_id: str | None = None

    def __post_init__(self) -> None:
        self.unique_id = (
            eti_util.sanitise_stableid(self.unique_id)
            if self.unique_id
            else f"{self.species}-{self.seqid}-{self.start}-{self.stop}"
        )

    @property
    def source(self) -> str | None:
        return self.unique_id


def load_genome(*, config: eti_config.InstalledConfig, species: str):
    """returns the genome with annotations"""
    genome_path = config.installed_genome(species) / SEQ_STORE_NAME
    storage = c3h5.load_seqs_data_unaligned(genome_path)
    dna = cogent3.get_moltype("dna")
    ann = eti_annots.Annotations(source=config.installed_genome(species))
    return cogent3.make_unaligned_seqs(
        storage,
        moltype=dna,
        annotation_db=ann,
        info={"species": species},
    )


def get_seqs_for_ids(
    *,
    config: eti_config.InstalledConfig,
    species: str,
    names: list[str],
    make_seq_name: typing.Callable | None = None,
) -> typing.Iterable[Sequence]:
    genome = load_genome(config=config, species=species)
    # is it possible to do batch query for all names?
    ann_db = genome.annotation_db

    for name in names:
        cds = list(
            ann_db.get_features_matching(
                name=name,
                biotype="cds",
                canonical=True,
            ),
        )
        if not cds:
            continue

        feature = genome.make_feature(feature=cds[0])
        seq = feature.get_slice()
        if callable(make_seq_name):
            seq.name = make_seq_name(feature)
        else:
            seq.name = f"{species}-{name}"
        seq.info["species"] = species
        seq.info["name"] = name
        # disconnect from annotation so the closure of the genome
        # does not cause issues when run in parallel
        seq.annotation_db = None
        yield seq

    del genome


def load_annotations_for_species(*, path: pathlib.Path) -> eti_annots.Annotations:
    """returns the annotation Db for species"""
    if not path.exists():
        eti_util.print_colour(
            text=f"{path.name!r} is missing",
            colour="red",
        )
        sys.exit(1)
    try:
        return eti_annots.Annotations(source=path)
    except FileNotFoundError:
        eti_util.print_colour(
            text=f"expected files not in {str(path)!r}",
            colour="red",
        )
        sys.exit(1)


def _get_all_gene_segments(
    *,
    annot_db: eti_annots.Annotations,
    limit: int | None,
    biotype: str | None,
) -> list[eti_annots.GeneData]:
    return list(annot_db.get_features_matching(biotype=biotype, limit=limit))


def _get_selected_gene_segments(
    *,
    annot_db: eti_annots.Annotations,
    limit: int | None,
    stableids: list[str],
    biotype: str | None,
) -> list[eti_annots.GeneData]:
    result = []
    for stable_id in stableids:
        record = list(
            annot_db.get_features_matching(
                biotype=biotype,
                stable_id=stable_id,
                limit=limit,
            ),
        )
        result.extend(record)
    return result


def get_gene_segments(
    *,
    annot_db: eti_annots.Annotations,
    limit: int | None = None,
    species: str | None = None,
    stableids: list[str] | None = None,
    biotype: str = "protein_coding",
) -> list[genome_segment]:
    """return genome segment information for genes

    Parameters
    ----------
    annot_db
        feature db
    limit
        limit number of records to
    species
        species name, overrides inference from annot_db.source
    """
    species = species or annot_db.source.parent.name
    records = (
        _get_selected_gene_segments(
            annot_db=annot_db,
            limit=limit,
            stableids=stableids,
            biotype=biotype,
        )
        if stableids
        else _get_all_gene_segments(annot_db=annot_db, limit=limit, biotype=biotype)
    )
    for i, record in enumerate(records):
        segment = genome_segment(
            species=species,
            start=record["start"],
            stop=record["stop"],
            strand=record["strand"],
            seqid=record["seqid"],
            unique_id=record["name"],
        )
        records[i] = segment
    return records


def get_gene_table_for_species(
    *,
    annot_db: eti_annots.Annotations,
    limit: int | None = None,
) -> Table:
    """
    returns gene data from Annotations

    Parameters
    ----------
    annot_db
        feature db
    limit
        limit number of records to
    """
    table = annot_db.genes.gene_table
    if limit:
        table = table[:limit]
    return table


def get_species_gene_summary(
    *,
    annot_db: eti_annots.Annotations,
    species_map: eti_species.SpeciesNameMap,
    species: str | None = None,
) -> Table:
    """
    returns the Table summarising data for species_name

    Parameters
    ----------
    annot_db
        feature db
    species_map
        map of common, latin and ensembl db names
    species
        species name, overrides inference from annot_db.source
    """
    # for now, just biotype
    species = species or annot_db.source.parent.name
    counts = annot_db.biotypes.count_distinct()
    try:
        common_name = species_map.get_species_name(species)
    except ValueError:
        common_name = species

    counts.title = f"{common_name} features"
    counts.format_column("count", lambda x: f"{x:,}")
    return counts


def get_species_repeat_summary(
    *,
    annot_db: eti_annots.Annotations,
    species_map: eti_species.SpeciesNameMap,
    species: str | None = None,
) -> Table:
    """
    returns the Table summarising repeat data for species_name

    Parameters
    ----------
    annot_db
        feature db
    species_map
        map of common, latin and ensembl db names
    species
        species name, overrides inference from annot_db.source
    """
    # for now, just biotype
    species = species or annot_db.source.parent.name
    counts = annot_db.repeats.count_distinct(repeat_class=True, repeat_type=True)
    try:
        common_name = species_map.get_species_name(species)
    except ValueError:
        common_name = species

    counts = counts.sorted(columns=["repeat_type", "count"])
    counts.title = f"{common_name} repeat"
    counts.format_column("count", lambda x: f"{x:,}")
    return counts
