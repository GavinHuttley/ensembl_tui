import pathlib
import re
import typing

from cogent3 import load_table, make_table
from cogent3.core.tree import PhyloNode

from ensembl_tui import _util as eti_util

if typing.TYPE_CHECKING:  # pragma: no cover
    from cogent3.core.table import Table

SPECIES_NAME = "species.tsv"
StrOrNone = str | None

# essential table columns
TABLE_COLUMNS = ["abbrev", "common_name", "genome_name", "db_prefix"]

# some species name records have an accession identifier, which contains digits
_accession = re.compile(r"\d")


def load_species(species_path: eti_util.PathType | None) -> "Table":
    """returns Table from species_path

    if species_path does not exist, defaults to default one"""
    if species_path is None or not pathlib.Path(species_path).exists():
        species_path = eti_util.get_resource_path(SPECIES_NAME)

    return load_table(species_path)


def make_species_map(species_path: eti_util.PathType | None) -> "SpeciesNameMap":
    """returns a SpeciesNameMap from the default species table"""
    species_table = load_species(species_path)
    return SpeciesNameMap.from_table(species_table)


def _genome_name_to_species_name(genome_name: str) -> str:
    """convert an ensembl db prefix to a species name"""
    name = genome_name.replace("_", " ").capitalize().split()
    return " ".join(p for p in name if not _accession.search(p))


def _latin_to_abbrev(genome_to_abrv: dict[str, str], latin_name: str) -> str | None:
    """convert a latin name to a abbreviation"""
    latin_name = latin_name.replace(" ", "_")
    found = []
    for genome_name, abrv in genome_to_abrv.items():
        if genome_name == latin_name:
            return abrv
        if genome_name.startswith(latin_name):
            found.append(abrv)
    if len(found) == 1:
        return found[0]

    return None


class SpeciesNameMap:
    """mapping between common names and latin names"""

    def __init__(
        self,
        *,
        abbrev_common: dict[str, str],
        abbrev_genome: dict[str, str],
        abbrev_db: dict[str, str],
    ) -> None:
        """provides mappings between abbrev:genome:db_prefix:common name"""
        # all queries will go via abbreviations
        # delete below
        self._db_prefix_to_common = {
            db: abbrev_common[abrv] for abrv, db in abbrev_db.items()
        }
        # common name mappings
        self._abrv_to_common = abbrev_common
        self._common_to_abrv = {c: a for a, c in abbrev_common.items()}
        # db prefix mappings
        self._abrv_to_db = abbrev_db
        self._db_to_abrv = {db: abrv for abrv, db in abbrev_db.items()}
        # species name mappings
        self._abrv_to_genome = abbrev_genome
        self._genome_to_abrv = {genome: abrv for abrv, genome in abbrev_genome.items()}

    def __str__(self) -> str:
        return str(self.to_table())

    def __repr__(self) -> str:
        return repr(self.to_table())

    def __contains__(self, item: str) -> bool:
        return bool(self._get_abbrev_for_name(item))

    def _repr_html_(self) -> str:
        table = self.to_table()
        return table._repr_html_()

    def _handle_optional_errors(self, name: str, level: str) -> None:
        msg = f"Unknown {name!r}"
        if level == "raise":
            raise ValueError(msg)
        if level == "warn":
            print(f"WARN: {msg}")

    def _get_abbrev_for_name(self, name: str) -> StrOrNone:
        if not isinstance(name, (str, bytes)):
            return False
        name = name.lower()
        if name in self._abrv_to_db:
            return name

        for attr in (self._common_to_abrv, self._db_to_abrv, self._genome_to_abrv):
            if name in attr:
                return attr[name]

        if " " in name and (abrv := _latin_to_abbrev(self._genome_to_abrv, name)):
            return abrv

        return None

    def get_common_name(
        self, name: str, level: typing.Literal["ignore", "warn", "raise"] = "raise"
    ) -> StrOrNone:
        """returns the common name"""
        if abrv := self._get_abbrev_for_name(name):
            return self._abrv_to_common[abrv]

        self._handle_optional_errors(name, level)
        return None

    def get_species_name(
        self, name: str, level: typing.Literal["ignore", "warn", "raise"] = "ignore"
    ) -> StrOrNone:
        """returns the latin name"""
        if abrv := self._get_abbrev_for_name(name):
            return _genome_name_to_species_name(self._abrv_to_genome[abrv])

        self._handle_optional_errors(name, level)
        return None

    def get_genome_name(
        self, name: str, level: typing.Literal["ignore", "warn", "raise"] = "ignore"
    ) -> str | None:
        """returns the Ensembl genome name"""
        if abrv := self._get_abbrev_for_name(name):
            return self._abrv_to_genome[abrv]

        self._handle_optional_errors(name, level)
        return None

    def get_ensembl_db_prefix(
        self, name: str, level: typing.Literal["ignore", "warn", "raise"] = "ignore"
    ) -> str | None:
        """returns the Ensembl db prefix"""
        if abrv := self._get_abbrev_for_name(name):
            return self._abrv_to_db[abrv]

        self._handle_optional_errors(name, level)
        return None

    def get_abbreviation(
        self, name: str, level: typing.Literal["ignore", "warn", "raise"] = "ignore"
    ) -> StrOrNone:
        """returns the abbreviation for the given name"""
        if abrv := self._get_abbrev_for_name(name):
            return abrv

        self._handle_optional_errors(name, level)
        return None

    def to_table(self) -> "Table":
        """returns cogent3 Table"""
        rows = [
            [abrv, self._abrv_to_common[abrv], self._abrv_to_genome[abrv], db_prefix]
            for abrv, db_prefix in self._abrv_to_db.items()
        ]
        return make_table(
            header=TABLE_COLUMNS,
            data=rows,
            space=2,
        ).sorted()

    @classmethod
    def from_table(cls, species_table: "Table") -> "SpeciesNameMap":
        """uses TABLE_COLUMNS from species_table to create a SpeciesNameMap"""
        from ensembl_tui._name import EnsemblDbName

        abrv_genome = {}
        abrv_common = {}
        abrv_db = {}
        for abrv, common, genome, db_prefix in species_table.to_list(
            columns=TABLE_COLUMNS
        ):
            abrv = abrv.strip().lower()
            abrv_genome[abrv] = genome.lower()
            abrv_common[abrv] = common.lower()
            abrv_db[abrv] = EnsemblDbName(db_prefix.lower()).prefix

        return cls(
            abbrev_common=abrv_common, abbrev_genome=abrv_genome, abbrev_db=abrv_db
        )

    def get_subset(self, names: list[str]) -> "SpeciesNameMap":
        """returns a species map subset for the given names"""
        selected = {self.get_genome_name(n, level="raise") for n in names}
        data = self.for_storage()
        subset = {k: v for k, v in data.items() if k == "header" or k in selected}
        return self.from_storage(subset)

    def for_storage(self) -> dict[str, str]:
        """creates a dict suitable for cfg storage"""
        delim = "\t"
        table = self.to_table()
        primary = "genome_name"
        header = [primary, *[c for c in table.header if c != primary]]
        result = {"header": delim.join(header)}
        for row in table.to_list(columns=header):
            result[row[0]] = delim.join(row[1:])
        return result

    @classmethod
    def from_storage(cls, d: dict[str, str]) -> "SpeciesNameMap":
        """creates a SpeciesNameMap from a dict created by for_storage"""
        delim = "\t"
        header = d.pop("header").split(delim)
        rows = [[key, *row.split(delim)] for key, row in d.items()]
        table = make_table(header=header, data=rows, space=2)
        return cls.from_table(table)


def species_from_ensembl_tree(
    tree: PhyloNode, species_map: SpeciesNameMap
) -> dict[str, str]:
    """get species identifiers from an Ensembl tree"""
    tip_names = tree.get_tip_names()
    selected_species = {}
    for tip_name in tip_names:
        name_fields = tip_name.lower().split("_")
        # produce parts of name starting with highly specific to
        # more general and look for matches
        for j in range(len(name_fields) + 1, 1, -1):
            n = "_".join(name_fields[:j])
            if n in species_map:
                selected_species[species_map.get_genome_name(n)] = tip_name
                break
        else:
            msg = f"cannot establish species for {'_'.join(name_fields)}"
            raise ValueError(msg)

    return selected_species


def make_unique_abbrevs(names: list[str]) -> dict[str, str]:
    """makes unique abbreviations from the species names"""
    abbrevs: dict[str, str] = {}
    for name in names:
        parts = [p for p in name.split("_") if not _accession.search(p)]
        parts = [parts[0][:3], *[p[:4] for p in parts[1:]], ""]
        i = 1
        while "-".join(p for p in parts if p) in abbrevs.values():
            i += 1
            parts[-1] = f"{i}"

        abbrevs[name] = "-".join(p for p in parts if p)

    return abbrevs
