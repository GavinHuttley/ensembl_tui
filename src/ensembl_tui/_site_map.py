import typing
from dataclasses import dataclass
from functools import cache

from cogent3.util.misc import extend_docstring_from

_ensembl_site_map = {}


class register_ensembl_site_map:
    """
    registration decorator for Ensembl site-map classes

    The registration key must be a string that of the domain name.

    Parameters
    ----------
    domain: str of domain name, must be unique
    """

    def __init__(self, domain: str):
        if not isinstance(domain, str):
            raise TypeError(f"{domain!r} is not a string")

        domain = domain.strip()
        if not domain:
            raise ValueError("cannot have empty string domain")

        assert domain not in _ensembl_site_map, (
            f"{domain!r} already in {list(_ensembl_site_map)}"
        )

        self._domain = domain

    def __call__(self, func):
        # pass through
        _ensembl_site_map[self._domain] = func
        return func


StrOrNone = typing.Union[str, type(None)]


@dataclass(slots=True)
class SiteMap:
    """records the locations of specific attributes relative to an Ensembl release"""

    site: str
    db_host: str
    db_port: int
    remote_path: str
    species_file_name: str
    _seqs_path: str = "fasta"
    _annotations_path: str = "mysql"
    _alignments_path: str | None = None
    _homologies_path: str | None = None
    _trees_path: str | None = None

    def get_seqs_path(self, ensembl_name: str, collection_name: str | None) -> str:
        """path to unmasked genome sequences"""
        # this needs to be self._seqs_path/ensembl_name/dna
        # except when it's part of a collection, in which case its
        # self._seqs_path/collection_name/ensembl_name/dna
        if collection_name:
            return f"{self._seqs_path}/{collection_name}/{ensembl_name}/dna"
        return f"{self._seqs_path}/{ensembl_name}/dna"

    def get_annotations_path(self, ensembl_name: str) -> str:
        return f"{self._annotations_path}/{ensembl_name}"

    @property
    def alignments_path(self) -> StrOrNone:
        return self._alignments_path

    @property
    def homologies_path(self) -> StrOrNone:
        return self._homologies_path

    @property
    def trees_path(self) -> StrOrNone:
        return self._trees_path

    def get_remote_release_path(self, release: str) -> str:
        return f"{self.remote_path}/release-{release}"


@extend_docstring_from(SiteMap)
@register_ensembl_site_map("main")
@register_ensembl_site_map("vertebrates")
@register_ensembl_site_map("ftp.ensembl.org")
def ensembl_main_sitemap() -> SiteMap:
    """the main Ensembl site map"""
    return SiteMap(
        site="ftp.ensembl.org",
        _alignments_path="maf/ensembl-compara/multiple_alignments",
        _homologies_path="tsv/ensembl-compara/homologies",
        _trees_path="compara/species_trees",
        db_host="ensembldb.ensembl.org",
        db_port=3306,
        remote_path="pub",
        species_file_name="species_EnsemblVertebrates.txt",
    )


@register_ensembl_site_map("metazoa")
@register_ensembl_site_map("ftp.ensemblgenomes.org")
def ensembl_metazoa_sitemap() -> SiteMap:
    """the metazoa Ensembl site map"""
    return SiteMap(
        site="ftp.ensemblgenomes.org",
        _alignments_path=None,
        _homologies_path="tsv/ensembl-compara/homologies",
        _trees_path=None,
        db_host="mysql-eg-publicsql.ebi.ac.uk",
        db_port=4157,
        remote_path="pub/metazoa",
        species_file_name="species_EnsemblMetazoa.txt",
    )


@register_ensembl_site_map("protists")
def ensembl_protists_sitemap() -> SiteMap:
    """the protists Ensembl site map"""
    return SiteMap(
        site="ftp.ensemblgenomes.org",
        _alignments_path=None,
        _homologies_path="tsv/ensembl-compara/homologies",
        _trees_path=None,
        db_host="mysql-eg-publicsql.ebi.ac.uk",
        db_port=4157,
        remote_path="pub/protists",
        species_file_name="species_EnsemblProtists.txt",
    )


@cache
def get_site_map(domain: str) -> SiteMap:
    """Returns a site map instance for the specified domain.

    Parameters
    ----------
    domain : str
        The Ensembl domain name (e.g., 'main', 'vertebrates', 'metazoa', 'protists').
        Use get_site_map_names() to see all available options.

    Returns
    -------
    SiteMap
        Site configuration containing FTP host, database connection info, and paths.

    Raises
    ------
    KeyError
        If the domain is not registered.
    """
    return _ensembl_site_map[domain]()


def get_site_map_names() -> list[str]:
    """Returns all registered Ensembl domain names."""
    return [n for n in _ensembl_site_map if not n.startswith("ftp")]
