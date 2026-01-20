import configparser
import fnmatch
import pathlib
import sys
import typing
from collections.abc import Generator, Sequence
from dataclasses import dataclass

from ensembl_tui import _site_map as eti_site_map
from ensembl_tui import _species as eti_species
from ensembl_tui import _util as eti_util

if typing.TYPE_CHECKING:  # pragma: no cover
    from cogent3.core.table import Table

INSTALLED_CONFIG_NAME = "installed.cfg"
DOWNLOADED_CONFIG_NAME = "downloaded.cfg"

_COMPARA_NAME: str = "compara"
_ALIGNS_NAME: str = "aligns"
_HOMOLOGIES_NAME: str = "homologies"
_GENOMES_NAME: str = "genomes"

_VERSION_SECTION = "software versions"
_SPECIES_MAP_SECTION = "species_map"


def make_relative_to(
    staging_path: pathlib.Path,
    install_path: pathlib.Path,
) -> pathlib.Path:
    assert staging_path.is_absolute() and install_path.is_absolute()

    for i, (s_part, i_part) in enumerate(
        zip(staging_path.parts, install_path.parts, strict=False),
    ):
        if s_part != i_part:
            break
    change_up = ("..",) * (len(staging_path.parts) - i)
    rel_path = change_up + install_path.parts[i:]
    return pathlib.Path(*rel_path)


@dataclass
class Config:
    domain: str  # User-specified domain (e.g., "main", "metazoa")
    release: str
    staging_path: pathlib.Path
    install_path: pathlib.Path
    species_dbs: dict[str, list[str]]
    align_names: Sequence[str]
    tree_names: Sequence[str]
    homologies: bool
    species_map: eti_species.SpeciesNameMap

    def __post_init__(self) -> None:
        self.staging_path = pathlib.Path(self.staging_path)
        self.install_path = pathlib.Path(self.install_path)

    def get_core_db_names(self) -> list[str]:
        names = [self.species_map.get_ensembl_db_prefix(n) for n in self.db_names]
        return [n for n in names if n]

    @property
    def staging_template_path(self) -> pathlib.Path:
        return self.staging_genomes / "coredb_templates"

    @property
    def db_names(self) -> Generator[str, None, None]:
        for species in self.species_dbs:
            yield self.species_map.get_ensembl_db_prefix(species)

    @property
    def staging_genomes(self) -> pathlib.Path:
        return self.staging_path / _GENOMES_NAME

    @property
    def install_genomes(self) -> pathlib.Path:
        return self.install_path / _GENOMES_NAME

    @property
    def staging_homologies(self) -> pathlib.Path:
        return self.staging_path / _COMPARA_NAME / _HOMOLOGIES_NAME

    @property
    def install_homologies(self) -> pathlib.Path:
        return self.install_path / _COMPARA_NAME / _HOMOLOGIES_NAME

    @property
    def staging_aligns(self) -> pathlib.Path:
        return self.staging_path / _COMPARA_NAME / _ALIGNS_NAME

    @property
    def install_aligns(self) -> pathlib.Path:
        return self.install_path / _COMPARA_NAME / _ALIGNS_NAME

    def to_dict(self, relative_paths: bool = True) -> dict[str, dict[str, str]]:
        """returns cfg as a dict"""
        if not self.species_dbs:
            msg = "no db names"
            raise ValueError(msg)

        if not relative_paths:
            staging_path = str(self.staging_path)
            install_path = str(self.install_path)
        else:
            staging_path = "."
            install_path = str(make_relative_to(self.staging_path, self.install_path))

        data = {
            "remote path": {"domain": str(self.domain)},
            "local path": {
                "staging_path": staging_path,
                "install_path": install_path,
            },
            "release": {"release": self.release},
            "compara": {},
        }

        if self.align_names:
            data["compara"]["align_names"] = "".join(self.align_names)
        if self.tree_names:
            data["compara"]["tree_names"] = "".join(self.tree_names)

        if self.homologies:
            data["compara"]["homologies"] = ""

        if not data["compara"]:
            data.pop("compara")

        for db_name in self.species_dbs:
            data[db_name] = {"db": "core"}

        return data

    def write(self) -> None:
        """writes a ini to staging_path/DOWNLOADED_CONFIG_NAME

        Notes
        -----
        Updates value for staging_path to '.', and install directories to be
        relative to staging_path.
        """
        parser = configparser.ConfigParser()
        cfg = self.to_dict()
        for section, settings in cfg.items():
            parser.add_section(section)
            for option, val in settings.items():
                parser.set(section, option=option, value=val)

        # add the species map section (maps between db name, common name,
        # abbrev etc..)
        parser.add_section(_SPECIES_MAP_SECTION)
        subset = self.species_map.get_subset(list(self.species_dbs))
        for label, value in subset.for_storage().items():
            parser.set(_SPECIES_MAP_SECTION, option=label, value=value)

        self.staging_path.mkdir(parents=True, exist_ok=True)
        with (self.staging_path / DOWNLOADED_CONFIG_NAME).open(mode="w") as out:
            parser.write(out, space_around_delimiters=True)


@dataclass
class InstalledConfig:
    release: str
    install_path: pathlib.Path
    software_versions: dict[str, str]
    species_map: eti_species.SpeciesNameMap

    def __hash__(self) -> int:
        return id(self)

    def __post_init__(self) -> None:
        self.install_path = pathlib.Path(self.install_path)

    @property
    def compara_path(self) -> pathlib.Path:
        return self.install_path / _COMPARA_NAME

    @property
    def homologies_path(self) -> pathlib.Path:
        return self.compara_path / _HOMOLOGIES_NAME

    @property
    def aligns_path(self) -> pathlib.Path:
        return self.compara_path / _ALIGNS_NAME

    @property
    def genomes_path(self) -> pathlib.Path:
        return self.install_path / _GENOMES_NAME

    def installed_genome(self, species: str) -> pathlib.Path:
        db_name = self.species_map.get_genome_name(species, level="raise")
        return self.genomes_path / db_name

    def list_genomes(self) -> list[str]:
        """returns list of installed genomes"""
        return [
            p.name for p in self.genomes_path.glob("*") if p.name in self.species_map
        ]

    def path_to_alignment(self, pattern: str, suffix: str) -> pathlib.Path | None:
        """returns the full path to alignment matching the name

        Parameters
        ----------
        pattern
            glob pattern for the Ensembl alignment name
        """
        if eti_util.contains_glob_pattern(pattern):
            align_dirs = [
                d
                for d in self.aligns_path.glob("*")
                if fnmatch.fnmatch(d.name, pattern)
            ]
        elif pattern:
            align_dirs = [d for d in self.aligns_path.glob("*") if pattern in d.name]
        else:
            align_dirs = None
        if not align_dirs:
            return None

        if len(align_dirs) > 1:
            msg = f"{pattern!r} matches too many directories in {self.aligns_path} {align_dirs}"
            raise ValueError(
                msg,
            )

        align_dir = align_dirs[0]
        if not list(align_dir.glob(f"*{suffix}")):
            msg = f"{align_dir} does not contain file with suffix {suffix}"
            raise FileNotFoundError(msg)
        return align_dir

    def get_version_table(self) -> "Table":
        """returns table of software versions used to make installation"""

        from cogent3 import make_table

        header = ["package", "version"]
        return make_table(
            header=header,
            data=sorted(self.software_versions.items()),
            index_name="package",
            title="Installation software versions:",
        )


def _get_dependency_versions() -> dict[str, str]:
    import re
    from importlib import metadata
    from importlib.util import find_spec

    # get the declared dependencies
    deps = {"ensembl_tui"}
    for pkg in metadata.requires("ensembl_tui"):
        if "extra" in pkg:
            continue
        if match := re.match(r"^[A-Za-z0-9_-]+", pkg):
            pkg = match.group(0).replace("-", "_")
            if find_spec(pkg):
                deps.add(pkg)

    return {pkg: metadata.version(pkg) for pkg in deps}


def write_installed_cfg(config: Config) -> eti_util.PathType:
    """writes an ini file under config.installed_path"""
    parser = configparser.ConfigParser()
    parser.add_section("release")
    parser.set("release", "release", config.release)
    # get the declared dependencies
    deps = _get_dependency_versions()
    parser.add_section(_VERSION_SECTION)
    for pkg, vers in deps.items():
        parser.set(_VERSION_SECTION, pkg, vers)

    # add the species map section (maps between db name, common name, abbrev etc..)
    parser.add_section(_SPECIES_MAP_SECTION)
    subset = config.species_map.get_subset(list(config.species_dbs))
    store_map = subset.for_storage()
    parser.set(_SPECIES_MAP_SECTION, "header", store_map.pop("header"))
    # now add the species, value the remainder comma separated
    for genome, value in store_map.items():
        parser.set(_SPECIES_MAP_SECTION, genome, value)

    outpath = config.install_path / INSTALLED_CONFIG_NAME
    outpath.parent.mkdir(parents=True, exist_ok=True)
    with outpath.open(mode="w") as out:
        parser.write(out)
    return outpath


def read_installed_cfg(path: eti_util.PathType) -> InstalledConfig:
    """reads an ini file under config.installed_path"""
    path = pathlib.Path(path).expanduser()
    parser = configparser.ConfigParser()
    path = (
        path if path.name == INSTALLED_CONFIG_NAME else (path / INSTALLED_CONFIG_NAME)
    )
    if not path.exists():
        eti_util.print_colour(f"{path!s} does not exist, exiting", colour="red")
        sys.exit(1)

    parser.read(path)
    release = parser.get("release", "release")
    if parser.has_section(_VERSION_SECTION):
        software_versions = dict(parser.items(_VERSION_SECTION))
    else:
        software_versions = {}

    if parser.has_section("species_map"):
        map_data = dict(parser.items("species_map"))
        sp_map = eti_species.SpeciesNameMap.from_storage(map_data)
    else:
        sp_map = eti_species.make_species_map(species_path=None)

    return InstalledConfig(
        release=release,
        install_path=path.parent,
        software_versions=software_versions,
        species_map=sp_map,
    )


def _standardise_path(
    path: eti_util.PathType,
    config_path: pathlib.Path,
) -> pathlib.Path:
    path = pathlib.Path(path).expanduser()
    return path if path.is_absolute() else (config_path / path).resolve()


def _pop_section(
    config: configparser.ConfigParser, section_name: str
) -> dict[str, str]:
    """Pop a section from config, returning its items as a dict"""
    # Get all items from the section
    data = dict(config.items(section_name))
    # Remove the section
    config.remove_section(section_name)

    return data


def _validate_and_resolve_domain(remote_section: dict[str, str]) -> tuple[str, str]:
    """Validate and resolve domain/host from remote path section.

    Handles backward compatibility by accepting both 'domain' and 'host'.
    If both present, 'domain' takes precedence with a deprecation warning.
    If only 'host' present, issues deprecation warning.

    Parameters
    ----------
    remote_section : dict
        The [remote path] section from config file

    Returns
    -------
    tuple[str, str]
        A tuple of (domain, host) where:
        - domain: The user-specified domain name (e.g., "main", "metazoa")
        - host: The resolved FTP hostname from site map (e.g., "ftp.ensembl.org")

    Raises
    ------
    SystemExit
        If validation fails (no domain/host specified, or invalid domain)
    """
    import warnings

    domain_value = remote_section.get("domain")
    host_value = remote_section.get("host")

    # Validation: at least one must be present
    if not domain_value and not host_value:
        msg = (
            "Config error: [remote path] section must contain either 'domain' or 'host'"
        )
        raise ValueError(msg)

    if host_value:
        msg = "The 'host' option in [remote path] is deprecated, use 'domain' instead."
        warnings.warn(
            msg,
            DeprecationWarning,
            stacklevel=3,
        )
        domain = domain_value or host_value
    else:
        domain = domain_value

    try:
        site_map = eti_site_map.get_site_map(domain)
    except KeyError as err:
        msg = f"Invalid domain '{domain}'. Available domains: {', '.join(sorted(eti_site_map.get_site_map_names()))}"
        raise ValueError(msg) from err

    # Resolve the actual FTP host from the site map
    site_map = eti_site_map.get_site_map(domain)
    host = site_map.site

    return domain, host


def read_config(
    *,
    config_path: pathlib.Path,
    species_map: eti_species.SpeciesNameMap | None = None,
    root_dir: pathlib.Path | None = None,
) -> Config:
    """returns ensembl release, local path, and db specifics from the provided
    config path"""
    from ensembl_tui._download import download_ensembl_tree, get_species_for_alignments

    if not config_path.exists():
        eti_util.print_colour(f"File not found {config_path.resolve()!s}", colour="red")
        sys.exit(1)

    parser = configparser.ConfigParser()

    with config_path.expanduser().open() as f:
        parser.read_file(f)

    if parser.has_section(_SPECIES_MAP_SECTION):
        # species map embedded in config takes precedence
        map_data = _pop_section(parser, _SPECIES_MAP_SECTION)
        sp_map = eti_species.SpeciesNameMap.from_storage(map_data)
    elif species_map is None:
        sp_map = eti_species.make_species_map(species_path=None)
    else:
        sp_map = species_map

    if root_dir is None:
        root_dir = config_path.parent

    release = _pop_section(parser, "release")["release"]
    remote_section = _pop_section(parser, "remote path")
    domain, host = _validate_and_resolve_domain(remote_section)
    site_map = eti_site_map.get_site_map(domain)
    # paths
    paths = _pop_section(parser, "local path")
    staging_path = _standardise_path(paths["staging_path"], root_dir)
    install_path = _standardise_path(paths["install_path"], root_dir)

    homologies = parser.has_option("compara", "homologies")
    align_names = []
    tree_names = []
    if parser.has_section("compara"):
        compara = _pop_section(parser, "compara")
        align_names = (
            [n.strip() for n in compara["align_names"].split(",")]
            if "align_names" in compara
            else []
        )
        tree_names = (
            [n.strip() for n in compara["tree_names"].split(",")]
            if "tree_names" in compara
            else []
        )

    species_dbs = {}
    for section in parser.sections():
        sec = _pop_section(parser, section)
        dbs = [db.strip() for db in sec["db"].split(",")]
        species_name = sp_map.get_genome_name(section, level="raise")
        species_dbs[species_name] = dbs

    # we also want homologies if we want alignments
    homologies = homologies or bool(align_names)

    if not species_dbs and (align_names or tree_names):
        found = set()
        if tree_names:
            # add all species in the tree to species_dbs
            for tree_name in tree_names:
                tree = download_ensembl_tree(
                    host=host,
                    release=release,
                    site_map=site_map,
                    tree_fname=tree_name,
                )
                if tree is None:
                    continue

                sp = set(
                    eti_species.species_from_ensembl_tree(tree, species_map=sp_map)
                )
                found |= sp

        if align_names:
            # add all species in the alignments to species_dbs
            sp = set(
                get_species_for_alignments(
                    host=host,
                    release=release,
                    site_map=site_map,
                    align_names=align_names,
                    species_map=sp_map,
                )
            )
            found |= sp

        species_dbs |= {n: ["core"] for n in found}

    return Config(
        domain=domain,
        release=release,
        staging_path=staging_path,
        install_path=install_path,
        species_dbs=species_dbs,
        align_names=align_names,
        tree_names=tree_names,
        homologies=homologies,
        species_map=sp_map,
    )
