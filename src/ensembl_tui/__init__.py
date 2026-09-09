"""Ensembl terminal user interface tools"""

from importlib.metadata import version
from warnings import filterwarnings

filterwarnings("ignore", message=".*MPI")
filterwarnings("ignore", message="Can't drop database.*")
filterwarnings("ignore", message="A worker stopped while some jobs.*")

# the version is declared in pyproject.toml, we read it back from the
# installed distribution metadata so there is a single source of truth
__version__ = version("ensembl_tui")
