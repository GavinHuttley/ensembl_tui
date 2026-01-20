"""Ensembl terminal user interface tools"""

from warnings import filterwarnings

filterwarnings("ignore", message=".*MPI")
filterwarnings("ignore", message="Can't drop database.*")
filterwarnings("ignore", message="A worker stopped while some jobs.*")

__version__ = "0.7.5"
