"""Ensembl terminal user interface tools"""

import os
from warnings import filterwarnings

filterwarnings("ignore", message=".*MPI")
filterwarnings("ignore", message="Can't drop database.*")
filterwarnings("ignore", message="A worker stopped while some jobs.*")

os.environ["COGENT3_NEW_TYPE"] = "1"
__version__ = "0.2.1"
