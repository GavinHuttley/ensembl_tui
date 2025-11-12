import os
import sys

import nox

_py_versions = range(10, 14)

# on python >= 3.12 this will improve speed of test coverage a lot
if sys.version_info >= (3, 12):
    os.environ["COVERAGE_CORE"] = "sysmon"


@nox.session(python=[f"3.{v}" for v in _py_versions])
def test(session):
    session.install("-e.[test]")
    session.run("pip", "list")
    session.chdir("tests")
    session.run(
        "pytest",
        "-s",
        "-x",
        *session.posargs,
    )
