import nox

_py_versions = range(11, 15)


@nox.session(python=[f"3.{v}" for v in _py_versions], venv_backend="uv")
def test(session):
    # on python >= 3.12 this will improve speed of test coverage a lot. Keyed
    # off the session interpreter, not the one running nox.
    if tuple(int(p) for p in session.python.split(".")) >= (3, 12):
        session.env["COVERAGE_CORE"] = "sysmon"
    session.install("-e", ".", "--group", "test")
    session.chdir("tests")
    session.run(
        "pytest",
        "-s",
        "-x",
        *session.posargs,
    )
