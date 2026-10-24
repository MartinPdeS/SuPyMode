# Contributing to SuPyMode

Contributions are welcome: focused bug fixes, numerical validation cases,
documentation corrections, and feature work.

SuPyMode requires Python 3.11 or newer, CMake, a C++ compiler, and pybind11.
After changing native sources, rebuild the editable installation before tests.

```bash
git clone https://github.com/MartinPdeS/SuPyMode.git
cd SuPyMode
make editable
make test
```

Use public Python APIs in tests, document physical conventions and numerical
limits, and keep gallery examples small enough for the documentation build.
Do not commit generated documentation, compiled artifacts, or caches. The release
helper maintains the tracked `SuPyMode/_version.py` alongside `pyproject.toml`,
the Conda recipe, and Zenodo metadata.

Before opening a pull request, run `make quality`, `make test`, and
`make release-check`. Keep changes focused and report platform/compiler details
for native changes.

Releases use semantic versioning. `make release patch`, `make release minor`,
and `make release major` create and push the release commit and exact tag.
`make tag VERSION=vX.Y.Z` creates a local tag without pushing.

## Shared release tooling

Install the pinned release tool in the interpreter used by Make, without
rebuilding this package:

```sh
python -m pip install "martinpdes-maintenance==0.2.2"
make release-check PYTHON=python
```

The same dependency is available through the project's `release` extra.
`tools/check_release.py`, `tools/next_release_version.py`, and
`tools/release_tag.py` are compatibility wrappers around `mps-maintenance`.
Package-specific metadata paths live in `[tool.mps-maintenance]` in
`pyproject.toml`. Release tooling is not a runtime dependency.

Preview the next release's metadata changes from a clean, full clone:

```sh
python tools/release_tag.py vX.Y.Z --dry-run
```

Fetch tags before selecting a version. Existing `make tag` and `make release`
commands keep their local-tag and push behavior respectively. The shared tool
requires an unused stable tag and a clean working tree. It stages only configured
release metadata and never pushes itself. If committing or tagging fails,
inspect `git status` and `git log -1` before recovering or retrying.

Package versions are explicit in `[project].version`. Builds do not derive
versions from Git tags or require setuptools-scm. Release commands update the
project version, package version file, Conda recipe, and Zenodo metadata together.
