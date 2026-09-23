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
Do not commit generated documentation, compiled artifacts, caches, or
`SuPyMode/_version.py`; SCM versioning generates that file for releases.

Before opening a pull request, run `make quality`, `make test`, and
`make release-check`. Keep changes focused and report platform/compiler details
for native changes.

Releases use semantic versioning. `make release patch`, `make release minor`,
and `make release major` create and push the release commit and exact tag.
`make tag VERSION=vX.Y.Z` creates a local tag without pushing.
