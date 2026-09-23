#!/usr/bin/env python3
"""Check that SuPyMode release metadata agrees across project files."""
import argparse
from pathlib import Path
import re
import subprocess
import sys
import tomllib

ROOT = Path(__file__).resolve().parents[1]
CONDA_RECIPE_PATH = ROOT / "conda.recipe" / "meta.yaml"
VERSION_FILE = ROOT / "SuPyMode" / "_version.py"


def normalized_version(version: str) -> str:
    """Return a version without an optional Git tag prefix."""
    return version.removeprefix("v")


def extract(pattern: str, path: Path, description: str) -> str:
    """Extract one version from a text file."""
    match = re.search(pattern, path.read_text(encoding="utf-8"), flags=re.MULTILINE)
    if match is None:
        raise RuntimeError(f"could not find {description} in {path.relative_to(ROOT)}")
    return match.group(1)


def versions() -> dict[str, str]:
    """Return release versions declared by each project artifact."""
    project = tomllib.loads((ROOT / "pyproject.toml").read_text(encoding="utf-8"))
    return {
        "pyproject.toml": str(project["project"]["version"]),
        "conda.recipe/meta.yaml": extract(
            r"^\s*version:\s*['\"]?([^'\"\s]+)", CONDA_RECIPE_PATH, "Conda version"
        ),
        "SuPyMode/_version.py": extract(
            r"^__version__\s*=\s*version\s*=\s*['\"]([^'\"]+)",
            VERSION_FILE,
            "generated package version",
        ),
    }


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--version")
    arguments = parser.parse_args()
    try:
        declared = versions()
    except (KeyError, RuntimeError) as error:
        print(f"release check failed: {error}", file=sys.stderr)
        return 1
    expected = normalized_version(arguments.version) if arguments.version else declared["pyproject.toml"]
    failures = []
    for source, version in declared.items():
        status = "OK" if version == expected else "MISMATCH"
        print(f"{status:8} {source}: {version}")
        if version != expected:
            failures.append(f"{source} declares {version}; expected {expected}")

    tag = f"v{expected}"
    tag_exists = subprocess.run(
        ["git", "rev-parse", "--verify", "--quiet", f"refs/tags/{tag}"],
        cwd=ROOT,
        check=False,
        stdout=subprocess.DEVNULL,
    ).returncode == 0
    print(f"{'OK' if tag_exists else 'INFO':8} Git tag: {tag} {'exists' if tag_exists else 'not created'}")

    if failures:
        for failure in failures:
            print(f"ERROR: {failure}", file=sys.stderr)
        return 1
    print(f"Release metadata consistently declares {expected}.")
    return 0

if __name__ == "__main__":
    raise SystemExit(main())
