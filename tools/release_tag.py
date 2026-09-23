#!/usr/bin/env python3
"""Create a SuPyMode release commit and annotated semantic-version tag."""
import argparse
import os
from pathlib import Path
import re
import subprocess
import sys

ROOT = Path(__file__).resolve().parents[1]
CONDA_RECIPE_PATH = ROOT / "conda.recipe" / "meta.yaml"
PYPROJECT_PATH = ROOT / "pyproject.toml"
VERSION_FILE = ROOT / "SuPyMode" / "_version.py"
TAG_PATTERN = re.compile(r"v(?P<version>(?:0|[1-9]\d*)\.(?:0|[1-9]\d*)\.(?:0|[1-9]\d*))$")

def run(*command: str, capture_output: bool = False, env: dict[str, str] | None = None) -> str:
    result = subprocess.run(command, cwd=ROOT, check=True, text=True, capture_output=capture_output, env=env)
    return result.stdout.strip() if capture_output else ""

def validate_tag(tag: str) -> str:
    """Return the PEP 440 version represented by a release tag."""
    match = TAG_PATTERN.fullmatch(tag)
    if match is None:
        raise ValueError("tag must use vMAJOR.MINOR.PATCH, for example v2.2.3")
    return match.group("version")


def update_conda_recipe(version: str) -> None:
    """Write the release version into the Conda recipe."""
    recipe = CONDA_RECIPE_PATH.read_text(encoding="utf-8")
    updated_recipe, replacements = re.subn(
        r"(?m)^(?P<prefix>\s*version:\s*)\S+(?P<suffix>\s*(?:#.*)?)$",
        rf"\g<prefix>{version}\g<suffix>",
        recipe,
        count=1,
    )
    if replacements != 1:
        raise RuntimeError("could not find exactly one package version in conda.recipe/meta.yaml")
    CONDA_RECIPE_PATH.write_text(updated_recipe, encoding="utf-8")


def update_pyproject(version: str) -> None:
    """Write the release version into the Python project metadata."""
    project = PYPROJECT_PATH.read_text(encoding="utf-8")
    updated_project, replacements = re.subn(
        r'(?ms)(^\[project\]\n.*?^version\s*=\s*)["\'][^"\']+["\']',
        rf'\g<1>"{version}"',
        project,
        count=1,
    )
    if replacements != 1:
        raise RuntimeError("could not find exactly one project version in pyproject.toml")
    PYPROJECT_PATH.write_text(updated_project, encoding="utf-8")


def generate_version_file(version: str) -> None:
    """Generate ``_version.py`` through the configured SCM-versioning tool."""
    environment = os.environ.copy()
    environment["SETUPTOOLS_SCM_PRETEND_VERSION"] = version
    run(sys.executable, "-m", "setuptools_scm", "--force-write-version-files", env=environment)
    if not VERSION_FILE.exists():
        raise RuntimeError("SCM versioning did not generate SuPyMode/_version.py")


def create_release(tag: str, version: str) -> None:
    """Write metadata, commit it, and create the annotated release tag."""
    update_conda_recipe(version)
    update_pyproject(version)
    generate_version_file(version)
    run("git", "add", "conda.recipe/meta.yaml", "pyproject.toml", "SuPyMode/_version.py")
    run("git", "commit", "-m", f"Release {tag}")
    run("git", "tag", "-a", tag, "-m", f"Release {tag}")


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("tag", help="annotated release tag, for example v2.2.3")
    tag = parser.parse_args().tag
    try:
        if run("git", "status", "--porcelain", capture_output=True):
            raise RuntimeError("working tree is not clean; commit or stash changes before creating a release tag")
        if run("git", "tag", "--list", tag, capture_output=True):
            raise RuntimeError(f"tag {tag} already exists")
        create_release(tag, validate_tag(tag))
    except (RuntimeError, ValueError, subprocess.CalledProcessError) as error:
        print(f"release aborted: {error}", file=sys.stderr)
        return 1
    print(f"created release commit and annotated tag {tag}")
    return 0

if __name__ == "__main__":
    raise SystemExit(main())
