"""Verify the installed release commands against repository metadata."""

import re
import subprocess
import sys
from pathlib import Path


REPOSITORY_ROOT = Path(__file__).resolve().parents[1]


def run_tool(name: str, *arguments: str) -> subprocess.CompletedProcess:
    """Run the installed release tool using the current test interpreter."""
    return subprocess.run(
        [sys.executable, "-m", "martinpdes_maintenance.cli", name, *arguments],
        cwd=REPOSITORY_ROOT,
        capture_output=True,
        text=True,
        check=False,
    )


def test_release_versions_are_aligned() -> None:
    """The shared check validates every configured release artifact."""
    result = run_tool("check")
    assert result.returncode == 0, result.stderr
    assert "pyproject.toml" in result.stdout
    assert ".zenodo.json" in result.stdout
    assert "conda.recipe/meta.yaml" in result.stdout
    assert "SuPyMode/_version.py" in result.stdout


def test_release_check_rejects_unexpected_version() -> None:
    """An explicitly incorrect expected version fails the compatibility command."""
    result = run_tool("check", "--version", "v999999.0.0")
    assert result.returncode == 1
    assert "MISMATCH" in result.stdout


def test_next_release_prints_semantic_tag() -> None:
    """Make receives exactly one usable tag from the shared next command."""
    result = run_tool("next", "patch")
    assert result.returncode == 0, result.stderr
    assert re.fullmatch(r"v[0-9]+\.[0-9]+\.[0-9]+\n", result.stdout)
