"""Tests for repository workflow helpers."""

import importlib.util
from pathlib import Path
from types import ModuleType


REPOSITORY_ROOT = Path(__file__).resolve().parents[1]
TOOLS_DIRECTORY = REPOSITORY_ROOT / "tools"


def load_tool_module(name: str) -> ModuleType:
    """Load a repository-only workflow helper without packaging it."""
    path = TOOLS_DIRECTORY / f"{name}.py"
    specification = importlib.util.spec_from_file_location(f"supymode_test_{name}", path)
    if specification is None or specification.loader is None:
        raise RuntimeError(f"could not load workflow helper: {path}")
    module = importlib.util.module_from_spec(specification)
    specification.loader.exec_module(module)
    return module


check_release = load_tool_module("check_release")
next_release_version = load_tool_module("next_release_version")


def test_release_versions_are_aligned() -> None:
    """Package release artifacts declare one version."""
    declared = check_release.versions()
    assert len(set(declared.values())) == 1


def test_normalized_version() -> None:
    """Release checks accept plain and Git-style version strings."""
    assert check_release.normalized_version("2.2.2") == "2.2.2"
    assert check_release.normalized_version("v2.2.2") == "2.2.2"


def test_latest_version_ignores_non_semantic_tags() -> None:
    """Only three-component semantic version tags determine the next release."""
    assert next_release_version.TAG_PATTERN.fullmatch("v2.0.3.3") is None