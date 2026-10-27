"""Keep explicit release metadata consistent without Git version discovery."""

import json
import re
import runpy

import pytest
from pathlib import Path

try:
    import tomllib
except ModuleNotFoundError:
    import tomli as tomllib


ROOT = Path(__file__).resolve().parents[1]


def conda_version(recipe):
    """Read a simple quoted or unquoted Conda package version scalar."""
    match = re.search(
        r"^  version: (?:(?P<quote>[\"'])(?P<quoted>[^\"']+)(?P=quote)|(?P<plain>[^\s#\"']+))\s*(?:#.*)?$",
        recipe,
        re.MULTILINE,
    )
    if not match:
        raise ValueError("Conda recipe has no valid version scalar")
    return match.group("quoted") or match.group("plain")


@pytest.mark.parametrize("value", ["2.4.12", '"2.4.12"', "'2.4.12'"])
def test_conda_version_accepts_yaml_string_styles(value):
    assert conda_version(f"package:\n  version: {value}\n") == "2.4.12"


def test_conda_version_rejects_mismatched_quotes():
    with pytest.raises(ValueError, match="valid version"):
        conda_version("package:\n  version: \"2.4.12'\n")


def test_explicit_release_metadata_agrees():
    """Source imports, distribution metadata, and Conda declare one version."""
    metadata = tomllib.loads((ROOT / "pyproject.toml").read_text())
    project = metadata["project"]
    assert "version" not in project.get("dynamic", [])
    version_file = ROOT / metadata["tool"]["mps-maintenance"]["version-file"]
    package = runpy.run_path(str(version_file))
    assert package["__version__"] == package["version"] == project["version"]
    assert package["version_tuple"] == tuple(
        int(part) for part in project["version"].split(".")
    )
    recipe = (ROOT / "conda.recipe/meta.yaml").read_text()
    assert conda_version(recipe) == project["version"]
    assert metadata["tool"]["mps-maintenance"]["version-method"] == "literal"
    dependencies = metadata["build-system"]["requires"] + project["dependencies"]
    assert not any("setuptools-scm" in item.replace("_", "-") for item in dependencies)

    zenodo = json.loads((ROOT / ".zenodo.json").read_text())
    assert zenodo["version"] == project["version"]
