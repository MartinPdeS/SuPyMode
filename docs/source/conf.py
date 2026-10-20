"""SuPyMode documentation using the shared MartinPdeS-docs configuration."""

import os
from pathlib import Path
import sys
import tomllib
from importlib.metadata import distribution as package_distribution

from martinpdes_docs import configure


ROOT = Path(__file__).resolve().parents[2]
DOCS_SOURCE = Path(__file__).resolve().parent
sys.path.insert(0, str(DOCS_SOURCE))
installed_site = str(package_distribution("SuPyMode").locate_file(""))
if installed_site in sys.path:
    sys.path.remove(installed_site)
sys.path.insert(0, installed_site)

examples_sections = ["workflow", "alpha_profile", "validation", "extra"]

globals().update(
    configure(
        project="SuPyMode",
        package="SuPyMode",
        repository="MartinPdeS/SuPyMode",
        logo="_static/logo.png",
        favicon="_static/favicon.png",
        version=os.getenv("tag", "latest"),
        html_theme_options={
            "logo": {
                "text": "",
                "image": "_static/logo.png",
                "alt_text": "SuPyMode logo",
            },
        },
        sphinx_gallery_conf={
            "examples_dirs": [f"../examples/{name}" for name in examples_sections],
            "gallery_dirs": [f"gallery/{name}" for name in examples_sections],
            "image_scrapers": ("matplotlib", "pyvista"),
            "default_thumb_file": str(
                Path(__file__).resolve().parent / "_static" / "gallery-default.png"
            ),
            "within_subsection_order": "FileNameSortKey",
            "abort_on_example_error": True,
        },
    )
)

extensions += [
    "sphinx.ext.mathjax",
    "pyvista.ext.plot_directive",
    "sphinx.ext.autosectionlabel",
]
autodoc_default_options.update(
    {"members": True, "exclude-members": "__annotations__", "member-order": "bysource"}
)
napoleon_use_ivar = True
napoleon_custom_sections = [("Extra Parameters", "params_style")]
autosectionlabel_maxdepth = 1
html_context.update({"doc_path": "docs/source", "default_mode": "dark"})
intersphinx_mapping = {"numpy": ("https://numpy.org/devdocs", None)}

metadata = tomllib.loads((ROOT / "pyproject.toml").read_text(encoding="utf-8"))
version = release = metadata["project"]["version"]

from _docstrings import normalize_plot_docstrings


def setup(app):
    """Prepare decorator-added parameter sections before Napoleon processes them."""
    app.connect("autodoc-process-docstring", normalize_plot_docstrings, priority=100)
