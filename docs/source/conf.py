"""SuPyMode documentation using the shared MartinPdeS-docs configuration."""

from pathlib import Path
import sys
import tomllib

from martinpdes_docs import configure


ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(Path(__file__).resolve().parent))

gallery_sections = ["workflow", "alpha_profile", "validation", "extra"]

globals().update(
    configure(
        project="SuPyMode",
        package="SuPyMode",
        repository="MartinPdeS/SuPyMode",
        logo="_static/logo.png",
        favicon="_static/favicon.png",
        prefer_installed_package=True,
        gallery_sections=gallery_sections,
        gallery_default_thumb="_static/gallery-default.png",
        gallery_options={
            "image_scrapers": ("matplotlib", "pyvista"),
            "abort_on_example_error": True,
        },
        html_theme_options={
            "logo": {
                "text": "",
                "image": "_static/logo.png",
                "alt_text": "SuPyMode logo",
            },
        },
    )
)

extensions += ["pyvista.ext.plot_directive"]
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
