"""SuPyMode documentation using the shared MartinPdeS-docs configuration."""

from pathlib import Path
import sys

from martinpdes_docs import configure


package_name = "SuPyMode"
docs_source = Path(__file__).resolve().parent
sys.path.insert(0, str(docs_source))

gallery_sections = ["workflow", "alpha_profile", "validation", "extra"]

globals().update(
    configure(
        project=package_name,
        package=package_name,
        version="auto",
        author="Martin Poinsinet de Sivry-Houle",
        repository=f"MartinPdeS/{package_name}",
        logo="_static/logo.png",
        favicon="_static/favicon.png",
        prefer_installed_package=True,
        docs_root=docs_source.parent,
        profile="scientific-gallery",
        gallery_sections=gallery_sections,
        gallery_default_thumb=str(docs_source / "_static" / "gallery-default.png"),
        gallery_options={
            "image_scrapers": ("matplotlib", "pyvista"),
            "abort_on_example_error": True,
        },
        extension_additions=["pyvista.ext.plot_directive"],
        logo_text="",
        logo_alt_text="SuPyMode logo",
        autodoc_default_options={
            "members": True,
            "exclude-members": "__annotations__",
            "member-order": "bysource",
        },
    )
)

napoleon_custom_sections = [("Extra Parameters", "params_style")]
autosectionlabel_maxdepth = 1


from _docstrings import normalize_plot_docstrings


def setup(app):
    """Prepare decorator-added parameter sections before Napoleon processes them."""
    app.connect("autodoc-process-docstring", normalize_plot_docstrings, priority=100)
