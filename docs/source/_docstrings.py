"""Normalize appended plotting-parameter documentation before Napoleon parses it."""

from inspect import cleandoc


def normalize_plot_docstrings(app, what, name, obj, options, lines):
    """Normalize independently indented sections appended by plotting decorators."""
    for index, line in enumerate(lines):
        if line.strip() != "Extra Parameters":
            continue
        main = cleandoc("\n".join(lines[:index]))
        extra = cleandoc("\n".join(lines[index:]))
        lines[:] = main.splitlines() + [""] + extra.splitlines() + [""]
        break
