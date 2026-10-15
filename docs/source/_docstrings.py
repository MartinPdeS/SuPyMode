"""Normalize appended plotting-parameter documentation before Napoleon parses it."""


def normalize_plot_docstrings(app, what, name, obj, options, lines):
    """Dedent the Extra Parameters section appended by MPSPlots decorators."""
    for index, line in enumerate(lines):
        if line.strip() != "Extra Parameters":
            continue
        indentation = len(line) - len(line.lstrip())
        if indentation:
            prefix = " " * indentation
            lines[index:] = [
                item[indentation:] if item.startswith(prefix) else item
                for item in lines[index:]
            ]
        break
