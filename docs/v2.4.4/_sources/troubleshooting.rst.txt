Troubleshooting
===============

Installation and imports
------------------------

Use the same Python interpreter for installation and execution::

    python -m pip install --upgrade SuPyMode
    python -c "import SuPyMode; print(SuPyMode.__version__)"

A missing native module often means that no compatible wheel was installed or
a source build did not finish. Confirm your Python version and platform before
changing the simulation. Source builds need the repository's compiler and native
library dependencies.

Plots on a server
-----------------

Use ``MPLBACKEND=Agg`` for Matplotlib on a headless machine. Save figures to files
when no interactive display is available.

Three-dimensional PyVista rendering additionally needs a working off-screen
rendering backend. A failure to render does not by itself establish a problem
with the calculated modes.

Convergence and mode labels
---------------------------

Compare mesh resolutions and taper sample counts. Increase the number of added
modes when tracking becomes ambiguous, and inspect fields rather than relying
only on automatic labels. See :doc:`theory` and the validation examples.
