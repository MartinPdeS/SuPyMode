#!/usr/bin/env python
# -*- coding: utf-8 -*-

import logging
from matplotlib.colors import LinearSegmentedColormap

try:
    from mayavi         import mlab
    from tvtk.tools     import visual
    from tvtk.api       import tvtk
    from tvtk.common    import configure_input_data
except ImportError:
    logging.warning('Mayavi package could not be loaded! Not 3D rendering available.')

CMAP   = 'jet'
RED    = (1,0,0)
BLACK  = (0,0,0)
BLUE   = (0,0,1)
YELLOW = (1,1,0)
WHITE  = (1,1,1)
GREY   = (0.3, 0.3, 0.3)


FieldMap = LinearSegmentedColormap.from_list('my_gradient', (
                 # Edit this gradient at https://eltos.github.io/gradient/#4C71FF-0025B3-000000-C7030D-FC4A53
                 (0.000, (0.298, 0.443, 1.000)),
                 (0.250, (0.000, 0.145, 0.702)),
                 (0.500, (0.000, 0.000, 0.000)),
                 (0.750, (0.780, 0.012, 0.051)),
                 (1.000, (0.988, 0.290, 0.325))))


# -
