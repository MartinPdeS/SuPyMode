#!/usr/bin/env python
# -*- coding: utf-8 -*-

import logging

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




# -
