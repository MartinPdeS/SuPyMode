"""
Propagation constant: DCFC
==========================
"""

# %%
# Imports
# ~~~~~~~
import numpy
from SuPyMode.workflow import Workflow, fiber_catalogue, Boundaries2D
from PyFiberModes import LP01
from PyFiberModes.fiber import load_fiber
from MPSPlots.render2D import SceneList

wavelength = 1550e-9
fiber_name = 'DCF1300S_33'

# %%
# Generating the fiber structure
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Here we define the cladding and fiber structure to model the problem
fiber_list = [
    fiber_catalogue.load_fiber(fiber_name, wavelength=wavelength)
]


# %%
# Defining the boundaries of the system
boundaries = [
    Boundaries2D(right='symmetric', bottom='symmetric'),
    Boundaries2D(right='symmetric', bottom='anti-symmetric')
]

# %%
# Generating the computing workflow
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Workflow class to define all the computation parameters before initializing the solver
workflow = Workflow(
    fiber_list=fiber_list,          # List of fiber to be added in the mesh, the order matters.
    fusion_degree=None,             # Degree of fusion of the structure if applicable.
    wavelength=wavelength,          # Wavelength used for the mode computation.
    resolution=50,                  # Number of point in the x and y axis [is divided by half if symmetric or anti-symmetric boundaries].
    x_bounds="left",                # Mesh x-boundary structure.
    y_bounds="top",                 # Mesh y-boundary structure.
    boundaries=boundaries,          # Set of symmetries to be evaluated, each symmetry add a round of simulation
    n_sorted_mode=6,                # Total computed and sorted mode.
    n_added_mode=4,                 # Additional computed mode that are not considered later except for field comparison [the higher the better but the slower].
    plot_geometry=True,             # Plot the geometry mesh before computation.
    debug_mode=0,                   # Print the iteration step for the solver plus some other important steps.
    auto_label=True,                # Auto labeling the mode. Label are not always correct and should be verified afterwards.
    itr_final=0.2,                  # Final value of inverse taper ratio to simulate
    index_scrambling=0,             # Scrambling of refractive index value in order to lift mode degeneracy [useful for some analysis]
    n_step=50
)

superset = workflow.get_superset()
itr_list = superset.itr_list

# %%
# Computing the analytical values using FiberModes solver.
dcf_fiber = load_fiber(
    fiber_name=fiber_name,
    wavelength=wavelength,
    add_air_layer=True
)

# %%
# Preparing the figure
figure = SceneList(unit_size=(12, 4))

ax = figure.append_ax(
    x_label='Inverse taper ratio',
    y_label='Effective index',
    show_legend=True,
    font_size=18,
    tick_size=15,
    legend_font_size=18
)

pyfibermodes_mode = LP01
supymode_mode = superset.LP01

analytical = numpy.empty(itr_list.shape)
for idx, itr in enumerate(itr_list):
    _fiber = dcf_fiber.scale(factor=itr)
    analytical[idx] = _fiber.get_effective_index(mode=pyfibermodes_mode)

ax.add_line(
    x=itr_list,
    y=analytical,
    label=pyfibermodes_mode,
    line_style='-',
    line_width=2,
    color='red',
    layer_position=1
)

ax.add_scatter(
    x=itr_list,
    y=supymode_mode.index.get_values(),
    label=supymode_mode,
    color='black',
    line_width=2,
    edge_color='blue',
    marker_size=80,
    line_style='-',
    layer_position=2
)

_ = figure.show()


# -
