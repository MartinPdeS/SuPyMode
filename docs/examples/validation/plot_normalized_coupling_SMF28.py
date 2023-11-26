"""
Normalized coupling: SMF28
==========================
"""

# %%
# Imports
# ~~~~~~~
import numpy
from SuPyMode.workflow import Workflow, fiber_catalogue, Boundaries2D
from MPSPlots.render2D import SceneList
from PyFiberModes.fiber import load_fiber
from PyFiberModes.__future__ import get_normalized_LP_coupling
from PyFiberModes import LP01, LP02
from SuPyMode import load_superset
wavelength = 1550e-9


# %%
# Generating the fiber structure
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Here we define the cladding and fiber structure to model the problem
fiber = fiber_catalogue.SMF28(wavelength=wavelength)


fiber = fiber.scale(factor=10)

fiber_list = [fiber]

# %%
# Defining the boundaries of the system
boundaries = [
    Boundaries2D(right='symmetric', top='symmetric'),
]


# %%
# Generating the computing workflow
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Workflow class to define all the computation parameters before initializing the solver
workflow = Workflow(
    fiber_list=fiber_list,          # List of fiber to be added in the mesh, the order matters.
    fusion_degree=None,             # Degree of fusion of the structure if applicable.
    wavelength=wavelength,          # Wavelength used for the mode computation.
    resolution=15,                 # Number of point in the x and y axis [is divided by half if symmetric or anti-symmetric boundaries].
    x_bounds=[-300e-6, 0],          # Mesh x-boundary structure.
    y_bounds=[-300e-6, 0],          # Mesh y-boundary structure.
    boundaries=boundaries,          # Set of symmetries to be evaluated, each symmetry add a round of simulation
    n_sorted_mode=4,                # Total computed and sorted mode.
    n_added_mode=8,                 # Additional computed mode that are not considered later except for field comparison [the higher the better but the slower].
    plot_geometry=False,             # Plot the geometry mesh before computation.
    debug_mode=True,                # Print the iteration step for the solver plus some other important steps.
    auto_label=True,                # Auto labeling the mode. Label are not always correct and should be verified afterwards.
    itr_final=0.2,                  # Final value of inverse taper ratio to simulate
    n_step=100,
    index_scrambling=0              # Scrambling of refractive index value in order to lift mode degeneracy [useful for some analysis]
)


superset = workflow.get_superset()
superset.save_instance(filename='_smf28_testing', directory='auto')
# superset.plot('normalized-coupling').show()

superset = load_superset(filename='_smf28_testing', directory='auto')

sub_sampling = 15
itr_list_sub = superset.itr_list[::sub_sampling]

# %%
# Analytical LP01 solution
# ~~~~~~~~~~~~~~~~~~~
# Computing the mode propagation constant using analytical solver.
smf28 = load_fiber(
    fiber_name='SMF28',
    wavelength=wavelength,
    add_air_layer=False
)

smf28 = smf28.scale(factor=10)

analytic = numpy.empty(itr_list_sub.shape)
for idx, itr in enumerate(itr_list_sub):
    _fiber = smf28.scale(factor=itr)
    analytic[idx] = get_normalized_LP_coupling(fiber=_fiber, mode_0=LP01, mode_1=LP02)

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

x = itr_list_sub
y = analytic
ax.add_line(
    x=x,
    y=y,
    label="Analytic",
    line_style='-',
    line_width=2,
    color='red',
)

x = itr_list_sub
y = abs(superset.LP01.normalized_coupling.get_values(superset.LP02))[::sub_sampling]
ax.add_scatter(
    x=x,
    y=y,
    label="SuPyMode",
    color='black',
    line_width=2,
    edge_color='blue',
    marker_size=80,
)

_ = figure.show()
