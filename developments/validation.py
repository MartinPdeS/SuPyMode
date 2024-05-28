"""
Propagation constant: DCFC
==========================
"""

# %%
# Imports
# ~~~~~~~
import numpy
from SuPyMode.workflow import Workflow, fiber_catalogue, Boundaries
from PyFiberModes.__future__ import get_normalized_LP_coupling
from PyFiberModes.fiber import load_fiber
from MPSPlots.render2D import SceneList
import PyFiberModes

wavelength = 1550e-9
fiber_name = 'SMF28'
scale_factor = 4


# %%
# Generating the fiber structure
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Here we define the cladding and fiber structure to model the problem
fiber = fiber_catalogue.load_fiber(fiber_name, wavelength=wavelength, remove_cladding=False)
fiber.structure_list[-1].scale(scale_factor)
fiber_list = [fiber]


# %%
# Defining the boundaries of the system
boundaries = [
    Boundaries(left='symmetric', top='symmetric'),
]


# %%
# Generating the computing workflow
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Workflow class to define all the computation parameters before initializing the solver
workflow = Workflow(
    fiber_list=fiber_list,          # List of fiber to be added in the mesh, the order matters.
    fusion_degree='auto',           # Degree of fusion of the structure if applicable.
    wavelength=wavelength,          # Wavelength used for the mode computation.
    resolution=150,                 # Number of point in the x and y axis [is divided by half if symmetric or anti-symmetric boundaries].
    x_bounds="right",                # Mesh x-boundary structure.
    y_bounds="bottom",                 # Mesh y-boundary structure.
    boundaries=boundaries,          # Set of symmetries to be evaluated, each symmetry add a round of simulation
    n_sorted_mode=3,                # Total computed and sorted mode.
    n_added_mode=3,                 # Additional computed mode that are not considered later except for field comparison [the higher the better but the slower].
    # plot_geometry=True,             # Plot the geometry mesh before computation.
    debug_mode=2,                   # Print the iteration step for the solver plus some other important steps.
    auto_label=True,                # Auto labeling the mode. Label are not always correct and should be verified afterwards.
    itr_final=0.3,                  # Final value of inverse taper ratio to simulate
    index_scrambling=0,             # Scrambling of refractive index value in order to lift mode degeneracy [useful for some analysis]
    n_step=100
)

superset = workflow.get_superset()
superset.label_supermodes('LP01', 'LP02', 'LP21')
# superset.label_supermodes('LP01', 'LP02', 'LP21', 'LP03', 'LP22', 'LP41')

superset.plot(plot_type='field').show()

itr_list = superset.itr_list[::8]

# %%
# Computing the analytical values using FiberModes solver.
pyfibermodes_fiber = load_fiber(
    fiber_name=fiber_name,
    wavelength=wavelength,
    add_air_layer=False
)

pyfibermodes_fiber = pyfibermodes_fiber.scale(scale_factor)

analytical = numpy.empty(itr_list.shape)
for idx, itr in enumerate(itr_list):
    print(idx, itr)
    _fiber = pyfibermodes_fiber.scale(factor=itr)
    analytical[idx] = get_normalized_LP_coupling(
        fiber=_fiber,
        mode_0=PyFiberModes.LP01,
        mode_1=PyFiberModes.LP02
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

ax.add_line(
    x=itr_list,
    y=abs(analytical),
    label='Analytical',
    line_style='-',
    line_width=2,
    color='red',
    layer_position=1
)

simulation = superset.LP01.normalized_coupling.get_values(superset.LP02).imag

ax.add_scatter(
    x=superset.itr_list,
    y=abs(simulation),
    label="SuPyMode",
    color='black',
    line_width=2,
    edge_color='blue',
    marker_size=80,
    line_style='-',
    layer_position=2
)

_ = figure.show()
