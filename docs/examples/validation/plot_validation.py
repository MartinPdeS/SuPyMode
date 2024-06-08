"""
Propagation constant: DCFC
==========================
"""

# %%
# Imports
# ~~~~~~~
import numpy
from SuPyMode.workflow import Workflow, fiber_catalogue, Boundaries
import PyFiberModes
from PyFiberModes.fiber import load_fiber
from PyFiberModes.__future__ import get_normalized_LP_coupling
from MPSPlots.render2D import SceneList
import itertools

wavelength = 1550e-9
fiber_name = 'test_multimode_fiber'
scale_factor = 4

# %%
# Generating the fiber structure
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Here we define the cladding and fiber structure to model the problem
fiber = fiber_catalogue.load_fiber(fiber_name, wavelength=wavelength, remove_cladding=False)
fiber_list = [fiber]


# %%
# Defining the boundaries of the system
boundaries = [
    Boundaries(right='symmetric', bottom='symmetric'),
]

# %%
# Generating the computing workflow
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Workflow class to define all the computation parameters before initializing the solver
workflow = Workflow(
    fiber_list=fiber_list,          # List of fiber to be added in the mesh, the order matters.
    fusion_degree='auto',           # Degree of fusion of the structure if applicable.
    wavelength=wavelength,          # Wavelength used for the mode computation.
    resolution=180,                 # Number of point in the x and y axis [is divided by half if symmetric or anti-symmetric boundaries].
    x_bounds="left",                # Mesh x-boundary structure.
    y_bounds="top",                 # Mesh y-boundary structure.
    boundaries=boundaries,          # Set of symmetries to be evaluated, each symmetry add a round of simulation
    n_sorted_mode=7,                # Total computed and sorted mode.
    n_added_mode=6,                 # Additional computed mode that are not considered later except for field comparison [the higher the better but the slower].
    plot_geometry=True,             # Plot the geometry mesh before computation.
    debug_mode=1,                   # Print the iteration step for the solver plus some other important steps.
    auto_label=True,                # Auto labeling the mode. Label are not always correct and should be verified afterwards.
    itr_final=0.4,                  # Final value of inverse taper ratio to simulate
    index_scrambling=0,             # Scrambling of refractive index value in order to lift mode degeneracy [useful for some analysis]
    n_step=100
)

superset = workflow.get_superset()

superset.label_supermodes('LP01', 'LP21', 'LP02', 'LP03', 's')

superset.plot(plot_type='field').show()

itr_list = superset.model_parameters.itr_list

# %%
# Computing the analytical values using FiberModes solver.
initial_fiber = load_fiber(
    fiber_name=fiber_name,
    wavelength=wavelength,
    add_air_layer=False
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


def get_index_pyfibermodes(mode, itr_list, fiber):
    analytical = numpy.empty(itr_list.shape)

    for idx, itr in enumerate(itr_list):
        tapered_fiber = fiber.scale(factor=itr)
        analytical[idx] = tapered_fiber.get_effective_index(mode=mode)

    return analytical


for idx, mode in enumerate(['LP01', 'LP02', 'LP03']):
    color = f"C{idx}"

    supymode_mode = getattr(superset, mode)
    ax.add_scatter(
        x=itr_list,
        y=supymode_mode.index.data,
        label=str(supymode_mode),
        color='black',
        line_width=2,
        edge_color=color,
        marker_size=80,
        line_style='-',
        layer_position=2
    )

    analytical = get_index_pyfibermodes(mode=getattr(PyFiberModes, mode), itr_list=itr_list, fiber=initial_fiber)

    ax.add_line(
        x=itr_list,
        y=analytical,
        label=str(mode),
        line_style='-',
        line_width=2,
        color=color,
        layer_position=1
    )

_ = figure.show()


# %%
# Preparing the figure
figure = SceneList(unit_size=(12, 4))


ax = figure.append_ax(
    x_label='Inverse taper ratio',
    y_label='Normalized coupling',
    show_legend=True,
    font_size=18,
    tick_size=15,
    legend_font_size=18
)

# %%
# Computing the analytical values using FiberModes solver.
initial_fiber = load_fiber(
    fiber_name=fiber_name,
    wavelength=wavelength,
    add_air_layer=False
)


def get_normalized_coupling_pyfibermodes(mode_0, mode_1, itr_list, initial_fiber):
    analytical = numpy.empty(itr_list.shape)

    for idx, itr in enumerate(itr_list):
        tapered_fiber = initial_fiber.scale(factor=itr)

        analytical[idx] = get_normalized_LP_coupling(fiber=tapered_fiber, mode_0=mode_0, mode_1=mode_1)

    return analytical


for idx, (mode_0, mode_1) in enumerate(itertools.combinations(['LP01', 'LP02', 'LP03'], 2)):
    color = f"C{idx}"

    analytical = get_normalized_coupling_pyfibermodes(
        mode_0=getattr(PyFiberModes, mode_0),
        mode_1=getattr(PyFiberModes, mode_1),
        itr_list=itr_list[::2],
        initial_fiber=initial_fiber
    )

    ax.add_line(
        x=itr_list[::2],
        y=abs(analytical),
        label='Analytical',
        line_style='-',
        line_width=2,
        color=color,
        layer_position=1
    )

    simulation = getattr(superset, mode_0).normalized_coupling.get_values(getattr(superset, mode_1))

    ax.add_scatter(
        x=superset.model_parameters.itr_list,
        y=abs(simulation),
        color='black',
        line_width=2,
        edge_color=color,
        marker_size=80,
        line_style='-',
        layer_position=2,
        label=mode_0 + '-' + mode_1
    )

_ = figure.show()

# -
