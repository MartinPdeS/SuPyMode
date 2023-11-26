"""
Propagation constant: SMF28
===========================
"""

# %%
# Imports
# ~~~~~~~
import numpy
from SuPyMode.workflow import Workflow, fiber_catalogue, Boundaries2D, configuration
from MPSPlots.render2D import SceneList
from PyFiberModes.fiber import load_fiber
from PyFiberModes import LP01

wavelength = 1550e-9

# %%
# Generating the fiber structure
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Here we define the cladding and fiber structure to model the problem
clad_structure = configuration.ring.FusedProfile_01x01

fiber_list = [fiber_catalogue.SMF28(wavelength=wavelength)]


# %%
# Defining the boundaries of the system
boundaries = [
    Boundaries2D(right='symmetric', bottom='symmetric'),
    Boundaries2D(right='symmetric', bottom='anti-symmetric')
]

# %%
# SuPyMode LP01 solution
# ~~~~~~~~~~~~~~~~~~~
# Computing the mode propagation constant using analytical solver.
workflow = Workflow(
    fiber_list=fiber_list,          # List of fiber to be added in the mesh, the order matters.
    clad_structure=clad_structure,  # Cladding structure, if None provided then no cladding is set.
    fusion_degree=None,             # Degree of fusion of the structure if applicable.
    wavelength=wavelength,          # Wavelength used for the mode computation.
    resolution=100,                 # Number of point in the x and y axis [is divided by half if symmetric or anti-symmetric boundaries].
    x_bounds="centering-left",      # Mesh x-boundary structure.
    y_bounds="centering-top",       # Mesh y-boundary structure.
    boundaries=boundaries,          # Set of symmetries to be evaluated, each symmetry add a round of simulation
    n_sorted_mode=6,                # Total computed and sorted mode.
    n_added_mode=4,                 # Additional computed mode that are not considered later except for field comparison [the higher the better but the slower].
    plot_geometry=False,            # Plot the geometry mesh before computation.
    debug_mode=True,                # Print the iteration step for the solver plus some other important steps.
    auto_label=True,                # Auto labeling the mode. Label are not always correct and should be verified afterwards.
    itr_final=0.05,                 # Final value of inverse taper ratio to simulate
    index_scrambling=0              # Scrambling of refractive index value in order to lift mode degeneracy [useful for some analysis]
)

superset = workflow.get_superset()
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

analytic = numpy.empty(itr_list_sub.shape)
for idx, itr in enumerate(itr_list_sub):
    _fiber = smf28.scale(factor=itr)
    analytic[idx] = _fiber.get_effective_index(mode=LP01)


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
    x=superset.itr_list,
    y=superset.LP01.index.get_values(),
    label='LP01: SuPyMode',
    line_style='-',
    line_width=2,
    color='red',
    layer_position=1
)

sub_sampling = 20
ax.add_scatter(
    x=itr_list_sub,
    y=analytic,
    label='LP01: analytical',
    color='black',
    line_width=2,
    edge_color='blue',
    marker_size=80,
    line_style='-',
    layer_position=2
)

_ = figure.show()

# -
