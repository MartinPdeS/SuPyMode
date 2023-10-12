"""
Validation: 4 for circular symmetric structure
==============================================
"""

# %%
# Imports
# ~~~~~~~
import numpy
from SuPyMode.tools.fibermodes_validation import FiberModeSolver
from SuPyMode.workflow import Workflow, fiber_catalogue, Boundaries2D, configuration
from MPSPlots.render2D import SceneList

wavelength = 1550e-9
mode_numbers = ['LP01', 'LP02', 'LP03', 'LP41_a']
fiber_list = [fiber_catalogue.SMF28(wavelength=wavelength)]

clad_structure = configuration.ring.FusedProfile_01x01

# %%
# Generating the computing workflow
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Workflow class to define all the computation parameters before initializing the solver
workflow = Workflow(
    fiber_list=fiber_list,          # List of fiber to be added in the mesh, the order matters.
    clad_structure=clad_structure,  # Cladding structure, if None provided then no cladding is set.
    fusion_degree=None,             # Degree of fusion of the structure if applicable.
    wavelength=wavelength,          # Wavelength used for the mode computation.
    resolution=60,                  # Number of point in the x and y axis [is divided by half if symmetric or anti-symmetric boundaries].
    x_bounds="centering-left",      # Mesh x-boundary structure.
    y_bounds="centering-top",       # Mesh y-boundary structure.
    boundaries=[                    # Set of symmetries to be evaluated, each symmetry add a round of simulation
        Boundaries2D(right='symmetric', bottom='symmetric'),
        Boundaries2D(right='symmetric', bottom='anti-symmetric')
    ],
    n_sorted_mode=6,                # Total computed and sorted mode.
    n_added_mode=4,                 # Additional computed mode that are not considered later except for field comparison [the higher the better but the slower].
    plot_geometry=True,             # Plot the geometry mesh before computation.
    debug_mode=False,               # Print the iteration step for the solver plus some other important steps.
    auto_label=True,                # Auto labeling the mode. Label are not always correct and should be verified afterwards.
    itr_final=0.05,                 # Final value of inverse taper ratio to simulate
    index_scrambling=0              # Scrambling of refractive index value in order to lift mode degeneracy [useful for some analysis]
)

superset = workflow.get_superset()

# %%
fibermode_solver = FiberModeSolver(wavelength=wavelength)

fibermodes_data_sets = fibermode_solver.get_beta_vs_itr(
    mode_numbers=[m[:4] for m in mode_numbers],
    itr_list=numpy.linspace(1.0, 0.1, 100),
    debug_mode=False
)

figure = SceneList(unit_size=(12, 4))

ax = figure.append_ax(
    x_label='Inverse taper ratio',
    y_label='Effective index',
    show_legend=True,
    font_size=18,
    tick_size=15,
    legend_font_size=18
)

for idx, (mode_number, data_set) in enumerate(zip(mode_numbers, fibermodes_data_sets)):
    color = f'C{idx}'

    not_nan_idx = numpy.where(~numpy.isnan(data_set.y))
    y_data = data_set.y[not_nan_idx]
    x_data = data_set.x[not_nan_idx]

    ax.add_line(
        x=x_data,
        y=y_data,
        label=mode_number,
        line_style='-',
        line_width=2,
        color=color,
        layer_position=1
    )

    sub_samnpling = 20
    ax.add_scatter(
        x=superset.itr_list[::sub_samnpling],
        y=getattr(superset, mode_number).beta.get_values()[::sub_samnpling],
        label=mode_number,
        color='black',
        line_width=2,
        edge_color=color,
        marker_size=80,
        line_style='-',
        layer_position=2
    )

_ = figure.show()

# -
