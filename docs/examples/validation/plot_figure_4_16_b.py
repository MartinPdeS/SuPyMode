"""
Adiabatic criterion: X. D. ~ B
==============================
"""

# %%
# Imports
# ~~~~~~~
from numpy import genfromtxt
from MPSPlots.render2D import SceneList
from SuPyMode.tools.directories import validation_data_path
from SuPyMode.workflow import configuration, Workflow, fiber_catalogue, Boundaries2D

# %%
# Generating the fiber structure
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Here we make use of the FiberFusing to generate a fiber structure that we use
# as the cladding.
wavelength = 1550e-9

# %%
# Creating the geometry rasterization
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# The cladding being defined we create the cores that are distributed at
# each (virtual) center of the cladding.
# All the components: air -- cladding -- cores, are inputed into a geometry class
# that will generate the mesh which will be used in the finit-difference matrix.
fiber_0 = fiber_catalogue.CustomFiber(wavelength=1550e-9)
fiber_1 = fiber_catalogue.CustomFiber(wavelength=1550e-9)

fiber_0.add_silica_pure_cladding()
fiber_1.add_silica_pure_cladding()

fiber_0.create_and_add_new_structure(name='inner-clad', radius=20.0e-6 / 2, NA=0.05)  # V_2 = 4.0
fiber_1.create_and_add_new_structure(name='inner-clad', radius=20.0e-6 / 2, NA=0.05)  # V_2 = 4.0

fiber_0.create_and_add_new_structure(name='core', radius=3.8e-6 / 2, NA=0.113)  # V_2 = 4.0
fiber_1.create_and_add_new_structure(name='core', radius=3.8e-6 / 2, NA=0.113)  # V_2 = 4.0

# %%
# Generating the fiber structure
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Here we define the cladding and fiber structure to model the problem
clad_structure = configuration.ring.FusedProfile_02x02


fiber_list = [fiber_0, fiber_1]

# %%
# Defining the boundaries of the system
boundaries = [
    Boundaries2D(right='symmetric', top='symmetric'),
    Boundaries2D(right='symmetric', top='anti-symmetric')
]

workflow = Workflow(
    fiber_list=fiber_list,          # List of fiber to be added in the mesh, the order matters.
    clad_structure=clad_structure,  # Cladding structure, if None provided then no cladding is set.
    fusion_degree=0.9,              # Degree of fusion of the structure if applicable.
    wavelength=wavelength,          # Wavelength used for the mode computation.
    resolution=40,                  # Number of point in the x and y axis [is divided by half if symmetric or anti-symmetric boundaries].
    x_bounds="centering-left",      # Mesh x-boundary structure.
    y_bounds="centering-bottom",    # Mesh y-boundary structure.
    boundaries=boundaries,          # Set of symmetries to be evaluated, each symmetry add a round of simulation.
    n_sorted_mode=4,                # Total computed and sorted mode.
    n_added_mode=2,                 # Additional computed mode that are not considered later except for field comparison [the higher the better but the slower].
    debug_mode=False,               # Print the iteration step for the solver plus some other important steps.
    auto_label=True                 # Auto labeling the mode. Label are not always correct and should be verified afterwards.
)

# %%
# Field plotting
workflow.superset.plot(plot_type='field').show()


# %%
# Preparing the figure
figure = SceneList(unit_size=(8, 5), title='SBB figure 4.16-(b)')

ax = figure.append_ax(
    y_scale='log',
    x_label='ITR',
    y_label='Adiabatic criterion',
    show_legend=True,
    y_limits=[1e-4, 1]
)

coupling_list = [
    (workflow.superset.LP01, workflow.superset.LP02, '-', 'blue'),
    (workflow.superset.LP01, workflow.superset.LP21_a, '-', 'red'),
    (workflow.superset.LP01, workflow.superset.LP41_a, '-', 'orange'),
    (workflow.superset.LP11_b, workflow.superset.LP31_b, '-', 'purple'),
    (workflow.superset.LP11_b, workflow.superset.LP12_b, '-', 'green'),
    (workflow.superset.LP11_b, workflow.superset.LP51_b, '-', 'turquoise')
]

for mode_0, mode_1, line_style, color in coupling_list:
    adiabatic = mode_0.adiabatic.get_values(mode_1)

    ax.add_line(
        x=mode_0.itr_list,
        y=adiabatic * 1e-6,
        label=f'{mode_0.stylized_label}-{mode_1.stylized_label}',
        line_style=line_style,
        color=color
    )

# %%
# Comparisons with the other datas.
data_directory = [
    (f"{validation_data_path}/SBB_figure_4_16_b/LP01-LP02.csv", 'blue'),
    (f"{validation_data_path}/SBB_figure_4_16_b/LP01-LP21.csv", 'red'),
    (f"{validation_data_path}/SBB_figure_4_16_b/LP01-LP41.csv", 'orange'),
    (f"{validation_data_path}/SBB_figure_4_16_b/LP11-LP31.csv", 'purple'),
    (f"{validation_data_path}/SBB_figure_4_16_b/LP11-LP12.csv", 'green'),
    (f"{validation_data_path}/SBB_figure_4_16_b/LP11-LP51.csv", 'turquoise')
]

for data_dir, color in data_directory:
    data = genfromtxt(data_dir, delimiter=',').T
    ax.add_line(
        x=data[0],
        y=data[1],
        line_style="--",
        color=color,
        line_width=2
    )

figure.show()

# -
