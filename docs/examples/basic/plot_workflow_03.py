"""
3x3 Coupler
===========
"""

# %%
# Importing the script dependencies
from SuPyMode.workflow import Workflow, configuration, fiber_catalogue, Boundaries
from SuPyMode.workflow import AlphaProfile

# %%
# Creating the fiber list for mesh
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# In this example we want to simulate a single fiber at wavelength 1550 nm.
wavelength = 1550e-9

# %%
# Generating the fiber structure
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Here we define the cladding and fiber structure to model the problem
clad_structure = configuration.ring.FusedProfile_03x03

fiber_list = [
    fiber_catalogue.load_fiber('DCF1300S_42', wavelength=wavelength),
    fiber_catalogue.load_fiber('DCF1300S_33', wavelength=wavelength),
    fiber_catalogue.load_fiber('DCF1300S_26', wavelength=wavelength),
]

# %%
# Defining the boundaries of the system
boundaries = [
    Boundaries(),
]

# # %%
# # Generating the computing workflow
# # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# # Workflow class to define all the computation parameters before initializing the solver
# workflow = Workflow(
#     fiber_list=fiber_list,          # List of fiber to be added in the mesh, the order matters.
#     clad_structure=clad_structure,  # Cladding structure, if None provided then no cladding is set.
#     fusion_degree=0.8,           # Degree of fusion of the structure if applicable.
#     wavelength=wavelength,          # Wavelength used for the mode computation.
#     resolution=30,                  # Number of point in the x and y axis [is divided by half if symmetric or anti-symmetric boundaries].
#     x_bounds="centering",           # Mesh x-boundary structure.
#     y_bounds="centering",           # Mesh y-boundary structure.
#     boundaries=boundaries,          # Set of symmetries to be evaluated, each symmetry add a round of simulation
#     n_sorted_mode=2,                # Total computed and sorted mode.
#     n_added_mode=5,                 # Additional computed mode that are not considered later except for field comparison [the higher the better but the slower].
#     # plot_geometry=True,             # Plot the geometry mesh before computation.
#     debug_mode=3,                   # Print the iteration step for the solver plus some other important steps.
#     auto_label=True,                # Auto labeling the mode. Label are not always correct and should be verified afterwards.
#     itr_final=0.05,                  # Final value of inverse taper ratio to simulate
#     clad_rotation=0,                # Rotate the geoemtry in the given angle in degree
#     index_scrambling=1e-7           # Scrambling of refractive index value in order to lift mode degeneracy [useful for some analysis]
# )

# superset = workflow.get_superset()

# # %%
# # Field computation: :math:`E_{i,j}`
# # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# _ = superset.plot(plot_type='field', slice_list = [], itr_list=[1.0, 0.3, 0.1]).show()

# # %%
# # Effective index: :math:`n^{eff}_{i,j}`
# # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# _ = superset.plot(plot_type='index').show()

# _ = superset.plot(plot_type='beta').show()

# # %%
# # Modal normalized coupling: :math:`C_{i,j}`
# # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# _ = superset.plot(plot_type='normalized-coupling').show()

profile = AlphaProfile(symmetric=True, add_end_of_taper_section=True)

# %%
# Adding a first taper segment with large initial heating length (i.e. slow reduction)
profile.add_taper_segment(
    alpha=0,
    initial_heating_length=8e-3,
    stretching_length=0.2e-3 * 40
)

# %%
# Adding a first taper segment with small initial heating length (i.e. steep reduction)
profile.add_taper_segment(
    alpha=0,
    initial_heating_length=2e-3,
    stretching_length=0.2e-3 * 20
)


profile.initialize()


figure = profile.plot(
    show_radius=True,
    show_adiabatic=True,
    show_taper_angle=True
)

figure.show()

# # %%
# # Adiabatic criterion: :math:`\tilde{C}_{i,j}`
# # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# _ = superset.plot(plot_type='adiabatic', add_profile=profile).show()

# -
