"""
3x3 Coupler
===========
"""

# %%
# Importing the script dependencies
from SuPyMode.workflow import Workflow, configuration, fiber_catalogue, Boundaries

# %%
# Creating the fiber list for mesh
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# In this example we want to simulate a single fiber at wavelength 1550 nm.
wavelength = 900e-9

# %%
# Generating the fiber structure
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Here we define the cladding and fiber structure to model the problem
clad_structure = configuration.ring.FusedProfile_03x03

custom_fiber_1 = fiber_catalogue.CustomFiber(wavelength=wavelength)
custom_fiber_1.add_silica_pure_cladding(radius=62.5e-6, name='outer-clad')

custom_fiber_1.create_and_add_new_structure(
    radius=11.5e-6,
    NA=0.2,
    name='inner-clad'
)

custom_fiber_1.create_and_add_new_structure(
    radius=2.0e-6,
    NA=0.13,
    name='core'
)


custom_fiber_2 = fiber_catalogue.CustomFiber(wavelength=wavelength)
custom_fiber_2.add_silica_pure_cladding(radius=62.5e-6, name='outer-clad')

custom_fiber_2.create_and_add_new_structure(
    radius=13.1e-6,
    NA=0.19,
    name='inner-clad'
)

custom_fiber_2.create_and_add_new_structure(
    radius=2.05e-6,
    NA=0.117,
    name='core'
)


custom_fiber_3 = fiber_catalogue.CustomFiber(wavelength=wavelength)
custom_fiber_3.add_silica_pure_cladding(radius=62.5e-6, name='outer-clad')

custom_fiber_3.create_and_add_new_structure(
    radius=14.9e-6,
    NA=0.2,
    name='inner-clad'
)

custom_fiber_3.create_and_add_new_structure(
    radius=2.75e-6,
    NA=0.14,
    name='core'
)

fiber_list = [
    custom_fiber_1,
    custom_fiber_2,
    custom_fiber_3
]


# %%
# Defining the boundaries of the system
boundaries = [
    Boundaries(),
]

# %%
# Generating the computing workflow
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Workflow class to define all the computation parameters before initializing the solver
workflow = Workflow(
    fiber_list=fiber_list,          # List of fiber to be added in the mesh, the order matters.
    clad_structure=clad_structure,  # Cladding structure, if None provided then no cladding is set.
    fusion_degree=0.8,           # Degree of fusion of the structure if applicable.
    wavelength=wavelength,          # Wavelength used for the mode computation.
    resolution=220,                  # Number of point in the x and y axis [is divided by half if symmetric or anti-symmetric boundaries].
    x_bounds="centering",           # Mesh x-boundary structure.
    y_bounds="centering",           # Mesh y-boundary structure.
    boundaries=boundaries,          # Set of symmetries to be evaluated, each symmetry add a round of simulation
    n_sorted_mode=5,                # Total computed and sorted mode.
    n_added_mode=5,                 # Additional computed mode that are not considered later except for field comparison [the higher the better but the slower].
    plot_geometry=True,             # Plot the geometry mesh before computation.
    air_padding_factor=1.1,
    debug_mode=2,                   # Print the iteration step for the solver plus some other important steps.
    auto_label=True,                # Auto labeling the mode. Label are not always correct and should be verified afterwards.
    itr_final=0.01,                  # Final value of inverse taper ratio to simulate
    clad_rotation=0,                # Rotate the geoemtry in the given angle in degree
    index_scrambling=1e-7           # Scrambling of refractive index value in order to lift mode degeneracy [useful for some analysis]
)

superset = workflow.get_superset()

# # %%
# # Field computation: :math:`E_{i,j}`
# # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# _ = superset.plot(plot_type='field', slice_list=[], itr_list=[1.0, 0.3, 0.01]).show()

# # %%
# # Effective index: :math:`n^{eff}_{i,j}`
# # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# _ = superset.plot(plot_type='index').show()

# # %%
# # Modal normalized coupling: :math:`C_{i,j}`
# # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# _ = superset.plot(plot_type='normalized-coupling').show()

# # %%
# # Adiabatic criterion: :math:`\tilde{C}_{i,j}`
# # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# _ = superset.plot(plot_type='adiabatic').show()

# superset.save_superset_instance()
superset.generate_pdf_report(filename='rodrigo_mames_las_vergas.pdf', directory='auto')

# -
