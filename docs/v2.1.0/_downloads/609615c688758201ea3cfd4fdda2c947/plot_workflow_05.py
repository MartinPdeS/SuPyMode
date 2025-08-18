"""
5x5 Coupler
===========
"""


# %%
# Importing the script dependencies
from SuPyMode.workflow import Workflow, fiber_loader, Boundaries, DomainAlignment, Profile, StructureType, BoundaryValue

from PyOptik import MaterialBank

# %%
# Creating the fiber list for mesh
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# In this example we want to simulate a single fiber at wavelength 1550 nm.
wavelength = 1550e-9


# %%
# Generating the fiber structure
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Here we define the cladding and fiber structure to model the problem
clad_refractive_index = MaterialBank.fused_silica.compute_refractive_index(wavelength)  # Refractive index of silica at the specified wavelength

clad_structure = Profile()

clad_structure.add_structure(
    structure_type=StructureType.CIRCULAR,
    number_of_fibers=7,
    fusion_degree=0.3,
    fiber_radius=62.5e-6,
    compute_fusing=True
)

clad_structure.refractive_index = clad_refractive_index

fiber_list = [
    fiber_loader.load_fiber('DCF1300S_33', clad_refractive_index=clad_refractive_index, position=clad_structure.cores[0]),
    fiber_loader.load_fiber('SMF28', clad_refractive_index=clad_refractive_index, position=clad_structure.cores[1]),
    fiber_loader.load_fiber('DCF1300S_33', clad_refractive_index=clad_refractive_index, position=clad_structure.cores[2]),
    fiber_loader.load_fiber('SMF28', clad_refractive_index=clad_refractive_index, position=clad_structure.cores[3]),
    fiber_loader.load_fiber('DCF1300S_33', clad_refractive_index=clad_refractive_index, position=clad_structure.cores[4]),
    fiber_loader.load_fiber('DCF1300S_33', clad_refractive_index=clad_refractive_index, position=clad_structure.cores[5]),
    fiber_loader.load_fiber('SMF28', clad_refractive_index=clad_refractive_index, position=clad_structure.cores[6]),
]


# %%
# Defining the boundaries of the system
boundaries = [
    Boundaries(right=BoundaryValue.SYMMETRIC)
]

# %%
# Generating the computing workflow
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Workflow class to define all the computation parameters before initializing the solver
workflow = Workflow(
    fiber_list=fiber_list,          # List of fiber to be added in the mesh, the order matters.
    clad_structure=clad_structure,  # Cladding structure, if None provided then no cladding is set.
    wavelength=wavelength,          # Wavelength used for the mode computation.
    resolution=80,                  # Number of point in the x and y axis [is divided by half if symmetric or anti-symmetric boundaries].
    x_bounds=DomainAlignment.LEFT,  # Mesh x-boundary structure.
    y_bounds=DomainAlignment.CENTERING, # Mesh y-boundary structure.
    boundaries=boundaries,          # Set of symmetries to be evaluated, each symmetry add a round of simulation
    n_sorted_mode=4,                # Total computed and sorted mode.
    n_added_mode=3,                 # Additional computed mode that are not considered later except for field comparison [the higher the better but the slower].
    plot_geometry=True,             # Plot the geometry mesh before computation.
    debug_mode=0,                   # Print the iteration step for the solver plus some other important steps.
    auto_label=False,               # Auto labeling the mode. Label are not always correct and should be verified afterwards.
    itr_final=0.1,                  # Final value of inverse taper ratio to simulate
    index_scrambling=1e-4           # Scrambling of refractive index value in order to lift mode degeneracy [useful for some analysis]
)

superset = workflow.get_superset()

# %%
# Field computation: :math:`E_{i,j}`
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
_ = superset.plot(plot_type='field', itr_list=[1.0, 0.1])

# %%
# Effective index: :math:`n^{eff}_{i,j}`
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
_ = superset.plot(plot_type='index')

# %%
# Modal normalized coupling: :math:`C_{i,j}`
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
_ = superset.plot(plot_type='normalized-coupling')

# %%
# Adiabatic criterion: :math:`\tilde{C}_{i,j}`
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
_ = superset.plot(plot_type='adiabatic')

# -
