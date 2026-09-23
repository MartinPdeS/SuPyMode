"""
4x4 Coupler
===========
"""

# %%
# Importing the script dependencies
from SuPyMode.workflow import Workflow, fiber_loader, Boundaries, DomainAlignment, Profile, StructureType, BoundaryValue

from PyOptik import MaterialBank

wavelength = 1550e-9


# %%
# Generating the fiber structure
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Here we define the cladding and fiber structure to model the problem
clad_refractive_index = MaterialBank.fused_silica.compute_refractive_index(wavelength)  # Refractive index of silica at the specified wavelength

clad_structure = Profile()

clad_structure.add_structure(
    structure_type=StructureType.CIRCULAR,
    number_of_fibers=4,
    fusion_degree=0.3,
    fiber_radius=62.5e-6
)

clad_structure.refractive_index = clad_refractive_index

fiber_list = [
    fiber_loader.load_fiber('DCF1300S_33', clad_refractive_index=clad_refractive_index, position=clad_structure.cores[0]),
    fiber_loader.load_fiber('DCF1300S_33', clad_refractive_index=clad_refractive_index, position=clad_structure.cores[1]),
    fiber_loader.load_fiber('DCF1300S_26', clad_refractive_index=clad_refractive_index, position=clad_structure.cores[2]),
    fiber_loader.load_fiber('DCF1300S_33', clad_refractive_index=clad_refractive_index, position=clad_structure.cores[3]),
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
    resolution=60,                  # Number of point in the x and y axis [is divided by half if symmetric or anti-symmetric boundaries].
    x_bounds=DomainAlignment.LEFT, # Mesh x-boundary structure.
    y_bounds=DomainAlignment.CENTERING, # Mesh y-boundary structure.
    boundaries=boundaries,          # Set of symmetries to be evaluated, each symmetry add a round of simulation
    n_sorted_mode=4,                # Total computed and sorted mode.
    n_added_mode=3,                 # Additional computed mode that are not considered later except for field comparison [the higher the better but the slower].
    debug_mode=0,                   # Print the iteration step for the solver plus some other important steps.
    auto_label=True,                # Auto labeling the mode. Label are not always correct and should be verified afterwards.
    itr_final=0.05,                 # Final value of inverse taper ratio to simulate
)

workflow.initialize_geometry(plot=True)  # Initialize the geometry and plot it

workflow.run_solver()  # Run the solver to compute the modes


# %%
# Field computation: :math:`E_{i,j}`
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
_ = workflow.plot(plot_type='field', itr_list=[1.0, 0.05])

# %%
# Effective index: :math:`n^{eff}_{i,j}`
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
_ = workflow.plot(plot_type='index')

# %%
# Modal normalized coupling: :math:`C_{i,j}`
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
_ = workflow.plot(plot_type='normalized-coupling')

# %%
# Adiabatic criterion: :math:`\tilde{C}_{i,j}`
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
_ = workflow.plot(plot_type='adiabatic')

# -
