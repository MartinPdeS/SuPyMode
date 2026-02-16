"""
1x1 Coupler
===========
"""

# %%
# Importing the script dependencies
from SuPyMode.workflow import (
    Workflow,
    fiber_loader,
    Boundaries,
    BoundaryValue,
    DomainAlignment,
)
from PyOptik import MaterialBank

wavelength = 1550e-9


# %%
# Generating the fiber structure
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Here we define the cladding and fiber structure to model the problem

clad_refractive_index = MaterialBank.fused_silica.compute_refractive_index(
    wavelength
)  # Refractive index of silica at the specified wavelength

fiber = fiber_loader.load_fiber(
    "SMF28", clad_refractive_index=clad_refractive_index, remove_cladding=False
)

# %%
# Defining the boundaries of the system
boundaries = [
    Boundaries(right=BoundaryValue.SYMMETRIC, top=BoundaryValue.SYMMETRIC),
    Boundaries(right=BoundaryValue.SYMMETRIC, top=BoundaryValue.ANTI_SYMMETRIC),
]

# %%
# Generating the computing workflow
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Workflow class to define all the computation parameters before initializing the solver
workflow = Workflow(
    fiber_list=[fiber],  # List of fiber to be added in the mesh, the order matters.
    wavelength=wavelength,  # Wavelength used for the mode computation.
    resolution=30,  # Number of point in the x and y axis [is divided by half if symmetric or anti-symmetric boundaries].
    x_bounds=DomainAlignment.LEFT,  # Mesh x-boundary structure.
    y_bounds=DomainAlignment.BOTTOM,  # Mesh y-boundary structure.
    boundaries=boundaries,  # Set of symmetries to be evaluated, each symmetry add a round of simulation
    n_sorted_mode=4,  # Total computed and sorted mode.
    n_added_mode=2,  # Additional computed mode that are not considered later except for field comparison [the higher the better but the slower].
    debug_mode=0,  # Print the iteration step for the solver plus some other important steps.
    auto_label=True,  # Auto labeling the mode. Label are not always correct and should be verified afterwards.
    itr_final=0.05,  # Final value of inverse taper ratio to simulate
    clad_rotation=0,  # Rotate the geoemtry in the given angle in degree
)

workflow.initialize_geometry(plot=False)  # Initialize the geometry and plot it

workflow.run_solver()  # Run the solver to compute the modes

_ = workflow.plot(plot_type="geometry")

# %%
# Field computation: :math:`E_{i,j}`
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
_ = workflow.plot(plot_type="field", itr_list=[1.0, 0.1])

# %%
# Effective index: :math:`n^{eff}_{i,j}`
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
_ = workflow.plot(plot_type="index")

# %%
# Modal normalized coupling: :math:`C_{i,j}`
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
_ = workflow.plot(plot_type="normalized-coupling")

# %%
# Adiabatic criterion: :math:`\tilde{C}_{i,j}`
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
_ = workflow.plot(plot_type="adiabatic")

# -
