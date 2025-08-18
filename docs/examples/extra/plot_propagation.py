"""
4x4 Coupler [linear]
====================
"""

# %%
# Importing the script dependencies
from SuPyMode.workflow import Workflow, fiber_loader, Boundaries, DomainAlignment, Profile, StructureType
from SuPyMode.workflow import AlphaProfile

wavelength = 1550e-9

clad_refractive_index = 1.4444


# %%
# Generating the fiber structure
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Here we define the cladding and fiber structure to model the problem
clad_structure = Profile()

clad_structure.add_structure(
    structure_type=StructureType.CIRCULAR,
    number_of_fibers=3,
    fusion_degree=0.3,
    fiber_radius=62.5e-6,
    compute_fusing=True
)

clad_structure.refractive_index = clad_refractive_index


fiber_list = [
    fiber_loader.load_fiber('DCF1300S_20', clad_refractive_index=clad_refractive_index, position=clad_structure.cores[0]),
    fiber_loader.load_fiber('DCF1300S_26', clad_refractive_index=clad_refractive_index, position=clad_structure.cores[1]),
    fiber_loader.load_fiber('DCF1300S_33', clad_refractive_index=clad_refractive_index, position=clad_structure.cores[2])
]

capillary_tube = Profile()

capillary_tube.add_center_fiber(
    fiber_radius=150e-6,
)

capillary_tube.refractive_index = clad_structure.refractive_index - 15e-3

# %%
# Generating the computing workflow
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Workflow class to define all the computation parameters before initializing the solver
workflow = Workflow(
    fiber_list=fiber_list,              # List of fiber to be added in the mesh, the order matters.
    clad_structure=clad_structure,      # Cladding structure, if None provided then no cladding is set.
    capillary_tube=capillary_tube,
    wavelength=wavelength,              # Wavelength used for the mode computation.
    resolution=30,                      # Number of point in the x and y axis [is divided by half if symmetric or anti-symmetric boundaries].
    x_bounds=DomainAlignment.CENTERING, # Mesh x-boundary structure.
    y_bounds=DomainAlignment.CENTERING, # Mesh y-boundary structure.
    boundaries=[Boundaries()],      # Set of symmetries to be evaluated, each symmetry add a round of simulation
    n_sorted_mode=3,                # Total computed and sorted mode.
    n_added_mode=6,                 # Additional computed mode that are not considered later except for field comparison [the higher the better but the slower].
    debug_mode=0,                   # Print the iteration step for the solver plus some other important steps. [Does not work properly on jupyter notebooks]
    auto_label=False,               # Auto labeling the mode. Label are not always correct and should be verified afterwards.
    itr_final=0.1,                  # Final value of inverse taper ratio to simulate
)

workflow.initialize_geometry(plot=True)  # Initialize the geometry and plot it

workflow.run_solver()  # Run the solver to compute the modes

profile = AlphaProfile(symmetric=False, add_end_of_taper_section=True)

# %%
# Adding a first taper segment with large initial heating length (i.e. slow reduction)
profile.add_taper_segment(
    alpha=0,
    initial_heating_length=0.1e-3,
    stretching_length=0.1e-3
)

profile.initialize()

propagation = workflow.superset.propagate(add_coupling=False, profile=profile, initial_amplitude=[1, 0, 0])

propagation.plot()
# -
