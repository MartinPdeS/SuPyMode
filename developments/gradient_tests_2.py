from SuPyMode.workflow import Workflow, configuration, fiber_catalogue, Boundaries


wavelength = 1550e-9

fiber_list = [
    fiber_catalogue.load_fiber('DCF1300S_20', wavelength=wavelength),
    fiber_catalogue.load_fiber('DCF1300S_26', wavelength=wavelength),
    fiber_catalogue.load_fiber('DCF1300S_33', wavelength=wavelength)
]

clad_profile = configuration.ring.FusedProfile_03x03

capillary_tube = fiber_catalogue.CapillaryTube(
    radius=150e-6,
    wavelength=wavelength,
    delta_n=-15e-3
)


workflow = Workflow(
    fiber_list=fiber_list,          # List of fiber to be added in the mesh, the order matters.
    clad_structure=clad_profile,    # Cladding structure, if None provided then no cladding is set.
    capillary_tube=capillary_tube,  # In this case we add a wrapping capillary tube around the fused structure.
    fusion_degree='auto',           # Degree of fusion of the structure if applicable.
    wavelength=wavelength,          # Wavelength used for the mode computation.
    resolution=100,                 # Number of point in the x and y axis [is divided by half if symmetric or anti-symmetric boundaries].
    x_bounds="centering",           # Mesh x-boundary structure.
    y_bounds="centering",           # Mesh y-boundary structure.
    boundaries=[Boundaries()],      # Set of symmetries to be evaluated, each symmetry add a round of simulation
    n_sorted_mode=6,                # Total computed and sorted mode.
    n_added_mode=6,                 # Additional computed mode that are not considered later except for field comparison [the higher the better but the slower].
    plot_geometry=True,             # Plot the geometry mesh before computation.
    debug_mode=1,                   # Define the degree of debug printout, from 0 to 3. [Does not work properly on jupyter notebooks]
    auto_label=True,                # Auto labeling the mode. Label are not always correct and should be verified afterwards.
    itr_final=0.1,                  # Final value of inverse taper ratio to simulate
    index_scrambling=0e-4           # Scrambling of refractive index value in order to lift mode degeneracy [useful for some analysis]
)

superset = workflow.get_superset()