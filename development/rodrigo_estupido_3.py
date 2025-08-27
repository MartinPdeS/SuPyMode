from SuPyMode.workflow import Workflow, configuration, fiber_catalogue, Boundaries, AlphaProfile

wavelength = 1550e-9


fiber_0 = fiber_catalogue.GenericFiber(wavelength=wavelength)
fiber_0.add_silica_pure_cladding()
fiber_0.create_and_add_new_structure(name='core', radius=9.0e-6 / 2, NA=0.14)

fiber_1 = fiber_catalogue.GenericFiber(wavelength=wavelength)
fiber_1.add_silica_pure_cladding()
fiber_1.create_and_add_new_structure(name='core', radius=11.0e-6 / 2, NA=0.14)

fiber_2 = fiber_catalogue.GenericFiber(wavelength=wavelength)
fiber_2.add_silica_pure_cladding()
fiber_2.create_and_add_new_structure(name='core', radius=13.0e-6 / 2, NA=0.14)

fiber_configuration_0 = [
    fiber_0, fiber_1, fiber_2
]


# fiber_0 = fiber_catalogue.GenericFiber(wavelength=wavelength)
# fiber_0.add_silica_pure_cladding()
# fiber_0.create_and_add_new_structure(name='core', radius=30.0e-6 / 2, NA=0.19)  # or 0.15
# fiber_0.create_and_add_new_structure(name='core', radius=9.0e-6 / 2, NA=0.14)

# fiber_1 = fiber_catalogue.GenericFiber(wavelength=wavelength)
# fiber_1.add_silica_pure_cladding()
# fiber_1.create_and_add_new_structure(name='core', radius=30.0e-6 / 2, NA=0.19)  # or 0.15
# fiber_1.create_and_add_new_structure(name='core', radius=11.0e-6 / 2, NA=0.14)

# fiber_2 = fiber_catalogue.GenericFiber(wavelength=wavelength)
# fiber_2.add_silica_pure_cladding()
# fiber_2.create_and_add_new_structure(name='core', radius=30.0e-6 / 2, NA=0.19)  # or 0.15
# fiber_2.create_and_add_new_structure(name='core', radius=13.0e-6 / 2, NA=0.14)

# fiber_configuration_1 = [
#     fiber_0, fiber_1, fiber_2
# ]


# fiber_0 = fiber_catalogue.GenericFiber(wavelength=wavelength)
# fiber_0.add_silica_pure_cladding()
# fiber_0.create_and_add_new_structure(name='core', radius=30.0e-6 / 2, NA=0.19)
# fiber_0.create_and_add_new_structure(name='core', radius=9.0e-6 / 2, NA=0.14)

# fiber_1 = fiber_catalogue.GenericFiber(wavelength=wavelength)
# fiber_1.add_silica_pure_cladding()
# fiber_1.create_and_add_new_structure(name='core', radius=30.0e-6 / 2, NA=0.19)
# fiber_1.create_and_add_new_structure(name='core', radius=11.0e-6 / 2, NA=0.14)

# fiber_2 = fiber_catalogue.GenericFiber(wavelength=wavelength)
# fiber_2.add_silica_pure_cladding()
# fiber_2.create_and_add_new_structure(name='core', radius=30.0e-6 / 2, NA=0.19)
# fiber_2.create_and_add_new_structure(name='core', radius=13.0e-6 / 2, NA=0.14)

# fiber_configuration_1 = [
#     fiber_catalogue.load_fiber('DCF1300S_42', wavelength=wavelength),
#     fiber_catalogue.load_fiber('DCF1300S_20', wavelength=wavelength),
#     fiber_catalogue.load_fiber('DCF1300S_26', wavelength=wavelength),
# ]


clad_profile = configuration.ring.FusedProfile_03x03

capillary_tube = fiber_catalogue.CapillaryTube(
    radius=150e-6,
    wavelength=wavelength,
    delta_n=-15e-3
)


workflow = Workflow(
    fiber_list=fiber_configuration_0,  # List of fiber to be added in the mesh, the order matters.
    clad_structure=clad_profile,    # Cladding structure, if None provided then no cladding is set.
    capillary_tube=capillary_tube,  # In this case we add a wrapping capillary tube around the fused structure.
    fusion_degree='auto',           # Degree of fusion of the structure if applicable.
    wavelength=wavelength,          # Wavelength used for the mode computation.
    resolution=140,                 # Number of point in the x and y axis [is divided by half if symmetric or anti-symmetric boundaries].
    x_bounds="centering",           # Mesh x-boundary structure.
    y_bounds="centering",           # Mesh y-boundary structure.
    boundaries=[Boundaries()],      # Set of symmetries to be evaluated, each symmetry add a round of simulation
    n_sorted_mode=4,                # Total computed and sorted mode.
    n_added_mode=3,                 # Additional computed mode that are not considered later except for field comparison [the higher the better but the slower].
    # plot_geometry=True,             # Plot the geometry mesh before computation.
    plot_field=True,
    debug_mode=1,                   # Define the degree of debug printout, from 0 to 3. [Does not work properly on jupyter notebooks]
    auto_label=True,                # Auto labeling the mode. Label are not always correct and should be verified afterwards.
    itr_final=0.05,                  # Final value of inverse taper ratio to simulate
    itr_initial=1.00,
    n_step=500,
    index_scrambling=0e-4           # Scrambling of refractive index value in order to lift mode degeneracy [useful for some analysis]
)

superset = workflow.get_superset()


profile = AlphaProfile(symmetric=False, add_end_of_taper_section=True)

profile.add_taper_segment(
    alpha=0,
    initial_heating_length=2e-3,
    stretching_length=0.2e-3 * 20
)

profile.initialize()

superset.plot_adiabatic(add_profile=profile)
