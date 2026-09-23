"""
1x1 Coupler -- Base Workflow
============================
"""

# %%
# Importing the script dependencies
from SuPyMode.workflow import (
    Profile,
    GenericFiber,
    Boundaries,
    BoundaryValue,
    DomainAlignment,
    BackGround,
    Geometry,
    SuPySolver,
    StructureType,
)

wavelength = 1550e-9

clad_structure = Profile()

clad_structure.add_structure(
    structure_type=StructureType.CIRCULAR,
    number_of_fibers=2,
    fusion_degree=0.9,
    fiber_radius=62.5e-6,
)

clad_structure.refractive_index = 1.4444

fiber_0 = GenericFiber(position=clad_structure.cores[0])

fiber_0.create_and_add_new_structure(
    name="core", refractive_index=1.4480, radius=4.2 * 1e-6
)

fiber_1 = GenericFiber(position=clad_structure.cores[1])

fiber_1.create_and_add_new_structure(
    name="core", refractive_index=1.4460, radius=6.2 * 1e-6
)

background = BackGround(refractive_index=1.0)

geometry = Geometry(
    x_bounds=DomainAlignment.LEFT,
    y_bounds=DomainAlignment.CENTERING,
    resolution=80,
    boundary_pad_factor=1.1,
)

geometry.add_structure(background, clad_structure, fiber_0, fiber_1)

geometry.initialize()

geometry.plot()

solver = SuPySolver(
    mesh=geometry.mesh,
    x=geometry.coordinate_system.x_vector,
    y=geometry.coordinate_system.y_vector,
    tolerance=1e-20,
    max_iteration=5000,
    accuracy=2,
    debug_mode=False,
    extrapolation_order=2,
)

solver.init_superset(
    wavelength=wavelength,
    n_step=500,
    itr_initial=1.0,
    itr_final=0.05,
)

boundary = Boundaries(right=BoundaryValue.SYMMETRIC)

solver.add_modes(
    n_added_mode=3,
    n_sorted_mode=2,
    boundaries=boundary,
    auto_label=True,
)

boundary = Boundaries(right=BoundaryValue.ANTI_SYMMETRIC)

solver.add_modes(
    n_added_mode=3,
    n_sorted_mode=2,
    boundaries=boundary,
    auto_label=True,
)

superset = solver.superset

_ = superset.plot(plot_type="geometry")

# %%
# Field computation: :math:`E_{i,j}`
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
_ = superset.plot(plot_type="field", itr_list=[1.0, 0.1])

# %%
# Effective index: :math:`n^{eff}_{i,j}`
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
_ = superset.plot(plot_type="index")

# %%
# Modal normalized coupling: :math:`C_{i,j}`
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
_ = superset.plot(plot_type="normalized-coupling")

# %%
# Adiabatic criterion: :math:`\tilde{C}_{i,j}`
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
_ = superset.plot(plot_type="adiabatic")

# -
