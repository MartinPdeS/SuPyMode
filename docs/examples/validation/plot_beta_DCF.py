"""
Propagation constant: DCFC
==========================
"""

# %%
# Imports
# ~~~~~~~
import numpy
from SuPyMode.workflow import Workflow, fiber_catalogue, Boundaries
from PyFiberModes import LP01
from PyFiberModes.fiber import load_fiber
import matplotlib.pyplot as plt

wavelength = 1550e-9
fiber_name = 'DCF1300S_33'

# %%
# Generating the fiber structure
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Here we define the cladding and fiber structure to model the problem
fiber_list = [
    fiber_catalogue.load_fiber(fiber_name, wavelength=wavelength)
]


# %%
# Defining the boundaries of the system
boundaries = [
    Boundaries(right='symmetric', bottom='symmetric'),
    Boundaries(right='symmetric', bottom='anti-symmetric')
]

# %%
# Generating the computing workflow
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Workflow class to define all the computation parameters before initializing the solver
workflow = Workflow(
    fiber_list=fiber_list,          # List of fiber to be added in the mesh, the order matters.
    fusion_degree='auto',           # Degree of fusion of the structure if applicable.
    wavelength=wavelength,          # Wavelength used for the mode computation.
    resolution=80,                  # Number of point in the x and y axis [is divided by half if symmetric or anti-symmetric boundaries].
    x_bounds="left",                # Mesh x-boundary structure.
    y_bounds="top",                 # Mesh y-boundary structure.
    air_padding_factor=4.0,
    boundaries=boundaries,          # Set of symmetries to be evaluated, each symmetry add a round of simulation
    n_sorted_mode=6,                # Total computed and sorted mode.
    n_added_mode=4,                 # Additional computed mode that are not considered later except for field comparison [the higher the better but the slower].
    plot_geometry=True,             # Plot the geometry mesh before computation.
    debug_mode=0,                   # Print the iteration step for the solver plus some other important steps.
    auto_label=True,                # Auto labeling the mode. Label are not always correct and should be verified afterwards.
    itr_final=0.2,                  # Final value of inverse taper ratio to simulate
    index_scrambling=0,             # Scrambling of refractive index value in order to lift mode degeneracy [useful for some analysis]
    n_step=50
)

superset = workflow.get_superset()
itr_list = superset.model_parameters.itr_list

# %%
# Computing the analytical values using FiberModes solver.
dcf_fiber = load_fiber(
    fiber_name=fiber_name,
    wavelength=wavelength,
    add_air_layer=True
)

# %%
# Preparing the figure
figure, ax = plt.subplots(1, 1)
ax.set(
    xlabel='Inverse taper ratio',
    ylabel='Effective index',
)

pyfibermodes_mode = LP01
supymode_mode = superset.LP01

analytical = numpy.empty(itr_list.shape)
for idx, itr in enumerate(itr_list):
    _fiber = dcf_fiber.scale(factor=itr)
    analytical[idx] = _fiber.get_effective_index(mode=pyfibermodes_mode)

ax.plot(
    itr_list,
    analytical,
    label=str(pyfibermodes_mode),
    linestyle='-',
    linewidth=2,
    color='red',
)

ax.scatter(
    itr_list,
    supymode_mode.index.data,
    label=str(supymode_mode),
    color='black',
    s=80,
    linestyle='-',
)

ax.legend()

plt.show()


# -
