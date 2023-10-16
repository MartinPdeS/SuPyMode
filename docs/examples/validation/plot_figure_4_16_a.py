"""
Validation: 1 with Xavier Daxhelet program
==========================================
"""


# %%
# Imports
# ~~~~~~~
from numpy import genfromtxt
from MPSPlots.render2D import SceneList
from SuPyMode.tools.directories import validation_data_path
from SuPyMode.workflow import configuration, Workflow, fiber_catalogue, Boundaries2D

# %%
# Generating the fiber structure
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Here we make use of the FiberFusing to generate a fiber structure that we use
# as the cladding.

clad_structure = configuration.ring.FusedProfile_02x02

workflow = Workflow(
    fiber_list=[fiber_catalogue.SMF28(wavelength=1550e-9), fiber_catalogue.SMF28(wavelength=1550e-9)],
    clad_structure=clad_structure,
    fusion_degree=0.9,
    wavelength=1550e-9,
    resolution=60,
    x_bounds="centering-left",
    y_bounds="centering-bottom",
    boundaries=[
        Boundaries2D(right='symmetric', top='symmetric'),
        Boundaries2D(right='symmetric', top='anti-symmetric')
    ],
    n_sorted_mode=4,
    n_added_mode=2,
    debug_mode=False,
    auto_label=True
)

# %%
# Field plotting
workflow.superset.plot(plot_type='field').show()

workflow.superset.save_instance("figure4_16")

figure = SceneList(unit_size=(12, 5), title='SBB figure 4.16')

ax = figure.append_ax(
    y_scale='log',
    x_label='ITR',
    y_label='Adiabatic criterion',
    show_legend=True,
    y_limits=[1e-4, 1]
)

coupling_list = [
    (workflow.superset.LP01, workflow.superset.LP02, '-', 'blue'),
    (workflow.superset.LP01, workflow.superset.LP21_a, '-', 'red'),
    (workflow.superset.LP01, workflow.superset.LP41_a, '-', 'orange'),
    (workflow.superset.LP11_b, workflow.superset.LP31_b, '-', 'purple'),
    (workflow.superset.LP11_b, workflow.superset.LP12_b, '-', 'green'),
    (workflow.superset.LP11_b, workflow.superset.LP51_b, '-', 'turquoise')
]

for mode_0, mode_1, line_style, color in coupling_list:
    adiabatic = mode_0.adiabatic.get_values(mode_1)

    ax.add_line(
        x=mode_0.itr_list,
        y=adiabatic * 1e-6,
        label=f'{mode_0.stylized_label}-{mode_1.stylized_label}',
        line_style=line_style,
        color=color
    )

# %%
# Comparisons with the other datas.
data_directory = [
    (f"{validation_data_path}/SBB_figure_4_16_a/LP01-LP02.csv", 'blue'),
    (f"{validation_data_path}/SBB_figure_4_16_a/LP01-LP21.csv", 'red'),
    (f"{validation_data_path}/SBB_figure_4_16_a/LP01-LP41.csv", 'orange'),
    (f"{validation_data_path}/SBB_figure_4_16_a/LP11-LP31.csv", 'purple'),
    (f"{validation_data_path}/SBB_figure_4_16_a/LP11-LP12.csv", 'green'),
    (f"{validation_data_path}/SBB_figure_4_16_a/LP11-LP51.csv", 'turquoise')
]

for data_dir, color in data_directory:
    data = genfromtxt(data_dir, delimiter=',').T

    ax.add_line(
        x=data[0],
        y=data[1],
        line_style="--",
        color=color,
        label=None,
        line_width=2
    )

figure.show()


# -
