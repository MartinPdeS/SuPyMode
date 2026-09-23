from SuPyMode.workflow import AlphaProfile
from SuPyMode.plotter import Plotter

profile = AlphaProfile(symmetric=True, add_end_of_taper_section=True)


# %%
# Adding a first taper segment with large initial heating length (i.e. slow reduction)
profile.add_taper_segment(
    alpha=0, initial_heating_length=8e-3, stretching_length=0.2e-3 * 20
)


# %%
# Adding a first taper segment with small initial heating length (i.e. steep reduction)
profile.add_taper_segment(
    alpha=0, initial_heating_length=2e-3, stretching_length=0.2e-3 * 20
)

profile.initialize()

from SuPyMode.plotter import Plotter

plotter = Plotter(kind="adiabatic")
plotter.add_superset(workflow.superset, mode_of_interest="all", combination="pairs")
plotter.add_profile(profile)
plotter.show(add_profiles=True)
