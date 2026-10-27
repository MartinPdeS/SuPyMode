"""
Compare heating lengths
=======================

For the same stretch, a longer hot zone produces a gentler radius reduction.
Distances and heating lengths are in metres; the radius is normalised to its
initial value. This example compares profiles, not optical mode solutions.
"""

import matplotlib.pyplot as plt
from SuPyMode.workflow import AlphaProfile

fig, ax = plt.subplots()
for heating_length in (4e-3, 8e-3, 12e-3):
    profile = AlphaProfile(
        initial_radius=1.0, symmetric=False, add_end_of_taper_section=False
    )
    profile.add_taper_segment(
        alpha=0.2, initial_heating_length=heating_length, stretching_length=4e-3
    )
    profile.initialize()
    ax.plot(
        profile.distance * 1e3, profile.itr_list, label=f"{heating_length * 1e3:g} mm"
    )
ax.set(
    xlabel="Distance (mm)",
    ylabel="Inverse taper ratio",
    title="Effect of heating length",
)
ax.legend(title="Heating length")
plt.show()
