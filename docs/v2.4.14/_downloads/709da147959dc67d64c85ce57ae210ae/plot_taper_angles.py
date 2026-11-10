"""
Inspect taper angle
===================

The taper-angle quantity is the magnitude of the sampled radius slope, ``|dr/dz|``.
Compare two heating lengths before choosing a profile for an optical calculation.
"""

import matplotlib.pyplot as plt
from SuPyMode.workflow import AlphaProfile

fig, ax = plt.subplots()
for heating in (4e-3, 10e-3):
    profile = AlphaProfile(
        initial_radius=62.5e-6,
        n_point=200,
        symmetric=False,
        add_end_of_taper_section=False,
    )
    profile.add_taper_segment(
        alpha=0.2, initial_heating_length=heating, stretching_length=4e-3
    )
    profile.initialize()
    ax.plot(profile.distance * 1e3, profile.taper_angle, label=f"{heating * 1e3:g} mm")
ax.set(xlabel="Distance (mm)", ylabel="Radius slope |dr/dz|")
ax.legend(title="Heating length")
plt.show()
