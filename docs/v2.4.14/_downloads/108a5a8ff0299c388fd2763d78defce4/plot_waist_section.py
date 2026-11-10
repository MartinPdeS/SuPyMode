"""
Add a constant waist
====================

A constant section after the taper models a uniform waist. The mirror option
produces an output transition with the same radius profile as the input.
"""

import matplotlib.pyplot as plt
from SuPyMode.workflow import AlphaProfile

profile = AlphaProfile(
    initial_radius=62.5e-6, symmetric=True, add_end_of_taper_section=False
)
profile.add_taper_segment(
    alpha=0.2, initial_heating_length=8e-3, stretching_length=4e-3
)
profile.add_constant_segment(length=2e-3)
profile.initialize()
fig, ax = plt.subplots()
ax.plot(profile.distance * 1e3, profile.radius * 1e6)
ax.set(xlabel="Distance (mm)", ylabel="Radius (µm)", title="Symmetric taper with waist")
plt.show()
