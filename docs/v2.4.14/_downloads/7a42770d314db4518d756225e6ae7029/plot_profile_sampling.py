"""
Inspect taper sampling
======================

Compare two sampling densities of the same taper. A coarse grid is useful for
exploration; increase the sampling and check convergence before using a profile
in a propagation calculation.
"""

import matplotlib.pyplot as plt
from SuPyMode.workflow import AlphaProfile

fig, ax = plt.subplots()
for points in (20, 200):
    profile = AlphaProfile(
        n_point=points, symmetric=False, add_end_of_taper_section=False
    )
    profile.add_taper_segment(
        alpha=0.2, initial_heating_length=8e-3, stretching_length=4e-3
    )
    profile.initialize()
    ax.plot(
        profile.distance * 1e3,
        profile.itr_list,
        marker="." if points == 20 else None,
        label=f"{points} points",
    )
ax.set(xlabel="Distance (mm)", ylabel="Inverse taper ratio")
ax.legend()
plt.show()
