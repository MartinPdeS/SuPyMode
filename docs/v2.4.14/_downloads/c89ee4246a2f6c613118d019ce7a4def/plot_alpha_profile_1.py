"""
Symmetric coupler z-profile
===========================
"""

# %%
# Importing the script dependencies
from SuPyMode.workflow import AlphaProfile
import matplotlib.pyplot as plt

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


figure, axes = plt.subplots(nrows=3, ncols=1, figsize=(10, 9))

axes[0].plot(profile.itr_list, profile.adiabatic)
axes[0].set_xlabel("Inverse taper ratio")
axes[0].set_ylabel("Adiabatic criterion")

axes[1].plot(profile.distance, profile.radius)
axes[1].set_xlabel("Z-propagation")
axes[1].set_ylabel("Radius")

axes[2].plot(profile.distance, profile.taper_angle)
axes[2].set_xlabel("Z-propagation")
axes[2].set_ylabel("Taper angle")
plt.tight_layout()
plt.show()

# -
