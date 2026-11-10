"""
Validate a taper against its radius formula
===========================================

For one unsymmetrised alpha segment, compare the assembled sampled radius
with the analytical segment expression. The numerical profile interpolates
segment samples; a small discrepancy is expected away from those samples.
This validates the taper geometry, not the optical mode solver.
"""

import matplotlib.pyplot as plt
import numpy as np
from SuPyMode.workflow import AlphaProfile

alpha, heating, stretch = 0.2, 8e-3, 4e-3
profile = AlphaProfile(
    initial_radius=1.0, n_point=401, symmetric=False, add_end_of_taper_section=False
)
profile.add_taper_segment(
    alpha=alpha, initial_heating_length=heating, stretching_length=stretch, n_point=501
)
profile.initialize()
z = profile.distance
reference = (1 + 2 * alpha * z / ((1 - alpha) * heating)) ** (-1 / (2 * alpha))
error = np.abs(profile.radius - reference)
np.testing.assert_allclose(profile.radius, reference, rtol=1e-6, atol=1e-8)
print(f"Maximum absolute radius error: {error.max():.3g}")
fig, axes = plt.subplots(2, 1, sharex=True)
axes[0].plot(z * 1e3, reference, label="Analytical")
axes[0].plot(z * 1e3, profile.radius, "--", label="Interpolated")
axes[0].set_ylabel("Normalised radius")
axes[0].legend()
axes[1].plot(z * 1e3, error)
axes[1].set(xlabel="Distance (mm)", ylabel="Absolute error")
fig.tight_layout()
plt.show()
