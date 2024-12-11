from SuPyMode.workflow import AlphaProfile
import matplotlib.pyplot as plt

figure, ax = plt.subplots(1, 1)

for n, (l, t) in enumerate([(2e-3, 100), (4e-3, 100), (8e-3, 100)]):

    profile = AlphaProfile(label=f'Initial length: {l*1e3:.0f}mm', symmetric=False, add_end_of_taper_section=False, line_color=f'C{n}')

    profile.add_taper_segment(
        alpha=0,
        initial_heating_length=l,
        stretching_length=0.2e-3 * t
    )

    profile.initialize()

    profile.plot(show_radius=True, show_adiabatic=False, show_taper_angle=False, show=False, ax=ax)


plt.show()