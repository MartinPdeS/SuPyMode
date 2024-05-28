#!/usr/bin/env python
# -*- coding: utf-8 -*-

import matplotlib.pyplot as plt

from MPSPlots import CMAP
from matplotlib.animation import FuncAnimation, PillowWriter


def get_3_figures():

    fig = plt.figure(figsize=(8, 6))

    gs = fig.add_gridspec(
        2, 2,
        width_ratios=(1, 1),
        height_ratios=(1, 1),
        left=0.1,
        right=0.9,
        bottom=0.1,
        top=0.9,
        wspace=0.5 * 2,
        hspace=0.15)

    ax0 = fig.add_subplot(gs[0, 0:])
    ax1 = fig.add_subplot(gs[1, 0])
    ax2 = fig.add_subplot(gs[1, 1])

    ax1.tick_params(
        axis='both',
        which='both',
        top=False,
        bottom=False,
        right=False,
        left=False,
        labelleft=False,
        labelbottom=False,
        grid_alpha=0
    )
    ax2.tick_params(
        axis='both',
        which='both',
        top=False,
        bottom=False,
        right=False,
        left=False,
        labelleft=False,
        labelbottom=False,
        grid_alpha=0
    )

    ax1.set_aspect('equal')
    ax2.set_aspect('equal')

    return fig, ax0, ax1, ax2


class ModePropagationGifCreator():
    def __init__(
            self,
            superset,
            profile,
            max_number_of_mode: int = None,
            dark_background: bool = True):

        self.superset = superset
        self.profile = profile
        self.dark_background = dark_background

        if max_number_of_mode is None:
            self.number_of_mode = len(superset.supermodes)
        else:
            self.number_of_mode = max_number_of_mode

        self.generate_figure()

    def generate_profile_ax(self):
        self.ax_profile = self.figure.add_subplot(self.grid_spec[0, :])
        self.ax_profile.set_xlabel('Propagation distance z')
        self.ax_profile.set_ylabel('Coupler profile')

        top = self.profile.radius
        bottom = -self.profile.radius

        self.ax_profile.plot(self.profile.distance, top, color='black')
        self.ax_profile.plot(self.profile.distance, bottom, color='black')
        self.ax_profile.fill_between(self.profile.distance, top, bottom, color='lightblue', alpha=0.8)

    def generate_field_ax(self):
        self.field_axes = []
        for mode in range(self.number_of_mode):
            ax = self.figure.add_subplot(self.grid_spec[1, mode])
            ax.set_aspect('equal')
            ax.set_xticks([])
            ax.set_yticks([])
            self.field_axes.append(ax)
            ax.set_title(self.superset[mode].stylized_name)

    def generate_figure(self, unit_size: tuple = (3, 6)) -> None:
        figure_size = (unit_size[0] * self.number_of_mode, unit_size[1])
        self.figure = plt.figure(figsize=figure_size)

        self.grid_spec = self.figure.add_gridspec(
            2,
            self.number_of_mode,
            left=0.1,
            right=0.95,
            bottom=0.1,
            top=0.9,
            wspace=0.5,
            hspace=0.35
        )

        self.generate_profile_ax()

        self.generate_field_ax()

    def populate_axes(self, z: float):
        itr = self.profile.master_interpolation_z_to_itr(z)
        slice_structure = self.superset.get_slice_structure(itr=itr, add_symmetries=True)

        self.ax_profile.set_title(f'Z-distance: {z:>5.3e}    ITR: {itr:>5.3f}')
        self.profile_line = self.ax_profile.axvline(z, linestyle='--', color='red')

        for ax, field in zip(self.field_axes, slice_structure.fields):
            ax.pcolormesh(field, cmap=CMAP.BKR)

    def update_axes(self, z: float):
        itr = self.profile.master_interpolation_z_to_itr(z)
        slice_structure = self.superset.get_slice_structure(itr=itr, add_symmetries=True)
        self.profile_line.set_xdata(z)
        self.ax_profile.set_title(f'Z-distance: {z:>5.3e}    ITR: {itr:>5.3f}')

        for ax, field in zip(self.field_axes, slice_structure.fields):
            ax.clear()
            ax.pcolormesh(field, cmap=CMAP.BKR)

    def make_animation(self, n_step: int = 20, dpi: float = 100, fps: int = 50):
        self.populate_axes(z=0e-3)

        factor = self.profile.length / n_step

        def animate(iteration):
            distance = iteration * factor
            print(f"iteration: {iteration} / {n_step}")

            self.update_axes(distance)

        ani = FuncAnimation(
            self.figure,
            animate,
            interval=40,
            blit=True,
            repeat=True,
            frames=n_step
        )

        ani.save(
            "propagation.gif",
            dpi=dpi,
            writer=PillowWriter(fps=fps)
        )
# -
