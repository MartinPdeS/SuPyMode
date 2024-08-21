#!/usr/bin/env python
# -*- coding: utf-8 -*-

from typing import NoReturn

import numpy
import pyvista

from SuPyMode.profiles import AlphaProfile
import matplotlib.pyplot as plt


class Propagation:
    def __init__(self, superset: object, distance: numpy.ndarray, profile: AlphaProfile, amplitudes: numpy.ndarray, z_to_itr: object):
        """
        Args:
            profile (AlphaProfile): The profile to propagate.
            amplitude: The modes amplitudes.
        """
        self.superset = superset
        self.supermodes = superset.supermodes
        self.z_to_itr = z_to_itr
        self.distance = distance
        self.profile = profile
        self.amplitudes = amplitudes
        self.z_to_itr = self.profile.get_itr_vs_distance_interpolation()

    def plot(self, sub_sampling: int = 5, show_energy: bool = True, show_amplitudes: bool = True, **kwargs: dict) -> NoReturn:
        """
        Plots the propagation of amplitudes over a given profile, optionally showing energy and amplitude plots.

        Args:
            sub_sampling (int): The factor for sub-sampling data for plotting. Defaults to 5.
            show_energy (bool): Whether to plot the energy of the modes. Defaults to True.
            show_amplitudes (bool): Whether to plot the real part of the amplitudes. Defaults to True.
            **kwargs (dict): Additional keyword arguments for solver.

        """
        fig, ax = plt.subplots()
        ax.set(
            xlabel='Propagation distance (z)',
            ylabel='Inverse taper ratio (ITR)'
        )

        x_values = self.distance[::sub_sampling]

        for idx, mode in enumerate(self.supermodes):
            color = f"C{idx}"
            y_energy = numpy.abs(self.amplitudes[idx, ::sub_sampling]) ** 2
            y_amplitude = self.amplitudes[idx, ::sub_sampling].real

            if show_energy:
                ax.plot(x_values, y_energy, label=f'{mode.stylized_label} Energy', linewidth=2.0, linestyle='-', color=color)
            if show_amplitudes:
                ax.plot(x_values, y_amplitude, label=f'{mode.stylized_label} Amplitude', linewidth=2.0, linestyle='--', color=color)

        if show_energy:
            total_energy = numpy.sqrt(numpy.sum(numpy.abs(self.amplitudes) ** 2, axis=0))[::sub_sampling]
            ax.plot(x_values, total_energy, label='Total Energy', linewidth=3.0, linestyle='--', color='black')

        ax.legend()
        plt.show()

    def get_field_combination(self, amplitudes: numpy.ndarray, itr: int) -> numpy.ndarray:

        shape = self.superset.supermodes[0].field.data[0].shape
        total_field = numpy.zeros(shape).astype(complex)
        slice_number = 0

        for mode, amplitude in zip(self.superset.supermodes, amplitudes):
            total_field += amplitude * mode.field.data[slice_number]

        return total_field

    def generate_gif(
            self, *,
            sub_sampling: int = 5,
            mutliplicative_factor: float = -100,
            delta_azimuth: float = 0,
            save_directory: str = 'new_figure.gif',
            colormap: str = 'bwr',
            **kwargs) -> NoReturn:
        """
        Generates a gif video of the mode propagation.

        Args:
            sub_sampling (int): Propagation undersampling factor for the video production.
            mutliplicative_factor (float): Multiplicative factor for scaling.
            save_directory (str): Directory to save the generated GIF.
            delta_azimuth (float): Azimuthal change per frame.
            **kwargs (dict): Additional keyword arguments for solver.

        Returns:
            Tuple: Propagation distances, amplitudes, and ITR list.
        """
        sub_amplitudes = self.amplitudes[:, ::sub_sampling]
        initial_amplitudes = self.amplitudes[:, 0]

        sub_distance = self.distance[::sub_sampling]
        sub_itr = self.z_to_itr(sub_distance)

        total_field = self.get_field_combination(amplitudes=initial_amplitudes, itr=1.0)

        x, y = numpy.mgrid[0: total_field.shape[0], 0: total_field.shape[1]]
        grid = pyvista.StructuredGrid(x, y, total_field)

        plotter = pyvista.Plotter(notebook=False, off_screen=True)
        plotter.open_gif(save_directory, fps=20)
        plotter.view_isometric()
        # plotter.set_background('black', top='white')

        plotter.add_mesh(
            grid,
            scalars=total_field,
            style='surface',
            show_edges=True,
            edge_color='k',
            colormap=colormap,
            show_scalar_bar=False,
            clim=[-100, 100]
        )

        pts = grid.points.copy()
        azimuth = 0
        for z, amplitudes, itr in zip(sub_distance, sub_amplitudes.T, sub_itr):
            print(f'itr: {itr}')
            plotter.camera.elevation = -20
            plotter.camera.azimuth = azimuth
            azimuth += delta_azimuth

            structure = self.superset.get_slice_structure(itr=itr, add_symmetries=True)
            total_field = structure.get_field_combination(amplitudes, Linf_normalization=True) * mutliplicative_factor

            pts[:, -1] = total_field.T.ravel()
            plotter.update_coordinates(pts, render=True)
            plotter.update_scalars(total_field.T.ravel(), render=False)
            plotter.add_title(f'ITR: {itr: .3f}\t  z: {z: .3e}', font='courier', color='w', font_size=20)

            plotter.write_frame()

        plotter.close()

# -
