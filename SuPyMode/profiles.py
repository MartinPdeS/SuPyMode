#!/usr/bin/env python
# -*- coding: utf-8 -*-

from typing import Callable

import numpy
import matplotlib.pyplot as plt

import matplotlib.ticker as ticker
from matplotlib.animation import FuncAnimation, PillowWriter
from MPSPlots.styles import mps


from SuPyMode.binary.interface_taper import AlphaProfile


class AlphaProfile(AlphaProfile):
    line_style = "--"
    line_color = "black"

    def single_plot(function) -> Callable:
        """
        Decorator to apply a standard plotting style to any plotting function within this class.
        It automatically sets line styles, colors, and labels from the class attributes.

        Parameters:
            function (Callable): The plotting function to decorate.

        Returns:
            Callable: A wrapped plotting function that integrates additional styling and annotations.
        """

        def wrapper(
            self, ax: plt.Axes, line_color: str = None, line_style: str = None, **kwargs
        ):
            line_style = self.line_style if line_style is None else line_style
            line_color = self.line_color if line_color is None else line_color

            x, y = function(
                self, ax=ax, line_color=line_color, line_style=line_style, **kwargs
            )

            ax.plot(
                x,
                y,
                label=self.label,
                linestyle=line_style,
                color=line_color,
                linewidth=2,
            )

        return wrapper

    @single_plot
    def render_itr_vs_z_on_ax(
        self, ax: plt.Axes, **kwargs
    ) -> tuple[numpy.ndarray, numpy.ndarray]:
        """
        Renders a plot of inverse taper ratio (ITR) versus z-distance onto a given axis. This method is typically
        used for visualizing how the ITR changes along the length of the taper.

        Parameters:
            ax (plt.Axes): The matplotlib axis on which to plot.
            line_style (str, optional): Line style for the plot, defaults to class attribute.
            line_color (str, optional): Line color for the plot, defaults to class attribute.

        Returns:
            tuple[numpy.ndarray, numpy.ndarray]
        """
        ax.set(
            ylim=[0, None],
            xlabel="Z-propagation [mm]",
            ylabel="Taper radius [mm]",
            yscale="linear",
        )

        scale_x = 1e-3
        ticks_x = ticker.FuncFormatter(lambda x, pos: "{0:g}".format(x / scale_x))
        ax.xaxis.set_major_formatter(ticks_x)

        return self.distance, self.itr_list

    @single_plot
    def render_taper_angle_vs_z_on_ax(
        self, ax: plt.Axes, **kwargs
    ) -> tuple[numpy.ndarray, numpy.ndarray]:
        """
        Plots the taper angle as a function of z-distance on a provided axis. Useful for understanding the geometric
        changes in the taper profile over its length.

        Parameters:
            ax (plt.Axes): The matplotlib axis on which to plot.
            line_style (str, optional): Specifies the style of the plot line, if different from the class default.
            line_color (str, optional): Specifies the color of the plot line, if different from the class default.

        Returns:
            tuple[numpy.ndarray, numpy.ndarray]
        """
        ax.set(
            xlabel="Z-propagation [mm]",
            ylabel="Taper angle [rad]",
            yscale="linear",
        )

        scale_x = 1e-3
        ticks_x = ticker.FuncFormatter(lambda x, pos: "{0:g}".format(x / scale_x))
        ax.xaxis.set_major_formatter(ticks_x)

        return self.distance, self.taper_angle

    @single_plot
    def render_adiabatic_factor_vs_z_on_ax(
        self, ax: plt.Axes, **kwargs
    ) -> tuple[numpy.ndarray, numpy.ndarray]:
        """
        Plots the adiabatic criterion versus z-distance on the specified axis. This plot helps assess the
        efficiency and effectiveness of the taper in maintaining adiabatic conditions throughout its course.

        Parameters:
            ax (plt.Axes): The matplotlib axis on which to plot.
            line_style (str, optional): Line style for the plot, can override default.
            line_color (str, optional): Line color for the plot, can override default.

        Returns:
            tuple[numpy.ndarray, numpy.ndarray]
        """

        ax.set(
            xlabel="Z-propagation [mm]",
            ylabel="Adiabatic criterion",
            yscale="linear",
        )

        scale_x = 1e-3
        ticks_x = ticker.FuncFormatter(lambda x, pos: "{0:g}".format(x / scale_x))
        ax.xaxis.set_major_formatter(ticks_x)

        return self.distance, self.adiabatic

    @single_plot
    def render_adiabatic_factor_vs_itr_on_ax(
        self, ax: plt.Axes, **kwargs
    ) -> tuple[numpy.ndarray, numpy.ndarray]:
        """
        Add adiabatic criterion vs ITR plot to axis

        :param      ax:   The axis on which to add the plot
        :type       ax:   plt.Axes
        """
        ax.set(
            xlabel="ITR",
            ylabel="Adiabatic criterion",
            yscale="linear",
        )

        return self.itr_list, self.adiabatic

    def plot(
        self,
        show_radius: bool = True,
        show_adiabatic: bool = True,
        show_taper_angle: bool = True,
    ) -> None:
        """
        Generates plots based on the current state of the taper profile. This can include plots of radius vs. z-distance,
        adiabatic factor vs. ITR, and taper angle vs. z-distance, based on the specified flags.

        Parameters:
            show_radius (bool): If True, includes a plot of radius vs. z-distance.
            show_adiabatic (bool): If True, includes a plot of adiabatic factor vs. ITR.
            show_taper_angle (bool): If True, includes a plot of taper angle vs. z-distance.
        """
        shows = [show_radius, show_adiabatic, show_taper_angle]
        number_of_plots = sum(shows)

        methods = [
            self.render_itr_vs_z_on_ax,
            self.render_taper_angle_vs_z_on_ax,
            self.render_adiabatic_factor_vs_itr_on_ax,
        ]

        with plt.style.context(mps):
            figure, axes = plt.subplots(
                number_of_plots, 1, figsize=(10, 3 * number_of_plots), squeeze=False
            )

            axes = axes.flatten()

            for ax, show, method in zip(axes, shows, methods):
                if show:
                    ax.set_title(f"Minimum ITR: {self.smallest_itr:.4f}")
                    method(ax=ax)

        plt.show()

    def generate_propagation_gif(
        self,
        output_directory: str = "./new_gif.gif",
        dpi: int = 100,
        fps: int = 20,
        number_of_frames: int = 200,
        dark_background: bool = True,
    ) -> None:
        """
        Generates an animated GIF of light propagation in a taper structure.

        Parameters:
            output_directory (str): Path where the GIF will be saved.
            dpi (int): Dots per inch for the output GIF.
            fps (int): Frames per second for the animation.
            number_of_frames (int): Total number of frames in the animation.
            dark_background (bool): If True, use a dark background for the GIF.

        Returns:
            None
        """
        figure, ax = plt.subplots(1, 1, figsize=(12, 6))
        ax.set_xlabel(
            "Propagation axis [mm]", color="white" if dark_background else "black"
        )
        style_context = "dark_background" if dark_background else "default"

        sub_sampling_factor = int(self.distance.size / number_of_frames)

        sub_distance = self.distance[::sub_sampling_factor] * 1e3
        sub_radius = self.radius[::sub_sampling_factor]
        sub_itr_list = self.itr_list[::sub_sampling_factor]

        with plt.style.context(style_context):

            def init_func() -> tuple:
                line_0 = ax.plot(sub_distance, sub_radius, color="black")
                line_1 = ax.plot(sub_distance, -sub_radius, color="black")

                line_2 = ax.fill_between(
                    sub_distance, +sub_radius, -sub_radius, color="lightblue", alpha=0.8
                )

                return [*line_0, *line_1, line_2]

            def animate(slice_number: int) -> tuple:
                position = sub_distance[slice_number]
                itr = sub_itr_list[slice_number]
                title = f"[slice: {slice_number} - ITR: {itr:.3f}]"

                if slice_number > 0:
                    ax.lines[-1].remove()

                line_0 = ax.set_title(title, color="white")

                line_1 = ax.axvline(position, linestyle="--", linewidth=2, color="red")

                return line_0, line_1

        animation = FuncAnimation(
            fig=figure,
            func=animate,
            init_func=init_func,
            blit=True,
            repeat=True,
            frames=sub_itr_list.size,
        )

        animation.save(output_directory, dpi=dpi, writer=PillowWriter(fps=fps))


# -
