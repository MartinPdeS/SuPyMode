#!/usr/bin/env python
# -*- coding: utf-8 -*-

from typing import Tuple, Callable, Optional
from pydantic.dataclasses import dataclass
from dataclasses import field
from pydantic import ConfigDict
import numpy
from scipy.interpolate import interp1d
import matplotlib.pyplot as plt

import matplotlib.ticker as ticker
from matplotlib.animation import FuncAnimation, PillowWriter
from SuPyMode.helper import simple_plot_helper


config_dict = ConfigDict(
    extra='forbid',
    strict=True,
    arbitrary_types_allowed=True,
    kw_only=True,
)


@dataclass(config=config_dict)
class TaperSection():
    """
    A class to represent a taper section in optical fiber simulations.

    Attributes:
        z_array  (np.ndarray): The array of longitudinal positions along the taper (z-coordinates).
        radius_array (np.ndarray): The array of taper radii corresponding to the z positions.
        heating_length_initial (float, optional): The initial heating length of the taper section.
        heating_length_final (float, optional): The final heating length of the taper section.

    Properties:
        z_initial (float): Returns the initial z position of the taper.
        z_final (float): Returns the final z position of the taper.
        radius_initial (float): Returns the initial radius at the start of the taper.
        radius_final (float): Returns the radius at the end of the taper.
        is_constant (bool): Determines if the taper's radius is constant throughout.
        interpolation (callable): Provides an interpolation function for the radius over z.
    """

    z_array: numpy.ndarray
    radius_array: numpy.ndarray
    heating_length_initial: Optional[float] = None
    heating_length_final: Optional[float] = None

    @property
    def z_initial(self) -> float:
        """ Returns the initial z-coordinate of the taper section. """
        return self.z_array[0]

    @property
    def is_constant(self) -> float:
        """ Checks if the taper section's radius remains constant over its length. """
        return self.radius_array[0] == self.radius_array[-1]

    @property
    def z_final(self) -> float:
        """ Returns the final z-coordinate of the taper section. """
        return self.z_array[-1]

    @property
    def radius_initial(self) -> float:
        """ Returns the initial radius of the taper section. """
        return self.radius_array[0]

    @property
    def radius_final(self) -> float:
        """ Returns the final radius of the taper section. """
        return self.radius_array[-1]

    @property
    def interpolation(self):
        """
        Provides an interpolation function for radius as a function of z-coordinate.

        Returns:
            interp1d: An interpolator that estimates the radius at any z within the bounds
                      of z_array, with extrapolation set to zero outside the bounds.
        """
        return interp1d(
            x=self.z_array,
            y=self.radius_array,
            bounds_error=False,
            fill_value=0
        )


@dataclass(config=config_dict)
class AlphaProfile():
    r"""
    Represents a Gaussian profile for an optical fiber coupler.

    Translation table from article to class:
        - :math:`rho_w` = radius_segment
        - :math:`rho_0` = initial_radius
        - :math:`l_w` = heating_length_segment
        - :math:`x_0` = stretching_length

    Attributes:
        initial_radius (float): Initial radius of the taper structure, defaults to 1.
        n_point (int): Number of points for differential equation resolution, recommended 200+.
        symmetric (bool): If true, the taper structure is considered symmetric about z.
        label (str): Label for the profile, used in plotting.
        add_end_of_taper_section (bool): If true, adds a constant section at the end of the taper.
        line_color (str): Line color for plots, not part of the main data model.
        line_style (str): Line style for plots, not part of the main data model.
    """
    initial_radius: Optional[float] = 1
    n_point: Optional[int] = 200
    symmetric: Optional[bool] = False
    label: Optional[str] = 'profile'
    add_end_of_taper_section: Optional[bool] = True
    line_color: Optional[str] = field(default='black', repr=False)
    line_style: Optional[str] = field(default='--', repr=False)

    def __post_init__(self):
        """
        Initialize the section list after the dataclass fields are set. This method
        prepares the profile object for further operations such as adding sections
        or computing properties.
        """
        self.section_list = []

    @property
    def first_section(self) -> TaperSection:
        """
        Retrieves the first taper section added to the profile.

        Returns:
            TaperSection: The first taper section object in the section list.
        """
        return self.section_list[0]

    @property
    def last_section(self) -> TaperSection:
        """
        Retrieves the last taper section added to the profile.

        Returns:
            TaperSection: The last taper section object in the section list.
        """
        return self.section_list[-1]

    def add_constant_segment(self, *, length: float, n_point: int = 100) -> None:
        """
        Adds a constant section at the specified length with a given number of points to the profile.

        Parameters:
            length (float): Length of the constant section to be added.
            n_point (int): Number of points along the section for detailed resolution, defaults to 100.

        Returns:
            None
        """
        section = self.get_constant_custom_section(
            length=length,
            rho=self.last_radius,
            start_z=self.last_z,
            n_point=n_point
        )

        self.section_list.append(section)

    def add_end_of_taper_segment(self, *, n_point: int = 100) -> None:
        """
        Adds a constant segment at the end of the taper if the last section is not constant. This method
        ensures the taper ends smoothly or extends as needed.

        Parameters:
            n_point (int): Number of points along the section, defaults to 100.

        Returns:
            None
        """
        if self.last_section.is_constant:
            return

        length = self.last_section.heating_length_final / 2

        section = self.get_constant_custom_section(
            length=length,
            radius=self.last_radius,
            start_z=self.last_z,
            n_point=n_point
        )

        self.section_list.append(section)

    def get_constant_custom_section(
            self, *,
            length: float,
            radius: float,
            start_z: float = 0,
            n_point: int = 100) -> None:
        """
        Creates a constant section with specified parameters and adds it to the profile.

        Parameters:
            length (float): Length of the constant section to be added.
            radius (float): Radius of the section.
            start_z (float): Starting z-position for the section, defaults to 0.
            n_point (int): Number of points for resolution within the section, defaults to 100.

        Returns:
            TaperSection: The newly created taper section.
        """
        z_array = numpy.linspace(start_z, length + start_z, n_point)

        radius_array = numpy.ones(n_point) * radius

        section = TaperSection(
            z_array=z_array,
            radius_array=radius_array,
        )

        return section

    def evaluate_adiabatic_factor(self, itr: numpy.ndarray) -> numpy.ndarray:
        """
        Evaluates the adiabatic factor for given inverse taper ratios (ITR).

        Parameters:
            itr (numpy.ndarray): Array of inverse taper ratios.

        Returns:
            numpy.ndarray: Array of adiabatic factors corresponding to the provided ITRs.
        """
        interpolation = interp1d(
            x=self.itr_list,
            y=self.adiabatic,
            bounds_error=False,
            fill_value=numpy.nan
        )

        return interpolation(itr)

    def evaluate_distance_vs_itr(self, distance: numpy.ndarray) -> numpy.ndarray:
        """
        Evaluates the function of distance versus inverse taper ratio using interpolation.

        Parameters:
            distance (numpy.ndarray): Array of distances at which the ITR needs to be evaluated.

        Returns:
            numpy.ndarray: Array of ITRs at the specified distances.
        """
        interpolation = interp1d(
            x=self.itr_list,
            y=self.distance,
            bounds_error=True,
        )

        return interpolation(distance)

    def compute_radius_from_segment(
            self, *,
            alpha: float,
            initial_heating_length: float,
            stretching_length: float,
            initial_radius: float,
            distance: numpy.ndarray) -> Tuple[numpy.ndarray, float, float]:
        """
        Computes the radius as a function of the distance for a specific segment,
        applying a tapering formula based on the provided parameters.

        Args:
            alpha (float): Rate at which the heating section's influence changes over time.
            initial_heating_length (float): Initial length of the heating section.
            initial_radius (float): Radius at the start of the segment.
            stretching_length (float): Total length over which the segment is elongated.
            distance (numpy.ndarray): Array representing the z-distance.

        Returns:
            Tuple[numpy.ndarray, float, float]: A tuple containing:
                - radius (numpy.ndarray): Computed radius at each point in 'distance'.
                - final_radius (float): Radius at the end of the segment.
                - final_heating_length (float): Total length of the heating section after stretching.

        Raises:
            ValueError: If input conditions are not physically or mathematically valid.
        """
        self.assert_conditions(
            alpha=alpha,
            stretching_length=stretching_length,
            initial_heating_length=initial_heating_length
        )

        term0 = 2 * alpha * distance
        term2 = (1 - alpha) * initial_heating_length
        term3 = -1 / (2 * alpha)

        radius = initial_radius * (1 + term0 / term2)**term3
        final_radius = initial_radius * (1 + alpha * stretching_length / initial_heating_length)**(-1 / (2 * alpha))
        final_heating_length = initial_heating_length + alpha * stretching_length

        if numpy.any(radius <= 0):
            raise ValueError("Computed radius values contain non-physical negative or zero values.")

        return radius, final_radius, final_heating_length

    def assert_conditions(
            self, *,
            alpha: float,
            stretching_length: float,
            initial_heating_length: float) -> None:
        """
        Validates conditions for computing the taper segment.

        Args:
            alpha (float): Alpha parameter, non-zero to avoid division by zero.
            stretching_length (float): Length over which the segment is elongated.
            initial_heating_length (float): Initial length of the heating section.

        Raises:
            ValueError: If any condition that ensures a physically viable profile is violated.
        """
        if initial_heating_length <= 0:
            raise ValueError("Initial heating length must be positive.")

        if alpha == 0:
            raise ValueError("Alpha must not be zero to avoid division by zero in formula.")

        if alpha < 0 and stretching_length >= initial_heating_length / abs(alpha):
            raise ValueError("Stretching length for negative alpha exceeds the physically viable limit.")

    def add_taper_custom_segment(
            self, *,
            alpha: float,
            initial_heating_length: float,
            initial_radius: float,
            stretching_length: float,
            start_z: float = 0,
            n_point: int = 100) -> None:
        """
        Adds a custom tapered section to the profile based on specified parameters. This method is useful for creating
        detailed and specific taper geometries within the optical fiber.

        Parameters:
            alpha (float): Rate at which the heating section's influence changes over time.
            initial_heating_length (float): Initial length of the heating section before any stretching.
            initial_radius (float): Initial radius at the start of the taper segment.
            stretching_length (float): Length over which the taper is stretched.
            start_z (float): Starting z-coordinate for the taper segment, defaults to 0.
            n_point (int): Number of points along the taper for resolution, defaults to 100.

        Returns:
            None
        """
        alpha = 0.01 if alpha == 0 else alpha

        z_0 = (1 - alpha) * stretching_length / 2

        distance = numpy.linspace(0, z_0, n_point)

        assert distance[0] == 0, "Computation of taper section takes z as a reference and thus has to start with 0."

        radius, final_radius, final_heating_length = self.compute_radius_from_segment(
            alpha=alpha,
            initial_heating_length=initial_heating_length,
            stretching_length=stretching_length,
            initial_radius=initial_radius,
            distance=distance
        )

        section = TaperSection(
            z_array=distance + start_z,
            radius_array=radius,
            heating_length_initial=initial_heating_length,
            heating_length_final=final_heating_length
        )

        self.section_list.append(section)

    def compute_adiabatic(self, distance: numpy.ndarray, radius: numpy.ndarray) -> numpy.ndarray:
        """
        Computes the adiabatic factor, a measure of how gradually a taper changes, which is crucial for ensuring minimal
        light loss due to mode conversion.

        .. math::
            f_c = \frac{1}{\rho} \frac{d \rho}{d z}

        Parameters:
            distance (numpy.ndarray): Array of distances along the taper.
            radius (numpy.ndarray): Array of radii corresponding to each distance.

        Returns:
            numpy.ndarray: Array of adiabatic factors for the given distances and radii.
        """
        dz = numpy.gradient(distance, axis=0, edge_order=2)

        ditr = numpy.gradient(numpy.log(radius), axis=0, edge_order=2)

        return abs(ditr / dz)

    def compute_taper_angle(self, distance: numpy.ndarray, radius: numpy.ndarray) -> numpy.ndarray:
        r"""
        Computes the taper angle for a given set of distances and corresponding radii in the taper structure. This angle
        can provide insights into the performance and behavior of the taper under different conditions.
        From Tapered single-mode fibres and devices. Part 1: Adiabaticity criteria.

        .. math::
            f_c = \frac{d \rho}{d z} = \Omega


        Parameters:
            distance (numpy.ndarray): Array of distances along the taper.
            radius (numpy.ndarray): Array of radii corresponding to each distance.

        Returns:
            numpy.ndarray: Array of taper angles calculated from the rate of change of the radius with respect to the distance.
        """
        d_z = numpy.gradient(distance, axis=0, edge_order=2)

        d_rho = numpy.gradient(radius, axis=0, edge_order=2)

        return abs(d_rho / d_z)

    @property
    def smallest_itr(self) -> float:
        """
        Retrieves the smallest inverse taper ratio (ITR) of the taper structure. The smallest ITR can indicate the
        tightest part of the taper, which is critical for applications requiring precise control over light propagation.

        Returns:
            float: The smallest ITR value found in the taper profile.
        """
        return self.itr_list.min()

    @property
    def last_z(self) -> float:
        """
        Retrieves the last, or maximum, z-coordinate computed for the taper sections, which represents the end point
        of the taper structure.

        Returns:
            float: The last z-coordinate value in the taper structure.
        """
        if len(self.section_list) == 0:
            return 0
        else:
            return self.last_section.z_final

    @property
    def total_length(self) -> float:
        """
        Computes the total length of the taper structure, including both taper and constant sections, providing an overall
        size of the taper profile.

        Returns:
            float: Total length of the taper structure.
        """
        return self.last_section.z_final

    @property
    def last_radius(self) -> float:
        """
        Retrieves the radius at the last computed z-position, which represents the end radius of the taper structure.

        Returns:
            float: Radius at the last z-coordinate in the taper profile.
        """
        if len(self.section_list) == 0:
            return self.initial_radius
        else:
            return self.last_section.radius_final

    def initialize(self) -> None:
        """
        Initializes or re-initializes the profile, typically called after making modifications to the taper sections. This method
        recalculates and updates internal parameters to reflect the current state of the taper structure.

        Returns:
            None
        """
        if self.add_end_of_taper_section:
            self.add_end_of_taper_segment()

        distance = numpy.linspace(0, self.last_z, self.n_point)
        radius = self.compute_radius_from_segment_from_interpolation(distance)
        itr_list = radius / self.initial_radius
        adiabatic = self.compute_adiabatic(distance=distance, radius=radius)
        taper_angle = self.compute_taper_angle(distance=distance, radius=radius)

        if self.symmetric:
            self._distance = numpy.linspace(distance[0], 2 * distance[-1], 2 * distance.size - 1)
            self._itr_list = numpy.r_[itr_list, itr_list[-2::-1]]
            self._radius = numpy.r_[radius, radius[-2::-1]]
            self._adiabatic = numpy.r_[adiabatic, adiabatic[-2::-1]]
            self._taper_angle = numpy.r_[taper_angle, taper_angle[-2::-1]]
        else:
            self._distance = distance
            self._radius = radius
            self._itr_list = itr_list
            self._adiabatic = adiabatic
            self._taper_angle = taper_angle

    @property
    def distance(self) -> numpy.ndarray:
        try:
            return self._distance
        except AttributeError:
            raise AttributeError('Profile has not been initialized yet. The user need to run the initialize() method first. ')

    @property
    def radius(self) -> numpy.ndarray:
        try:
            return self._radius
        except AttributeError:
            raise AttributeError('Profile has not been initialized yet. The user need to run the initialize() method first. ')

    @property
    def itr_list(self) -> numpy.ndarray:
        try:
            return self._itr_list
        except AttributeError:
            raise AttributeError('Profile has not been initialized yet. The user need to run the initialize() method first. ')

    @property
    def adiabatic(self) -> numpy.ndarray:
        try:
            return self._adiabatic
        except AttributeError:
            raise AttributeError('Profile has not been initialized yet. The user need to run the initialize() method first. ')

    @property
    def taper_angle(self) -> numpy.ndarray:
        try:
            return self._taper_angle
        except AttributeError:
            raise AttributeError('Profile has not been initialized yet. The user need to run the initialize() method first. ')

    def compute_radius_from_segment_from_interpolation(self, z: numpy.ndarray) -> numpy.ndarray:
        """
        Computes the radius at specified z-distances based on interpolation from existing taper sections, providing a continuous
        profile of the radius along the taper.

        Parameters:
            z (numpy.ndarray): Array of z-distances at which to evaluate the radius.

        Returns:
            numpy.ndarray: Array of radius values interpolated along the given z-distances.
        """
        radius = numpy.zeros(z.size)

        for section in self.section_list:
            evaluation = section.interpolation(z)
            idx_non_null = evaluation != 0
            radius[idx_non_null] = evaluation[idx_non_null]

        return radius

    def add_taper_segment(
            self, *,
            alpha: float,
            initial_heating_length: float,
            stretching_length: float,
            initial_radius: float = None,
            n_point: int = 100) -> None:
        """
        Adds a tapered section following the previously defined section, using provided taper parameters to define the new section's
        geometry and behavior.

        Parameters:
            alpha (float): Rate at which the heating section's influence changes over time.
            initial_heating_length (float): Initial length of the heating section.
            stretching_length (float): Total length over which the segment is stretched.
            initial_radius (float): Radius at the start of the segment, defaults to the last radius if not provided.
            n_point (int): Number of points along the section for resolution, defaults to 100.

        Returns:
            None
        """
        return self.add_taper_custom_segment(
            alpha=alpha,
            initial_heating_length=initial_heating_length,
            initial_radius=self.last_radius,
            stretching_length=stretching_length,
            start_z=self.last_z,
            n_point=n_point
        )

    def get_itr_vs_distance_interpolation(self) -> Callable:
        """
        Generates an interpolation function for inverse taper ratio (ITR) as a function of distance.
        This allows for quick lookups of ITR at arbitrary distances along the taper.

        Returns:
            Callable: A function that interpolates ITR based on given distances.
        """
        return interp1d(
            x=self.distance,
            y=self.itr_list,
            bounds_error=False,
            fill_value=0
        )

    def single_plot(function) -> Callable:
        """
        Decorator to apply a standard plotting style to any plotting function within this class.
        It automatically sets line styles, colors, and labels from the class attributes.

        Parameters:
            function (Callable): The plotting function to decorate.

        Returns:
            Callable: A wrapped plotting function that integrates additional styling and annotations.
        """

        def wrapper(self, ax: plt.Axes, line_color: str = None, line_style: str = None, **kwargs):
            line_style = self.line_style if line_style is None else line_style
            line_color = self.line_color if line_color is None else line_color

            x, y = function(self, ax=ax, line_color=line_color, line_style=line_style, **kwargs)

            ax.plot(x, y, label=self.label, linestyle=line_style, color=line_color, linewidth=2)

        return wrapper

    @single_plot
    def render_itr_vs_z_on_ax(self, ax: plt.Axes, **kwargs) -> tuple[numpy.ndarray, numpy.ndarray]:
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
            xlabel='Z-propagation [mm]',
            ylabel='Inverse taper ratio [ITR]',
            yscale="linear",
        )

        scale_x = 1e-3
        ticks_x = ticker.FuncFormatter(lambda x, pos: '{0:g}'.format(x / scale_x))
        ax.xaxis.set_major_formatter(ticks_x)

        ax.legend()

        return self.distance, self.itr_list

    @single_plot
    def render_taper_angle_vs_z_on_ax(self, ax: plt.Axes, **kwargs) -> tuple[numpy.ndarray, numpy.ndarray]:
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
            xlabel='Z-propagation [mm]',
            ylabel='Taper angle [rad]',
            yscale="linear",
        )

        scale_x = 1e-3
        ticks_x = ticker.FuncFormatter(lambda x, pos: '{0:g}'.format(x / scale_x))
        ax.xaxis.set_major_formatter(ticks_x)

        ax.legend()

        return self.distance, self.taper_angle

    @single_plot
    def render_adiabatic_factor_vs_z_on_ax(self, ax: plt.Axes, **kwargs) -> tuple[numpy.ndarray, numpy.ndarray]:
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
            xlabel='Z-propagation [mm]',
            ylabel='Adiabatic criterion',
            yscale="linear",
        )

        scale_x = 1e-3
        ticks_x = ticker.FuncFormatter(lambda x, pos: '{0:g}'.format(x / scale_x))
        ax.xaxis.set_major_formatter(ticks_x)

        ax.legend()

        return self.distance, self.adiabatic

    @single_plot
    def render_adiabatic_factor_vs_itr_on_ax(self, ax: plt.Axes, **kwargs) -> tuple[numpy.ndarray, numpy.ndarray]:
        """
        Add adiabatic criterion vs ITR plot to axis

        :param      ax:   The axis on which to add the plot
        :type       ax:   plt.Axes
        """
        return self.itr_list, self.adiabatic

    @simple_plot_helper
    def plot(
            self,
            ax: plt.Axes,
            show_radius: bool = True,
            show_adiabatic: bool = True,
            show_taper_angle: bool = True) -> None:
        """
        Generates plots based on the current state of the taper profile. This can include plots of radius vs. z-distance,
        adiabatic factor vs. ITR, and taper angle vs. z-distance, based on the specified flags.

        Parameters:
            show_radius (bool): If True, includes a plot of radius vs. z-distance.
            show_adiabatic (bool): If True, includes a plot of adiabatic factor vs. ITR.
            show_taper_angle (bool): If True, includes a plot of taper angle vs. z-distance.
        """
        ax.set_title(f'Minimum ITR: {self.smallest_itr:.4f}')

        if show_radius:
            self.render_itr_vs_z_on_ax(ax=ax)

        if show_taper_angle:
            self.render_taper_angle_vs_z_on_ax(ax=ax)

        if show_adiabatic:
            self.render_adiabatic_factor_vs_itr_on_ax(ax=ax)

    def generate_propagation_gif(
            self,
            output_directory: str = './new_gif.gif',
            dpi: int = 100,
            fps: int = 20,
            number_of_frames: int = 200,
            dark_background: bool = True) -> None:
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
        ax.set_xlabel('Propagation axis [mm]', color='white' if dark_background else 'black')
        style_context = "dark_background" if dark_background else "default"

        sub_sampling_factor = int(self.distance.size / number_of_frames)

        sub_distance = self.distance[::sub_sampling_factor] * 1e3
        sub_radius = self.radius[::sub_sampling_factor]
        sub_itr_list = self.itr_list[::sub_sampling_factor]

        with plt.style.context(style_context):
            def init_func() -> tuple:
                line_0 = ax.plot(sub_distance, sub_radius, color='black')
                line_1 = ax.plot(sub_distance, -sub_radius, color='black')

                line_2 = ax.fill_between(
                    sub_distance,
                    +sub_radius,
                    -sub_radius,
                    color='lightblue',
                    alpha=0.8
                )

                return [*line_0, *line_1, line_2]

            def animate(slice_number: int) -> tuple:
                position = sub_distance[slice_number]
                itr = sub_itr_list[slice_number]
                title = f'[slice: {slice_number} - ITR: {itr:.3f}]'

                if slice_number > 0:
                    ax.lines[-1].remove()

                line_0 = ax.set_title(title, color='white')

                line_1 = ax.axvline(position, linestyle='--', linewidth=2, color='red')

                return line_0, line_1

        animation = FuncAnimation(
            fig=figure,
            func=animate,
            init_func=init_func,
            blit=True,
            repeat=True,
            frames=sub_itr_list.size
        )

        animation.save(
            output_directory,
            dpi=dpi,
            writer=PillowWriter(fps=fps)
        )

# -
