#!/usr/bin/env python
# -*- coding: utf-8 -*-

from typing import Tuple
import numpy
from scipy.interpolate import interp1d
import matplotlib.pyplot as plt

from SuPyMode import representation
from MPSPlots.render2D import SceneList, Axis
from matplotlib.animation import FuncAnimation, PillowWriter

from dataclasses import dataclass, field


@dataclass
class TaperSection():
    """
    A class to represent a taper section in optical fiber simulations.

    Attributes:
        z_array (np.ndarray): The array of longitudinal positions along the taper (z-coordinates).
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
    heating_length_initial: float = None
    heating_length_final: float = None

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


@dataclass
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
    initial_radius: float = 1
    n_point: int = 200
    symmetric: bool = False
    label: str = 'profile'
    add_end_of_taper_section: bool = True
    line_color: str = field(default='black', repr=False)
    line_style: str = field(default='--', repr=False)

    def __post_init__(self):
        self.section_list = []

    @property
    def first_section(self) -> TaperSection:
        """ Returns the first taper section added to the profile. """
        return self.section_list[0]

    @property
    def last_section(self) -> TaperSection:
        """ Returns the last taper section added to the profile. """
        return self.section_list[-1]

    def add_constant_segment(self, *, length: float, n_point: int = 100) -> None:
        """ Adds a constant section at the specified length with a given number of points. """
        section = self.get_constant_custom_section(
            length=length,
            rho=self.last_radius,
            start_z=self.last_z,
            n_point=n_point
        )

        self.section_list.append(section)

    def add_end_of_taper_segment(self, *, n_point: int = 100) -> None:
        """ Adds a constant segment at the end of the taper if the last section is not constant. """
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
        Add the constant section which length, radius and start position is to be provided

        :param      length:   Length of the constant section to be added
        :type       length:   float
        :param      radius:   Radius of the constant section to be added
        :type       radius:   float
        :param      start_z:  Initial z-position of the constant section to be added
        :type       start_z:  float
        :param      n_point:  The number of point where wo which evaluate that segment
        :type       n_point:  int

        :returns:   No returns
        :rtype:     None
        """
        z_array = numpy.linspace(start_z, length + start_z, n_point)

        radius_array = numpy.ones(n_point) * radius

        section = TaperSection(
            z_array=z_array,
            radius_array=radius_array,
        )

        return section

    def evaluate_adiabatic_factor(self, itr: numpy.ndarray) -> numpy.ndarray:
        interpolation = interp1d(
            x=self.itr_list,
            y=self.adiabatic,
            bounds_error=False,
            fill_value=numpy.nan
        )

        return interpolation(itr)

    def evaluate_distance_vs_itr(self, distance: numpy.ndarray) -> numpy.ndarray:
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
        Add a tapered section for a given alpha, initial_heating_length, initial_radius, stretching_length and starting z position

        :param      alpha:                  Alpha parameter which represent how the heating section changes in time
        :type       alpha:                  float
        :param      initial_heating_length: Initial length of the heating section
        :type       initial_heating_length: float
        :param      initial_radius:         The initial radius of the segment to be added
        :type       initial_radius:         float
        :param      stretching_length:      The total elongated lenght of the current segment to be added
        :type       stretching_length:      float
        :param      n_point:                The number of point where wo which evaluate that segment
        :type       n_point:                int

        :returns:   No returns
        :rtype:     None
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

    def compute_distance(self) -> numpy.ndarray:
        """
        Returns the distance array of the profile

        :returns:   The distance.
        :rtype:     numpy.ndarray
        """
        distance = numpy.linspace(0, self.last_z, self.n_point)

        self._distance = distance

    def compute_radius(self) -> numpy.ndarray:
        """
        Returns the array of radius [vs z-distance] for the taper structure

        :returns:   The ITR array
        :rtype:     numpy.ndarray
        """
        radius = self.compute_radius_from_segment_from_interpolation(self.distance)

        self._radius = radius

    def compute_itr_list(self) -> numpy.ndarray:
        """
        Returns the array of ITR value [vs z-distance] for the taper structure

        :returns:   The ITR array
        :rtype:     numpy.ndarray
        """
        itr_list = self.radius / self.initial_radius

        self._itr_list = itr_list

    def compute_adiabatic(self) -> numpy.ndarray:
        """
        Returns the array of adiabatc factor [vs ITR] for the taper structure

        .. math::
            f_c = \frac{1}{\rho} \frac{d \rho}{d z}

        :returns:   The adiabatic factor
        :rtype:     numpy.ndarray
        """
        dz = numpy.gradient(self.distance, axis=0, edge_order=2)

        ditr = numpy.gradient(numpy.log(self.radius), axis=0, edge_order=2)

        adiabatic = abs(ditr / dz)

        self._adiabatic = adiabatic

    def compute_taper_angle(self, symmetric: bool = None) -> numpy.ndarray:
        r"""
        Returns the array of taper angle for the taper structure
        From Tapered single-mode fibres and devices. Part 1: Adiabaticity criteria.
        Compute the adiabatic factor defined as:

        .. math::
            f_c = \frac{d \rho}{d z} = \Omega

        :returns:   The taper angle array
        :rtype:     numpy.ndarray
        """
        d_z = numpy.gradient(self.distance, axis=0, edge_order=2)

        d_rho = numpy.gradient(self.radius, axis=0, edge_order=2)

        taper_angle = abs(d_rho / d_z)

        self._taper_angle = taper_angle

    @property
    def smallest_itr(self) -> float:
        """
        Returns the smallest itr of the taper structure

        :returns:   Smallest itr value
        :rtype:     float
        """
        return self.itr_list.min()

    @property
    def last_z(self) -> float:
        """
        Returns the last, or equivalently the largest propagation distance computed

        :returns:   The z-distance
        :rtype:     float
        """
        if len(self.section_list) == 0:
            return 0
        else:
            return self.last_section.z_final

    @property
    def total_length(self) -> float:
        """
        Returns the total length of the component comprising the taper and constants sections.

        :returns:   { description_of_the_return_value }
        :rtype:     float
        """
        return self.last_section.z_final

    @property
    def last_radius(self) -> float:
        """
        Retunrs the radius value of the last z-position

        :returns:   The radius value
        :rtype:     float
        """
        if len(self.section_list) == 0:
            return self.initial_radius
        else:
            return self.last_section.radius_final

    def initialize(self) -> None:
        symmetric = self.symmetric

        self.symmetric = False

        if self.add_end_of_taper_section:
            self.add_end_of_taper_segment()

        self.compute_distance()
        self.compute_radius()
        self.compute_itr_list()
        self.compute_taper_angle()
        self.compute_adiabatic()

        self.symmetric = symmetric

    @property
    def distance(self):
        if self.symmetric:
            return numpy.linspace(self._distance[0], 2 * self._distance[-1], 2 * self._distance.size - 1)

        return self._distance

    @property
    def radius(self):
        if self.symmetric:
            return numpy.r_[self._radius, self._radius[-2::-1]]

        return self._radius

    @property
    def itr_list(self):
        if self.symmetric:
            return numpy.r_[self._itr_list, self._itr_list[-2::-1]]

        return self._itr_list

    @property
    def taper_angle(self):
        if self.symmetric:
            return numpy.r_[self._taper_angle, self._taper_angle[-2::-1]]

        return self._taper_angle

    @property
    def adiabatic(self):
        if self.symmetric:
            return numpy.r_[self._adiabatic, self._adiabatic[-2::-1]]

        return self._adiabatic

    def compute_radius_from_segment_from_interpolation(self, z: numpy.ndarray) -> numpy.ndarray:
        """
        Gets the radius of the component from all the interpolation segment.

        :param      z:    The distance array to which evaluate the radius of the component
        :type       z:    numpy.ndarray

        :returns:   The evaluated radius of the component vs the distance
        :rtype:     numpy.ndarray
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
        Add a tapered section following the previous one for a given alpha, initial_heating_length, stretching_length.

        :param      alpha:                  Alpha parameter which represent how the heating section changes in time
        :type       alpha:                  float
        :param      initial_heating_length: Initial length of the heating section
        :type       initial_heating_length: float
        :param      stretching_length:      The total elongated lenght of the current segment to be added
        :type       stretching_length:      float
        :param      n_point:                The number of point where wo which evaluate that segment
        :type       n_point:                int
        """
        return self.add_taper_custom_segment(
            alpha=alpha,
            initial_heating_length=initial_heating_length,
            initial_radius=self.last_radius,
            stretching_length=stretching_length,
            start_z=self.last_z,
            n_point=n_point
        )

    def get_itr_vs_distance_interpolation(self):
        return interp1d(
            x=self.distance,
            y=self.itr_list,
            bounds_error=False,
            fill_value=0
        )

    def single_plot(function):
        def wrapper(self, ax: Axis, line_color: str = None, line_style: str = None, **kwargs):
            line_style = self.line_style if line_style is None else line_style
            line_color = self.line_color if line_color is None else line_color

            x, y = function(
                self,
                ax=ax,
                line_color=line_color,
                line_style=line_style,
                **kwargs
            )

            ax.add_line(
                x=x,
                y=y,
                label=self.label,
                line_style=line_style,
                color=line_color
            )

        return wrapper

    @single_plot
    def render_itr_vs_z_on_ax(
            self,
            ax: Axis,
            line_style: str = None,
            line_color: str = None) -> None:
        """
        Add plot onto axis, the plots is ITR vs Z-distance

        :param      ax:   The axis on which to add the plot
        :type       ax:   Axis
        """
        ax.set_style(
            show_legend=False,
            y_limits=[0, None],
            x_label='Z-propagation [mm]',
            y_label='Inverse taper ratio [ITR]',
            x_scale_factor=1e3,
            y_scale="linear",
            line_width=2
        )

        return self.distance, self.itr_list

    @single_plot
    def render_taper_angle_vs_z_on_ax(
            self,
            ax: Axis,
            line_style: str = None,
            line_color: str = None) -> None:
        """
        Add plot onto axis, the plots is ITR vs Z-distance

        :param      ax:   The axis on which to add the plot
        :type       ax:   Axis
        """
        ax.set_style(
            show_legend=False,
            y_limits=[0, None],
            y_label='Taper angle [rad]',
            x_label='Z-propagation [mm]',
            x_scale_factor=1e3,
            y_scale="linear",
            line_width=2
        )

        return self.distance, self.taper_angle

    @single_plot
    def render_adiabatic_factor_vs_z_on_ax(
            self,
            ax: Axis,
            line_style: str = None,
            line_color: str = None) -> None:
        """
        Add plot onto axis, the plots is adiabatic criterion vs Z-distance

        :param      ax:   The axis on which to add the plot
        :type       ax:   Axis
        """
        ax.set_style(
            y_scale='log',
            y_label='Adiabatic criterion',
            x_label='z-distance'
        )

        return self.distance, self.adiabatic

    @single_plot
    def render_adiabatic_factor_vs_itr_on_ax(
            self,
            ax: Axis,
            line_style: str = None,
            line_color: str = None) -> None:
        """
        Add adiabatic criterion vs ITR plot to axis

        :param      ax:   The axis on which to add the plot
        :type       ax:   Axis
        """
        ax.set_style(**representation.adiabatic.Adiabatic.plot_style)

        return self.itr_list, self.adiabatic

    def plot(
            self,
            show_radius: bool = True,
            show_adiabatic: bool = True,
            show_taper_angle: bool = True) -> SceneList:
        """
        Generate two plots: ITR vs z distance and adiabatic criterion vs ITR

        :param      show_radius:       If True, plot shows radius as function of ITR
        :type       show_radius:       bool
        :param      show_adiabatic:    If True, plot shows adiabatic criterion as function of ITR
        :type       show_adiabatic:    bool
        :param      show_taper_angle:  If True, plot shows taper angle as function of ITR
        :type       show_taper_angle:  bool

        :returns:   The scene list.
        :rtype:     SceneList
        """
        self.initialize()

        figure = SceneList(
            title=f'Minimum ITR: {self.smallest_itr:.4f}',
            ax_orientation='vertical',
            unit_size=(8, 3)
        )

        if show_radius:
            ax = figure.append_ax()
            self.render_itr_vs_z_on_ax(ax=ax)

        if show_taper_angle:
            ax = figure.append_ax()
            self.render_taper_angle_vs_z_on_ax(ax=ax)

        if show_adiabatic:
            ax = figure.append_ax()
            self.render_adiabatic_factor_vs_itr_on_ax(ax=ax)

        figure.annotate_axis(position=(-0.15, 1.15))

        return figure

    def generate_propagation_gif(
            self,
            output_directory: str = './new_gif.gif',
            dpi: int = 100,
            fps: int = 20,
            number_of_frames: int = 200,
            dark_background: bool = True) -> None:
        """
        Genrates gif of the propagation of light into the taper structure

        :param      output_directory:  The output directory
        :type       output_directory:  str
        :param      dpi:               The dpi
        :type       dpi:               int
        :param      fps:               The fps [frame per seconde]
        :type       fps:               int
        :param      number_of_frames:  The number of frame
        :type       number_of_frames:  int
        :param      dark_background:   If True the background is black
        :type       dark_background:   bool

        :returns:   No returns
        :rtype:     None
        """
        self.initialize()

        figure, ax = plt.subplots(1, 1, figsize=(12, 6))

        ax.set_xlabel('Propagation axis [mm]', color='white')

        sub_sampling_factor = int(self.distance.size / number_of_frames)

        sub_distance = self.distance[::sub_sampling_factor] * 1e3
        sub_radius = self.radius[::sub_sampling_factor]
        sub_itr_list = self.itr_list[::sub_sampling_factor]

        if dark_background:
            style = plt.style.context("dark_background")
            ax.tick_params(colors='white', direction='out')
        else:
            style = plt.style.context('default')

        with style:
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

                line_1 = ax.axvline(
                    position,
                    linestyle='--',
                    linewidth=2,
                    color='red'
                )

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
