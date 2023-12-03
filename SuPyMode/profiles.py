#!/usr/bin/env python
# -*- coding: utf-8 -*-

import numpy
from scipy.interpolate import interp1d
import matplotlib.pyplot as plt

from SuPyMode import representation
from MPSPlots.render2D import SceneList, Axis
from matplotlib.animation import FuncAnimation, PillowWriter

from dataclasses import dataclass


@dataclass
class AlphaProfile():
    initial_radius: float = 1
    """ Initial radius of the taper structure """
    n_point: int = 200
    """ Number of point for solving the differential equation for the ITR vs distance. Keep it high [200+]"""
    symmetric: bool = False
    """ Bolean to defined if the taper structure is z-symmetric """
    label: str = 'profile'
    """ Label of the profile, shown as label of plots """
    line_color: str = 'black'
    """ Color of the lines for the plots """
    line_style: str = '--'
    """ Style of the lines for the plots """

    """
    Class represent the fiber structure coupler z-profile.
    This particular class is set to a Gaussian profile.
    Translation table:
        - rho_w = radius_segment
        - rho_0 = initial_radius
        - l_w = heating_length_segment
        - x_0 = stretching_length
    """

    def __post_init__(self):
        self.segment_interpolation = []
        self.z_segment = []
        self.radius_segments = []
        self.heating_length_segment = []

    def symmetrize_array(self, array: numpy.ndarray) -> numpy.ndarray:
        """
        Returns a symmetric array. The symmetry take place about the last point of the array.

        :param      array:  The input array
        :type       array:  numpy.ndarray

        :returns:   The symmetrized array
        :rtype:     numpy.ndarray
        """
        symmetric_array = numpy.r_[array, array[-2::-1]]

        return symmetric_array

    def get_distance(self, symmetric: bool = None) -> numpy.ndarray:
        if symmetric is None:
            symmetric = self.symmetric

        if symmetric:
            return numpy.linspace(0, 2 * self.last_z, 2 * self.n_point - 1)
        else:
            return numpy.linspace(0, self.last_z, self.n_point)

    def add_constant_segment(self, *, length: float, n_point: int = 100) -> None:
        """
        Add the constant section following the last section which length is to be evaluated.

        :param      length:   Length of the constant section to be added
        :type       length:   float
        :param      n_point:  The number of point where wo which evaluate that segment
        :type       n_point:  int
        """
        return self.add_constant_custom_section(
            length=length,
            rho=self.last_radius,
            start_z=self.last_z,
            n_point=n_point
        )

    def add_end_of_taper_segment(
            self, *,
            length: float = None,
            n_point: int = 100) -> None:
        """
        Add the constant section which length equal the final length of the
        heating section.

        :param      n_point:  The number of point where wo which evaluate that segment
        :type       n_point:  int
        """
        if length is None:
            length = self.last_heating_length / 2

        return self.add_constant_custom_section(
            length=length,
            radius=self.last_radius,
            start_z=self.last_z,
            n_point=n_point
        )

    def add_constant_custom_section(
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
        """
        z = numpy.linspace(0, length, n_point)

        radius = numpy.ones(n_point) * radius

        interpolation = interp1d(
            x=z + start_z,
            y=radius,
            bounds_error=False,
            fill_value=0
        )

        self.segment_interpolation.append(interpolation)
        self.z_segment.append(length + start_z)

    def evaluate_adiabatic_factor(self, itr: numpy.ndarray) -> numpy.ndarray:
        x_axis = self.get_itr_list()
        y_axis = self.get_adiabatic()

        interpolation = interp1d(
            x=x_axis,
            y=y_axis,
            bounds_error=False,
            fill_value=numpy.nan
        )

        return interpolation(itr)

    def evaluate_distance_vs_itr(self, distance: numpy.ndarray) -> numpy.ndarray:
        x_axis = self.get_itr_list()
        y_axis = self.get_distance()

        interpolation = interp1d(
            x=x_axis,
            y=y_axis,
            bounds_error=True,
        )

        return interpolation(distance)

    def get_radius_from_segment(
            self, *,
            alpha: float,
            initial_heating_length: float,
            stretching_length: float,
            initial_radius: float,
            distance: numpy.ndarray) -> tuple:
        """
        Gets the radius as a fonction of the distance for a specific segment.

        :param      alpha:                  Alpha parameter which represent how the heating section changes in time
        :type       alpha:                  float
        :param      initial_heating_length: Initial length of the heating section
        :type       initial_heating_length: float
        :param      initial_radius:         The initial radius of the segment to be added
        :type       initial_radius:         float
        :param      stretching_length:      The total elongated lenght of the current segment to be added
        :type       stretching_length:      float
        :param      distance:               Array representing the z-distance.
        :type       distance:               numpy.ndarray

        :returns:   The radius, final radius and final heating length
        :rtype:     tuple
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

        assert not numpy.any(radius < 0), "Negative radius value are not physical"

        return radius, final_radius, final_heating_length

    def assert_conditions(
            self, *,
            alpha: float,
            stretching_length: float,
            initial_heating_length: float) -> None:
        """
        Assert a few condition of viability of the proposed recipe.

        :param      alpha:                   The alpha
        :type       alpha:                   float
        :param      stretching_length:       The stretching length
        :type       stretching_length:       float
        :param      initial_heating_length:  The initial heating length
        :type       initial_heating_length:  float

        :returns:   No returns
        :rtype:     None
        """
        assert initial_heating_length > 0, "The initial heat lenght initial_heating_length cannot be negative!"

        if alpha < 0:
            assert stretching_length < initial_heating_length / abs(alpha), "Condition: x0 < initial_heating_length / |alpha| is not respected! see Birks article in the references!"

    def add_custom_segment(self, distance: numpy.ndarray, radius: numpy.ndarray) -> None:
        end_of_segment = distance[-1]
        final_radius_of_segment = radius[-1]

        interpolation = interp1d(
            x=distance,
            y=radius,
            bounds_error=False,
            fill_value=0
        )

        self.segment_interpolation.append(interpolation)
        self.z_segment.append(end_of_segment)
        self.radius_segments.append(final_radius_of_segment)

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
        """
        alpha = 0.01 if alpha == 0 else alpha

        z_0 = (1 - alpha) * stretching_length / 2
        end_of_segment = z_0 + start_z
        distance = numpy.linspace(0, z_0, n_point)

        assert distance[0] == 0, "Computation of taper section takes z as a reference and thus has to start with 0."

        radius, final_radius, final_heating_length = self.get_radius_from_segment(
            alpha=alpha,
            initial_heating_length=initial_heating_length,
            stretching_length=stretching_length,
            initial_radius=initial_radius,
            distance=distance
        )

        interpolation = interp1d(
            x=distance + start_z,
            y=radius,
            bounds_error=False,
            fill_value=0
        )

        self.segment_interpolation.append(interpolation)
        self.z_segment.append(end_of_segment)
        self.radius_segments.append(final_radius)
        self.heating_length_segment.append(final_heating_length)

    def interpret_symmetric(function):
        def wrapper(self, symmetric: bool = None) -> numpy.ndarray:
            if symmetric is None:
                symmetric = self.symmetric

            array: numpy.ndarray = function(self, symmetric=symmetric)

            if symmetric:
                return self.symmetrize_array(array)
            else:
                return array

            return array

        return wrapper

    @interpret_symmetric
    def get_radius(self, symmetric: bool = None) -> numpy.ndarray:
        """
        Returns the array of radius [vs z-distance] for the taper structure

        :returns:   The ITR array
        :rtype:     numpy.ndarray
        """
        distance = self.get_distance(symmetric=False)

        radius = self.get_radius_from_segment_from_interpolation(distance)

        return radius

    @interpret_symmetric
    def get_itr_list(self, symmetric: bool = None) -> numpy.ndarray:
        """
        Returns the array of ITR value [vs z-distance] for the taper structure

        :returns:   The ITR array
        :rtype:     numpy.ndarray
        """
        radius = self.get_radius(symmetric=False)

        itr_list = radius / self.initial_radius

        return itr_list

    @interpret_symmetric
    def get_adiabatic(self, symmetric: bool = None) -> numpy.ndarray:
        """
        Returns the array of adiabatc factor [vs ITR] for the taper structure

        .. math::
            f_c = \frac{1}{\rho} \frac{d \rho}{d z}

        :returns:   The adiabatic factor
        :rtype:     numpy.ndarray
        """
        distance = self.get_distance(symmetric=False)
        radius = self.get_radius(symmetric=False)

        dz = numpy.gradient(distance, axis=0, edge_order=2)

        ditr = numpy.gradient(numpy.log(radius), axis=0, edge_order=2)

        adiabatic = abs(ditr / dz)

        return adiabatic

    @interpret_symmetric
    def get_taper_angle(self, symmetric: bool = None) -> numpy.ndarray:
        r"""
        Returns the array of taper angle for the taper structure
        From Tapered single-mode fibres and devices. Part 1: Adiabaticity criteria.
        Compute the adiabatic factor defined as:

        .. math::
            f_c = \frac{d \rho}{d z} = \Omega

        :returns:   The taper angle array
        :rtype:     numpy.ndarray
        """
        distance = self.get_distance(symmetric=False)
        radius = self.get_radius(symmetric=False)

        d_z = numpy.gradient(distance, axis=0, edge_order=2)

        d_rho = numpy.gradient(radius, axis=0, edge_order=2)

        taper_angle = abs(d_rho / d_z)

        return taper_angle

    @property
    def smallest_itr(self) -> float:
        """
        Returns the smallest itr of the taper structure

        :returns:   Smallest itr value
        :rtype:     float
        """
        itr_list = self.get_itr_list()
        return itr_list.min()

    @property
    def last_z(self) -> float:
        """
        Returns the last, or equivalently the largest propagation distance computed

        :returns:   The z-distance
        :rtype:     float
        """
        if len(self.z_segment) == 0:
            return 0
        else:
            return self.z_segment[-1]

    @property
    def total_length(self) -> float:
        """
        Returns the total length of the component comprising the taper and constants sections.

        :returns:   { description_of_the_return_value }
        :rtype:     float
        """
        return self.get_distance()[-1]

    @property
    def last_radius(self) -> float:
        """
        Retunrs the radius value of the last z-position

        :returns:   The radius value
        :rtype:     float
        """
        if len(self.radius_segments) == 0:
            return self.initial_radius
        else:
            return self.radius_segments[-1]

    @property
    def last_heating_length(self) -> float:
        """
        Retunrs the heating lenght of the last z-position.

        :returns:   THe heating length
        :rtype:     float
        """
        return self.heating_length_segment[-1]

    def get_radius_from_segment_from_interpolation(self, z: numpy.ndarray) -> numpy.ndarray:
        """
        Gets the radius of the component from all the interpolation segment.

        :param      z:    The distance array to which evaluate the radius of the component
        :type       z:    numpy.ndarray

        :returns:   The evaluated radius of the component vs the distance
        :rtype:     numpy.ndarray
        """
        self.add_end_of_taper_segment()

        radius = numpy.zeros(z.size)

        for interpolation in self.segment_interpolation:
            evaluation = interpolation(z)
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
            x=self.get_distance(),
            y=self.get_itr_list(),
            bounds_error=False,
            fill_value=0
        )

    def get_length_scale(self, core_radius: float) -> numpy.ndarray:
        r"""
        Compute the adiabatic factor defined as:

        .. math::
            f_c = \frac{\rho}{\Omega}

        :returns:   The amplitudes as a function of the distance in the coupler
        :rtype:     numpy.ndarray
        """
        return core_radius / self.get_taper_angle()

    def single_plot(function):
        def wrapper(self, *args, ax: Axis, line_color: str = None, line_style: str = None, **kwargs):
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
        distance = self.get_distance()
        itr_list = self.get_itr_list()

        ax.set_style(
            show_legend=False,
            y_limits=[0, None],
            x_label='Z-propagation [mm]',
            y_label='Inverse taper ratio [ITR]',
            x_scale_factor=1e3,
            y_scale="linear",
            line_width=2
        )

        return distance, itr_list

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
        distance = self.get_distance()
        taper_angle = self.get_taper_angle()

        ax.set_style(
            show_legend=False,
            y_limits=[0, None],
            y_label='Taper angle [rad]',
            x_label='Z-propagation [mm]',
            x_scale_factor=1e3,
            y_scale="linear",
            line_width=2
        )

        return distance, taper_angle

    @single_plot
    def render_taper_length_scale_vs_itr_on_ax(
            self,
            ax: Axis,
            line_style: str = None,
            line_color: str = None) -> None:
        """
        Add ITR vs Z-distance plot to axis

        :param      ax:   The axis on which to add the plot
        :type       ax:   Axis
        """
        itr_list = self.get_itr_list()
        length_scale = self.get_length_scale(core_radius=self.initial_radius)

        ax.set_style(
            show_legend=False,
            y_label='Coupler length-scale',
            x_label='Inverse taper ratio [ITR]',
            x_scale_factor=1e3,
            y_scale="linear",
            line_width=2
        )

        return itr_list, length_scale

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
        distance = self.get_distance()
        adiabatic = self.get_adiabatic()

        ax.set_style(
            y_scale='log',
            y_label='Adiabatic criterion',
            x_label='z-distance'
        )

        return distance, adiabatic

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
        ax.set_style(**representation.adiabatic.ax_style)

        itr_list = self.get_itr_list()
        adiabatic = self.get_adiabatic()

        return itr_list, adiabatic

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

        figure, ax = plt.subplots(1, 1, figsize=(12, 6))
        ax.set_xlabel('z-distance', color='white')

        sub_sampling_factor = numpy.ceil(self.distance.size / number_of_frames).astype(int)

        sub_distance = self.distance[::sub_sampling_factor]
        sub_radius = self.radius[::sub_sampling_factor]
        sub_itr_list = self.itr_list[::sub_sampling_factor]

        ax.plot(sub_distance, sub_radius, color='black')
        ax.plot(sub_distance, -sub_radius, color='black')

        ax.fill_between(
            sub_distance,
            sub_radius,
            -sub_radius,
            color='lightblue',
            alpha=0.8
        )

        ax.axvline(sub_distance[0], linestyle='--', color='red')

        if dark_background:
            style = plt.style.context("dark_background")
            ax.tick_params(colors='white', direction='out')
        else:
            style = plt.stryle.context('default')

        with style:
            def animate(i):
                ax.lines[-1].remove()
                line0 = ax.set_title(f'[slice: {i} - ITR: {sub_itr_list[i]:.3f}]', color='white')
                line = ax.axvline(sub_distance[i], linestyle='--', color='red')
                return line0, line

            animation = FuncAnimation(
                fig=figure,
                func=animate,
                blit=True,
                repeat=True,
                frames=number_of_frames
            )

            animation.save(
                output_directory,
                dpi=dpi,
                writer=PillowWriter(fps=fps)
            )

# -
