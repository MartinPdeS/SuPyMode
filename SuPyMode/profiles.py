#!/usr/bin/env python
# -*- coding: utf-8 -*-

import numpy
from scipy.interpolate import interp1d
import matplotlib.pyplot as plt

from SuPyMode.tools import plot_style
from MPSPlots.render2D import SceneList, Axis
from matplotlib.animation import FuncAnimation, PillowWriter

from dataclasses import dataclass


@dataclass
class AlphaProfile():
    initial_radius: float = 1
    """ Initial radius of the taper structure """
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
        self.taper_segment = []
        self.z_segment = []
        self.radius_segments = []
        self.heating_length_segment = []

    def symmetrize_array(self, array: numpy.ndarray) -> numpy.ndarray:
        return numpy.r_[array, array[::-1]]

    def symmetrize_distance(self, distance: numpy.ndarray) -> numpy.ndarray:
        dz = abs(distance[0] - distance[1])
        return numpy.arange(2 * distance.size) * dz

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

    def add_end_of_taper_segment(self, *, length: float = None, n_point: int = 100) -> None:
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

    def add_constant_custom_section(self, *,
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
            z + start_z,
            radius,
            bounds_error=False,
            fill_value=0
        )

        self.taper_segment.append(interpolation)
        self.z_segment.append(length + start_z)

    def evaluate_adiabatic_factor(self, itr: numpy.ndarray) -> numpy.ndarray:
        interpolation = interp1d(
            x=self._itr_list,
            y=self._adiabatic_factor,
            bounds_error=False,
            fill_value=numpy.nan
        )

        return interpolation(itr)

    def evaluate_distance_vs_itr(self, distance: numpy.ndarray) -> numpy.ndarray:
        interpolation = interp1d(
            x=self._itr_list,
            y=self._distance_array,
            bounds_error=True,
        )

        return interpolation(distance)

    def get_radius_from_segment(self, *,
            alpha: float,
            initial_heating_length: float,
            stretching_length: float,
            initial_radius: float,
            distance: numpy.ndarray) -> tuple:
        """
        Gets the radius as a fonction of the distance for a specific segment.,

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

    def assert_conditions(self, *,
                          alpha: float,
                          stretching_length: float,
                          initial_heating_length: float) -> None:

        assert initial_heating_length > 0, "The initial heat lenght initial_heating_length cannot be negative!"

        if alpha < 0:
            assert stretching_length < initial_heating_length / abs(alpha), "Condition: x0 < initial_heating_length / |alpha| is not respected! see Birks article in the references!"

    def add_custom_segment(self, distance: numpy.ndarray, radius: numpy.ndarray) -> None:
        end_of_segment = distance[-1]
        final_radius_of_segment = radius[-1]

        interpolation = interp1d(
            distance,
            radius,
            bounds_error=False,
            fill_value=0
        )

        self.taper_segment.append(interpolation)
        self.z_segment.append(end_of_segment)
        self.radius_segments.append(final_radius_of_segment)

    def add_taper_custom_segment(self, *,
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
            distance + start_z,
            radius,
            bounds_error=False,
            fill_value=0
        )

        self.taper_segment.append(interpolation)
        self.z_segment.append(end_of_segment)
        self.radius_segments.append(final_radius)
        self.heating_length_segment.append(final_heating_length)

    @property
    def distance(self) -> numpy.ndarray:
        if self.symmetric:
            return self._symmetric_distance_array
        else:
            return self._distance_array

    @property
    def radius(self) -> numpy.ndarray:
        """
        Returns the array of radius [vs z-distance] for the taper structure

        :returns:   The ITR array
        :rtype:     numpy.ndarray
        """
        if self.symmetric:
            return self._symmetric_radius_array
        else:
            return self._radius_array

    @property
    def itr_list(self) -> numpy.ndarray:
        """
        Returns the array of ITR value [vs z-distance] for the taper structure

        :returns:   The ITR array
        :rtype:     numpy.ndarray
        """
        if self.symmetric:
            return self._symmetric_itr_list
        else:
            return self._itr_list

    @property
    def adiabatic(self) -> numpy.ndarray:
        """
        Returns the array of adiabatc factor [vs ITR] for the taper structure

        :returns:   The adiabatic factor
        :rtype:     numpy.ndarray
        """
        if self.symmetric:
            return self._symmetric_adiabatic_factor_array
        else:
            return self._adiabatic_factor

    @property
    def taper_angle(self) -> numpy.ndarray:
        """
        Returns the array of taper angle for the taper structure

        :returns:   The taper angle array
        :rtype:     numpy.ndarray
        """
        if self.symmetric:
            return self._symmetric_taper_angle_array
        else:
            return self._taper_angle_array

    @property
    def smallest_itr(self) -> float:
        """
        Returns the smallest itr of the taper structure

        :returns:   Smallest itr value
        :rtype:     float
        """
        return numpy.min(self.itr_list)

    @property
    def last_z(self) -> float:
        """
        Retunrs the last, or equavalently the largest propagation distance computed

        :returns:   The z-distance
        :rtype:     float
        """
        if len(self.z_segment) == 0:
            return 0
        else:
            return self.z_segment[-1]

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

    def get_radius_from_segment_from_interpolation(self, z: numpy.ndarray):
        radius = numpy.zeros(z.size)

        for interpolation in self.taper_segment:
            radius += interpolation(z)

        return radius

    def add_taper_segment(self, *,
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
        if initial_radius is None:
            initial_radius = self.last_radius

        return self.add_taper_custom_segment(
            alpha=alpha,
            initial_heating_length=initial_heating_length,
            initial_radius=initial_radius,
            stretching_length=stretching_length,
            start_z=self.last_z,
            n_point=n_point
        )

    def initialize(self, n_point: int = 400) -> None:
        """
        Initialize all the computation including z, radius, distance, length, adiabatic criterion

        :param      n_point:  The number of point of the z-linespace to evaluate all the parameters.
        :type       n_point:  int
        """
        self._distance_array = numpy.linspace(0, self.last_z, n_point)
        self._radius_array = self.get_radius_from_segment_from_interpolation(self._distance_array)
        self._itr_list = self._radius_array / self.initial_radius
        self._adiabatic_factor = self.get_adiabatic_factor()
        self._taper_angle_array = self.get_taper_angle()

        self.symmetrize_parameters()
        self.length = self.distance[-1]

        self.master_interpolation_z_to_itr = interp1d(
            self.distance,
            self.itr_list,
            bounds_error=False,
            fill_value=0
        )

    def symmetrize_parameters(self):
        self._symmetric_distance_array = self.symmetrize_distance(self._distance_array)
        self._symmetric_radius_array = self.symmetrize_array(self._radius_array)
        self._symmetric_itr_list = self.symmetrize_array(self._itr_list)
        self._symmetric_adiabatic_factor_array = self.symmetrize_array(self._adiabatic_factor)
        self._symmetric_taper_angle_array = self.symmetrize_array(self._taper_angle_array)

    def get_adiabatic_factor(self) -> numpy.ndarray:
        r"""
        Compute the adiabatic factor defined as:
        .. math::
          f_c = \frac{1}{\rho} \frac{d \rho}{d z}

        :returns:   The amplitudes as a function of the distance in the coupler
        :rtype:     numpy.ndarray
        """
        dz = numpy.gradient(self._distance_array, axis=0, edge_order=2)

        ditr = numpy.gradient(numpy.log(self._radius_array), axis=0, edge_order=2)

        return abs(ditr / dz)

    def get_taper_angle(self) -> numpy.ndarray:
        r"""
        From Tapered single-mode fibres and devices. Part 1: Adiabaticity criteria.
        Compute the adiabatic factor defined as:
        .. math::
          f_c = \frac{d \rho}{d z} = \Omega

        :returns:   The amplitudes as a function of the distance in the coupler
        :rtype:     numpy.ndarray
        """
        d_z = numpy.gradient(self._distance_array, axis=0, edge_order=2)

        d_rho = numpy.gradient(self._radius_array, axis=0, edge_order=2)

        return abs(d_rho / d_z)

    def get_length_scale(self, core_radius: float) -> numpy.ndarray:
        r"""
        Compute the adiabatic factor defined as:
        .. math::
          f_c = \frac{\rho}{\Omega}

        :returns:   The amplitudes as a function of the distance in the coupler
        :rtype:     numpy.ndarray
        """
        return core_radius / self.get_taper_angle()

    def _render_itr_vs_z_on_ax_(self, ax: Axis, **artist_kwargs) -> None:
        """
        Add plot onto axis, the plots is ITR vs Z-distance

        :param      ax:   The axis on which to add the plot
        :type       ax:   Axis
        """
        ax.set_style(**plot_style.z_profile)

        ax.add_line(
            x=self.distance,
            y=self.radius / self.initial_radius,
            label=self.label,
            **artist_kwargs
        )

    def _render_taper_angle_vs_z_on_ax_(self, ax: Axis, **artist_kwargs) -> None:
        """
        Add plot onto axis, the plots is ITR vs Z-distance

        :param      ax:   The axis on which to add the plot
        :type       ax:   Axis
        """
        ax.set_style(**plot_style.taper_angle)

        ax.add_line(
            x=self.distance,
            y=self.taper_angle,
            label=self.label,
            **artist_kwargs
        )

    def _render_taper_length_scale_vs_itr_on_ax_(self, ax: Axis, core_radius: float, **artist_kwargs) -> None:
        """
        Add ITR vs Z-distance plot to axis

        :param      ax:   The axis on which to add the plot
        :type       ax:   Axis
        """
        ax.set_style(**plot_style.taper_angle)

        length_scale = self.get_length_scale(core_radius=core_radius)

        ax.add_line(
            x=self.itr_list,
            y=length_scale,
            label='Coupler length-scale',
            **artist_kwargs
        )

    def _render_adiabatic_factor_vs_z_on_ax_(self, ax: Axis, **artist_kwargs) -> None:
        """
        Add plot onto axis, the plots is adiabatic criterion vs Z-distance

        :param      ax:   The axis on which to add the plot
        :type       ax:   Axis
        """
        ax.y_scale = 'log'
        ax.y_label = 'Adiabatic criterion'
        ax.x_label = 'z-distance'

        ax.add_line(
            x=self._distance_array,
            y=self._adiabatic_factor,
            label=self.label,
            **artist_kwargs
        )

    def _render_adiabatic_factor_vs_itr_on_ax_(self, ax: Axis, **artist_kwargs) -> None:
        """
        Add adiabatic criterion vs ITR plot to axis

        :param      ax:   The axis on which to add the plot
        :type       ax:   Axis
        """
        ax.set_style(**plot_style.adiabatic)

        ax.add_line(
            x=self.itr_list,
            y=self.adiabatic,
            label=self.label,
            **artist_kwargs
        )

    def plot(self,
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
        itr = self._symmetric_radius_array / self.initial_radius

        figure = SceneList(
            title=f'Minimum ITR: {itr.min():.4f}',
            ax_orientation='vertical',
            unit_size=(8, 2)
        )

        if show_radius:
            ax = figure.append_ax(y_limits=[0, None])
            self._render_itr_vs_z_on_ax_(ax)

        if show_taper_angle:
            ax = figure.append_ax(y_limits=[0, None])
            self._render_taper_angle_vs_z_on_ax_(ax)

        if show_adiabatic:
            ax = figure.append_ax()
            self._render_adiabatic_factor_vs_itr_on_ax_(ax)

        figure.annotate_axis(position=(-0.15, 1.15))
        return figure

    def generate_propagation_gif(self,
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


if __name__ == '__main__':
    profile = AlphaProfile(
        initial_radius=62.5e-6,
        symmetric=False,
        label='test profile',
        line_color='red'
    )

    profile.add_taper_segment(
        alpha=0,
        initial_heating_length=5e-3,
        stretching_length=20e-3,
        n_point=100
    )

    profile.initialize()

    profile.plot().show()

# -
