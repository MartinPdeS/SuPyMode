#!/usr/bin/env python
# -*- coding: utf-8 -*-

# Built-in imports
import numpy
from dataclasses import dataclass

# Third-party imports
from MPSPlots.render2D import SceneList

# Local imports
from PyFiberModes import FiberFactory
from FiberFusing.fiber import catalogue as fiber_catalogue
from SuPyMode.tools.analytics.supermode import Supermode


@dataclass
class DataSet:
    x: numpy.ndarray
    """ Array reprenting the x-axis """
    y: numpy.ndarray
    """ Array reprenting the y-axis """
    x_label: str = ''
    """ String for the x-axis of the plot"""
    y_label: str = ''
    """ String for the y-axis of the plot"""
    title: str = ''
    y_scale: str = 'linear'

    def plot(self, **kwargs):
        figure = SceneList(title='', unit_size=(12, 4))

        ax = figure.append_ax(
            x_label=self.x_label,
            y_label=self.y_label,
            y_scale=self.y_scale
        )

        ax.add_line(x=self.x, y=self.y, label=self.label)

        return figure


@dataclass
class Superset():
    wavelength: float
    """ Wavelenght at which compute the propagation modes """

    def __post_init__(self):
        self.wavenumber = 2 * numpy.pi / self.wavelength

    def get_mode(self, mode_number: str, wavelength: float, itr: float) -> Supermode:
        fiber = self.get_smf28_model(itr=itr)

        return Supermode(
            fiber=fiber,
            mode_number=mode_number,
            wavelength=self.wavelength
        )

    def get_smf28_model(self, itr: float = 1.0):
        return self.get_custom_fiber_model(
            custom_fiber_class=fiber_catalogue.SMF28,
            itr=itr
        )

    def get_custom_fiber_model(self, custom_fiber_class, itr: float = 1):
        factory = FiberFactory()
        factory.addLayer()

        custom_fiber_instance = custom_fiber_class(wavelength=self.wavelength)

        custom_fiber_instance.scale(factor=itr)

        for structure in custom_fiber_instance.structure_list[::-1]:
            factory.addLayer(
                name=structure.name,
                radius=structure.radius,
                index=structure.index
            )

        fiber = factory[0]

        fiber.supymode_fiber = custom_fiber_instance

        fiber.radius_boundary = numpy.max(custom_fiber_instance.boundaries)

        return fiber

    def get_normalized_mode_field(self,
            mode,
            r_space: numpy.ndarray,
            normalization: str = 'max') -> numpy.ndarray:
        """
        Gets the normalized mode field.

        :param      mode:           The mode
        :type       mode:           { type_description }
        :param      r_space:        The r space
        :type       r_space:        numpy.ndarray
        :param      normalization:  The normalization
        :type       normalization:  str

        :returns:   The normalized mode field.
        :rtype:     numpy.ndarray
        """
        if (mode.l, mode.m) == (0, 1):
            field = abs(mode.evaluate_field_at_r(r_space))
        else:
            field = (mode.evaluate_field_at_r(r_space))

        match normalization.lower():
            case 'max':
                norm = abs(field).max()
            case 'l2':
                dr = abs(r_space[1] - r_space[0])
                norm = numpy.trapz(numpy.square(field) * r_space, dx=dr, axis=0) * 2 * numpy.pi
                norm = numpy.sqrt(norm)
            case 'center':
                idx_center = numpy.argmin(abs(r_space))
                center_value = field[idx_center]
                norm = center_value
            case 'cmt':
                dr = abs(r_space[1] - r_space[0])
                norm = numpy.trapz(numpy.square(field) * r_space, dx=dr, axis=0)
                norm = numpy.sqrt(norm)
            case 'scalar_coupling':  # equation page 229: Bures
                dr = abs(r_space[1] - r_space[0])
                norm = 0.5 * numpy.trapz(numpy.square(field) * r_space, dx=dr, axis=0) * 2 * numpy.pi
                norm = numpy.sqrt(norm)

        return field / norm

    def get_mode_field_vs_r(self,
            itr: float,
            mode_number: Supermode,
            resolution: int = 100,
            normalization: str = 'max') -> float:
        """
        Gets the mode field as a function of r.

        :param      itr:            The itr
        :type       itr:            float
        :param      mode_number:    The mode number
        :type       mode_number:    Supermode
        :param      resolution:     The resolution
        :type       resolution:     int
        :param      normalization:  The normalization
        :type       normalization:  str

        :returns:   The mode field vs r.
        :rtype:     float
        """
        fiber = self.get_smf28_model(itr=itr)

        bound = fiber.radius_boundary * 1.2

        r_space = numpy.linspace(0, bound, resolution)

        mode = Supermode(
            fiber=fiber,
            mode_number=mode_number,
            wavelength=self.wavelength
        )

        field_r_mode = self.get_normalized_mode_field(
            mode=mode,
            r_space=r_space,
            normalization=normalization
        )

        data_set = DataSet(
            x=r_space,
            y=field_r_mode,
            y_label='Field amplitude',
            x_label='Radial distance',
            title=mode_number
        )

        data_set.fiber = fiber

        return data_set

    def get_adiabatic_vs_itr(self,
            mode_number_0: str,
            mode_number_1: str,
            itr_list: list,
            resolution: int,
            debug_mode: bool = True) -> DataSet:
        """
        Gets the adiabatic criterion figure vs itr.

        :param      mode_number_0:  The mode number 0
        :type       mode_number_0:  str
        :param      mode_number_1:  The mode number 1
        :type       mode_number_1:  str
        :param      itr_list:       The itr list
        :type       itr_list:       list
        :param      resolution:     The resolution
        :type       resolution:     int

        :returns:   The adiabatic vs itr.
        :rtype:     DataSet
        """
        array = numpy.zeros(itr_list.size)

        for idx, itr in enumerate(itr_list):
            if debug_mode:
                print(f'{itr:.3f = }\t{idx = }', end='\r')

            fiber = self.get_smf28_model(itr=itr)

            mode_0 = Supermode(
                fiber=fiber,
                mode_number=mode_number_0,
                wavelength=self.wavelength
            )

            mode_1 = Supermode(
                fiber=fiber,
                mode_number=mode_number_1,
                wavelength=self.wavelength
            )

            coupling = self.get_normalized_coupling(
                fiber=fiber,
                mode_0=mode_0,
                mode_1=mode_1,
                resolution=resolution
            )

            delta_beta = abs(mode_0.beta - mode_1.beta)

            array[idx] = delta_beta / coupling

        data_set = DataSet(
            x=itr_list,
            y=array,
            y_label='adiabatic criterion',
            x_label='ITR',
            y_scale='log',
            title=mode_number_0 + ':' + mode_number_1
        )

        return data_set

    def get_normalized_coupling(self,
            fiber,
            mode_0: Supermode,
            mode_1: Supermode,
            resolution: int = 200) -> float:
        """
        Gets the normalized coupling between two supermodes.

        :param      fiber:       The fiber
        :type       fiber:       FiberFactory
        :param      mode_0:      The msuperode 0
        :type       mode_0:      Supermode
        :param      mode_1:      The supermode 1
        :type       mode_1:      Supermode
        :param      resolution:  The resolution
        :type       resolution:  int

        :returns:   The normalized coupling.
        :rtype:     float
        """
        bound = fiber.radius_boundary * 1.2

        r_space = numpy.linspace(0, bound, resolution)

        field_mode_0_vs_r = self.get_normalized_mode_field(
            mode=mode_0,
            r_space=r_space,
            normalization='l2'
        )

        field_mode_1_vs_r = self.get_normalized_mode_field(
            mode=mode_1,
            r_space=r_space,
            normalization='l2'
        )

        fields_term = 0
        for structure in fiber.supymode_fiber.structure_list:
            current_layer = structure
            if structure.name == 'air':
                previous_layer = current_layer
                continue

            idx_radius = numpy.argmin(abs(r_space - current_layer.radius))
            field_mode_0_at_radius = field_mode_0_vs_r[idx_radius]
            field_mode_1_at_radius = field_mode_1_vs_r[idx_radius]

            fields_term += structure.radius**2 * (current_layer.index**2 - previous_layer.index**2) * (field_mode_0_at_radius * field_mode_1_at_radius)

            previous_layer = current_layer

        # Equation 7.39 Jacques Bures
        term0 = 0.5
        term0 *= self.wavenumber ** 2 / numpy.sqrt(mode_0.beta * mode_1.beta)
        term0 *= 1 / abs(mode_0.beta - mode_1.beta)

        coupling = abs(term0 * fields_term)

        return coupling

    def get_overlap_integral(self,
            fiber,
            mode_0: Supermode,
            mode_1: Supermode,
            resolution: int = 200,
            normalization: str = 'l2') -> float:

        bound = fiber.radius_boundary * 1.2

        r_space = numpy.linspace(0, bound, resolution)

        field_mode_0_vs_r = self.get_normalized_mode_field(
            mode=mode_0,
            r_space=r_space,
            normalization=normalization
        )

        field_mode_1_vs_r = self.get_normalized_mode_field(
            mode=mode_1,
            r_space=r_space,
            normalization=normalization
        )

        dr = abs(r_space[1] - r_space[0])

        overlap_integral = numpy.trapz(field_mode_0_vs_r * field_mode_1_vs_r * r_space, dx=dr, axis=0) * 2 * numpy.pi

        return overlap_integral

    def get_overlap_integral_vs_itr(self,
            mode_number_0: str,
            mode_number_1: str,
            itr_list: numpy.ndarray,
            resolution: int,
            debug_mode: bool = True) -> DataSet:
        """
        Gets the overlap integral of two supermodes vs itr.

        :param      mode_number_0:  The mode number 0
        :type       mode_number_0:  str
        :param      mode_number_1:  The mode number 1
        :type       mode_number_1:  str
        :param      itr_list:       The itr list
        :type       itr_list:       numpy.ndarray
        :param      resolution:     The resolution
        :type       resolution:     int

        :returns:   The overlap integral vs itr.
        :rtype:     DataSet
        """
        array = numpy.zeros(itr_list.size)

        for idx, itr in enumerate(itr_list):
            if debug_mode:
                print(f'{itr:.3f = }\t{idx = }', end='\r')
            fiber = self.get_smf28_model(itr=itr)

            mode_0 = Supermode(
                fiber=fiber,
                mode_number=mode_number_0,
                wavelength=self.wavelength
            )

            mode_1 = Supermode(
                fiber=fiber,
                mode_number=mode_number_1,
                wavelength=self.wavelength
            )

            overlap_integral = self.get_overlap_integral(
                fiber=fiber,
                mode_0=mode_0,
                mode_1=mode_1,
                resolution=resolution
            )

            array[idx] = overlap_integral

        data_set = DataSet(
            x=itr_list,
            y=array,
            y_label='Normalized coupling',
            x_label='ITR',
            title=mode_number_0 + ':' + mode_number_1
        )

        return data_set

    def get_beta(self, mode_number: str, itr: float) -> float:
        """
        Gets the propagation constant beta.

        :param      mode_number:  The mode number
        :type       mode_number:  str
        :param      itr:          The itr
        :type       itr:          float

        :returns:   The beta.
        :rtype:     float
        """
        fiber = self.get_smf28_model(itr=itr)

        mode = Supermode(
            fiber=fiber,
            mode_number=mode_number,
            wavelength=self.wavelength
        )

        return mode.beta

    def get_beta_vs_itr(self,
            mode_number: str,
            itr_list: numpy.ndarray,
            debug_mode: bool = True) -> DataSet:
        """
        Gets the beta vs itr.

        :param      mode_number:  The mode number
        :type       mode_number:  str
        :param      itr_list:     The itr list where beta is evaluated
        :type       itr_list:     numpy.ndarray

        :returns:   The beta vs itr.
        :rtype:     DataSet
        """
        array = numpy.zeros(itr_list.size)

        for idx, itr in enumerate(itr_list):
            if debug_mode:
                print(f'\t{idx = }\t{itr = :.3f}', end='\r')
            fiber = self.get_smf28_model(itr=itr)

            mode = Supermode(
                fiber=fiber,
                mode_number=mode_number,
                wavelength=self.wavelength
            )

            try:
                propagation_constant = mode.beta
            except ValueError:
                propagation_constant = numpy.nan

            array[idx] = propagation_constant

        data_set = DataSet(
            x=itr_list,
            y=array,
            y_label='Propagation constant',
            x_label='ITR',
            title=mode_number
        )

        return data_set

    def get_normalized_coupling_vs_itr(self,
            mode_number_0: str,
            mode_number_1: str,
            itr_list: numpy.ndarray,
            resolution: int,
            debug_mode: bool = True) -> DataSet:
        """
        Gets the normalized coupling vs itr.

        :param      mode_number_0:  The mode number 0
        :type       mode_number_0:  str
        :param      mode_number_1:  The mode number 1
        :type       mode_number_1:  str
        :param      itr_list:       The itr list
        :type       itr_list:       { type_description }
        :param      resolution:     The resolution
        :type       resolution:     int

        :returns:   The normalized coupling vs itr.
        :rtype:     DataSet
        """
        array = numpy.zeros(itr_list.size)

        for idx, itr in enumerate(itr_list):
            fiber = self.get_smf28_model(itr=itr)

            mode_0 = Supermode(
                fiber=fiber,
                mode_number=mode_number_0,
                wavelength=self.wavelength
            )

            mode_1 = Supermode(
                fiber=fiber,
                mode_number=mode_number_1,
                wavelength=self.wavelength
            )

            if mode_0.l != mode_1.l:
                coupling = 0
            else:
                try:
                    coupling = self.get_normalized_coupling(
                        fiber=fiber,
                        mode_0=mode_0,
                        mode_1=mode_1,
                        resolution=resolution
                    )
                except ValueError:
                    coupling = numpy.nan

            array[idx] = coupling

        return DataSet(
            x=itr_list,
            y=array,
            y_label='Normalized coupling',
            x_label='ITR',
            title=mode_number_0 + ':' + mode_number_1
        )


# -
