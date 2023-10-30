#!/usr/bin/env python
# -*- coding: utf-8 -*-

# Built-in imports
import numpy
from dataclasses import dataclass
from scipy.constants import epsilon_0 as e0, mu_0

# Third-party imports
from MPSPlots.render2D import SceneList

# Local imports
from PyFiberModes import Wavelength, FiberFactory, Mode, field as FieldClass
from FiberFusing.fiber import catalogue as fiber_catalogue


def get_fibermodes_custom_fiber(custom_fiber, wavelength: float):
    factory = FiberFactory()
    factory.addLayer()

    for structure in custom_fiber.structure_list[::-1]:
        factory.addLayer(
            name=structure.name,
            radius=structure.radius,
            index=structure.index
        )

    fiber = factory[0]

    fiber.supymode_fiber = custom_fiber

    fiber.radius_boundary = numpy.max(custom_fiber.boundaries)

    return fiber


def get_fibermodes_smf28_taper(wavelength: float, itr: float):
    custom_fiber = fiber_catalogue.SMF28(wavelength=wavelength)
    custom_fiber.scale(factor=itr)

    return get_fibermodes_custom_fiber(
        custom_fiber=custom_fiber,
        wavelength=wavelength,
    )


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

        ax = figure.append_ax(x_label=self.x_label, y_label=self.y_label, y_scale=self.y_scale)

        ax.add_line(x=self.x, y=self.y, label=self.label)

        return figure


@dataclass
class FibermodeSupermode():
    mode_number: str
    """ Mode number [starts with LP] """
    fiber: FiberFactory
    """ Fiber type to which compute the mode """
    wavelength: float
    """ Wavelenght to which compute the mode """

    def __post_init__(self):
        l, m = self.mode_number[2:]
        self.l, self.m = int(l), int(m)
        self.mode = Mode('LP', self.l, self.m)
        self.wavenumber = 2 * numpy.pi / self.wavelength

    @property
    def beta(self):
        """ Returns propgation constant of the supermode """
        return self.wavenumber * self.neff

    @property
    def neff(self) -> float:
        """ Returns effective refractive index """
        return self.fiber.neff(self.mode, self.wavelength)

    @property
    def norm_factor(self) -> float:
        """ Returns the norm factor of the supermode """
        factor = 0.5 * numpy.sqrt(e0 / mu_0)

        return factor * (self.beta / self.wavenumber)

    def get_fiber_instance(self):
        return get_fibermodes_custom_fiber

    def get_field(self, bound: float, resolution: int = 200) -> numpy.ndarray:
        fields = FieldClass.Field(
            fiber=self.fiber,
            mode=self.mode,
            wl=self.wavelength,
            r=bound,
            np=resolution
        )

        return fields.Ex()

    def get_norm2_field(self, field: numpy.ndarray) -> float:
        """ Returns the L2 norm of the supermode field """
        field = numpy.square(field)
        sum_field = field.sum(axis=0).sum(axis=0)
        return numpy.sqrt(sum_field)

    def get_l2_normalized_field(self, bound: float, resolution: int = 200) -> numpy.ndarray:
        """ Returns a L2 normalized field array """
        field_object = FieldClass.Field(
            fiber=self.fiber,
            mode=self.mode,
            wl=self.wavelength,
            r=bound,
            np=resolution
        )

        field_array = field_object.Ex()

        norm_l2 = self.get_norm2_field(field_array)

        normalized_field = field_array / numpy.sqrt(norm_l2)

        return normalized_field

    def evaluate_field_at_r(self, r_space: numpy.ndarray) -> numpy.ndarray:
        """ Returns array corresponding to the supermode field evaluated at the r_space position """
        wavelength = Wavelength(self.wavelength)

        array = numpy.zeros(r_space.size)

        for idx, r in enumerate(r_space):
            er, hr = self.fiber._rfield(self.mode, wavelength, r)
            array[idx] = er[0]

        return array.squeeze()


class FibermodesSuperset():
    def __init__(self, wavelength: float):
        self.wavelength = wavelength

    @property
    def wavenumber(self) -> float:
        return 2 * numpy.pi / self.wavelength

    def get_mode(self, mode_number: str, wavelength: float, itr: float) -> FibermodeSupermode:
        fiber = get_fibermodes_smf28_taper(wavelength=wavelength, itr=itr)

        return FibermodeSupermode(
            fiber=fiber,
            mode_number=mode_number,
            wavelength=self.wavelength
        )

    def get_normalized_mode_field(self, mode, r_space, normalization: str = 'max'):
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
            mode_number: FibermodeSupermode,
            resolution: int = 100,
            normalization: str = 'max') -> float:
        """
        Gets the mode field as a function of r.

        :param      itr:            The itr
        :type       itr:            float
        :param      mode_number:    The mode number
        :type       mode_number:    FibermodeSupermode
        :param      resolution:     The resolution
        :type       resolution:     int
        :param      normalization:  The normalization
        :type       normalization:  str

        :returns:   The mode field vs r.
        :rtype:     float
        """
        fiber = get_fibermodes_smf28_taper(wavelength=self.wavelength, itr=itr)

        bound = fiber.radius_boundary * 1.2

        r_space = numpy.linspace(0, bound, resolution)

        mode = FibermodeSupermode(
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

            fiber = get_fibermodes_smf28_taper(wavelength=self.wavelength, itr=itr)

            mode_0 = FibermodeSupermode(
                fiber=fiber,
                mode_number=mode_number_0,
                wavelength=self.wavelength
            )

            mode_1 = FibermodeSupermode(
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
            mode_0: FibermodeSupermode,
            mode_1: FibermodeSupermode,
            resolution: int = 200) -> float:
        """
        Gets the normalized coupling between two supermodes.

        :param      fiber:       The fiber
        :type       fiber:       FiberFactory
        :param      mode_0:      The msuperode 0
        :type       mode_0:      FibermodeSupermode
        :param      mode_1:      The supermode 1
        :type       mode_1:      FibermodeSupermode
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
            mode_0: FibermodeSupermode,
            mode_1: FibermodeSupermode,
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
            fiber = get_fibermodes_smf28_taper(wavelength=self.wavelength, itr=itr)

            mode_0 = FibermodeSupermode(
                fiber=fiber,
                mode_number=mode_number_0,
                wavelength=self.wavelength
            )

            mode_1 = FibermodeSupermode(
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
        fiber = get_fibermodes_smf28_taper(wavelength=self.wavelength, itr=itr)

        mode = FibermodeSupermode(
            fiber=fiber,
            mode_number=mode_number,
            wavelength=self.wavelength
        )

        return mode.beta

    def get_beta_vs_itr(self, mode_number: str, itr_list: numpy.ndarray, debug_mode: bool = True) -> DataSet:
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
            fiber = get_fibermodes_smf28_taper(wavelength=self.wavelength, itr=itr)

            mode = FibermodeSupermode(
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
            fiber = get_fibermodes_smf28_taper(wavelength=self.wavelength, itr=itr)

            mode_0 = FibermodeSupermode(
                fiber=fiber,
                mode_number=mode_number_0,
                wavelength=self.wavelength
            )

            mode_1 = FibermodeSupermode(
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


class FiberModeSolver():
    def __init__(self, wavelength: float):
        self.wavelength = wavelength
        self.solver = FibermodesSuperset(wavelength=wavelength)

    def get_field_distribution(self,
            itr: float,
            resolution: int,
            mode_numbers: list,
            normalization: str = 'max') -> list:
        """
        Gets the fibermode field distribution.

        :param      wavelength:     The wavelength
        :type       wavelength:     float
        :param      itr:            The itr
        :type       itr:            float
        :param      resolution:     The resolution
        :type       resolution:     int
        :param      mode_numbers:   The mode numbers
        :type       mode_numbers:   list
        :param      normalization:  The normalization
        :type       normalization:  str

        :returns:   The fibermode field distributions.
        :rtype:     list[DataSet]
        """
        data_sets = []
        for mode_number in mode_numbers:
            print(f'{mode_number = }', end='\r')
            data_set = self.solver.get_mode_field_vs_r(
                mode_number=mode_number,
                itr=itr,
                resolution=resolution,
                normalization=normalization
            )

            data_sets.append(data_set)

        return data_sets

    def plot_field_distribution(self,
            itr: float,
            resolution: int,
            mode_numbers: list,
            normalization: str = 'max',
            add_interfaces: bool = True) -> SceneList:
        """
        Gets the fibermode field distribution.

        :param      wavelength:     The wavelength
        :type       wavelength:     float
        :param      itr:            The itr
        :type       itr:            float
        :param      resolution:     The resolution
        :type       resolution:     int
        :param      mode_numbers:   The mode numbers
        :type       mode_numbers:   list
        :param      normalization:  The normalization
        :type       normalization:  str

        :returns:   The figure.
        :rtype:     SceneList
        """
        figure = SceneList(title='Fibermodes field distribution')
        ax = figure.append_ax(
            show_legend=True,
            line_width=2,
            x_label='Radial distance',
            y_label='Normalized amplitude'
        )

        data_sets = self.get_field_distribution(
            mode_numbers=mode_numbers,
            resolution=resolution,
            normalization=normalization,
            itr=itr
        )

        for idx, (mode, fibermode_data_set) in enumerate(zip(mode_numbers, data_sets)):
            color = f'C{idx}'
            ax.add_line(
                x=fibermode_data_set.x,
                y=fibermode_data_set.y,
                color=color,
                label=fibermode_data_set.title,
            )

        if add_interfaces:
            for structure in data_sets[0].fiber.supymode_fiber.structure_list:
                if structure.name == 'air':
                    continue

                ax.add_vertical_line(
                    x=structure.radius,
                    y_max=ax.get_y_max(),
                    y_min=ax.get_y_min(),
                    color='black',
                    line_style='--'
                )

        return figure

    def get_normalized_coupling(self,
            itr_list: numpy.ndarray,
            resolution: int,
            mode_couples: list,
            debug_mode: bool = True) -> list:
        """
        Gets the fibermode coupling.

        :param      itr_list:      The itr list
        :type       itr_list:      list
        :param      resolution:    The resolution for mode field
        :type       resolution:    int
        :param      mode_couples:  The mode couples to compute
        :type       mode_couples:  list

        :returns:   The fibermode coupling.
        :rtype:     list
        """
        data_sets = []
        for mode_couple in mode_couples:
            if debug_mode:
                print(f'{mode_couple = }', end='\r')

            coupling_data_set = self.solver.get_normalized_coupling_vs_itr(
                mode_number_0=mode_couple[0],
                mode_number_1=mode_couple[1],
                itr_list=itr_list,
                resolution=resolution,
                debug_mode=debug_mode
            )

            data_sets.append(coupling_data_set)

        return data_sets

    def plot_normalized_coupling(self,
            itr_list: numpy.ndarray,
            resolution: int,
            mode_couples: list) -> SceneList:
        """
        Plot the fibermode coupling.

        :param      itr_list:      The itr list
        :type       itr_list:      list
        :param      resolution:    The resolution for mode field
        :type       resolution:    int
        :param      mode_couples:  The mode couples to compute
        :type       mode_couples:  list

        :returns:   The fibermode coupling.
        :rtype:     SceneList
        """
        figure = SceneList(title='Fibermodes normalized coupling')
        ax = figure.append_ax(
            show_legend=True,
            line_width=2,
            x_label='Inverse taper ratio [ITR]',
            y_label='Normalized coupling'
        )

        data_sets = self.get_normalized_coupling(
            resolution=resolution,
            mode_couples=mode_couples,
            itr_list=itr_list
        )

        for idx, (mode_couple, fibermode_data_set) in enumerate(zip(mode_couples, data_sets)):
            color = f'C{idx}'
            ax.add_line(
                x=fibermode_data_set.x,
                y=fibermode_data_set.y,
                color=color,
                label=fibermode_data_set.title,
            )

        return figure

    def get_beta_vs_itr(self,
            mode_numbers: list,
            itr_list: numpy.ndarray,
            debug_mode: bool = True) -> list:
        """
        Gets the propagation constant [beta] vs itr.

        :param      mode_numbers:  The mode numbers
        :type       mode_numbers:  list
        :param      itr_list:      The itr list
        :type       itr_list:      numpy.ndarray

        :returns:   The beta vs itr.
        :rtype:     list
        """
        data_sets = []
        for mode_number in mode_numbers:
            if debug_mode:
                print(f'{mode_number = }', end='\r')
            data_set = self.solver.get_beta_vs_itr(
                mode_number=mode_number,
                itr_list=itr_list,
            )

            data_sets.append(data_set)

        return data_sets

    def plot_beta_vs_itr(self, mode_numbers: list, itr_list: numpy.ndarray) -> SceneList:
        """
        Plot the propagation constant [beta] vs itr.

        :param      mode_numbers:  The mode numbers
        :type       mode_numbers:  list
        :param      itr_list:      The itr list
        :type       itr_list:      numpy.ndarray

        :returns:   The figure
        :rtype:     SceneList
        """
        figure = SceneList(title='Fibermodes propagation constant')

        ax = figure.append_ax(
            show_legend=True,
            line_width=2,
            y_label=r'Propagation constant [\beta]',
            x_label='Inverse taper ratio [ITR]'
        )

        data_sets = self.get_beta_vs_itr(
            mode_numbers=mode_numbers,
            itr_list=itr_list
        )

        for idx, (mode_couple, data_set) in enumerate(zip(mode_numbers, data_sets)):
            color = f'C{idx}'
            ax.add_line(
                x=data_set.x,
                y=data_set.y,
                color=color,
                label=data_set.title,
            )

        return figure


if __name__ == '__main__':
    fibermode_solver = FiberModeSolver(wavelength=1550e-9)

    # fibermode_solver.plot_normalized_coupling(
    #     resolution=500,
    #     mode_couples=[('LP01', 'LP02'), ('LP01', 'LP03'), ('LP02', 'LP03')],
    #     itr_list=numpy.linspace(1.0, 0.1, 100)
    # ).show()

    fibermode_solver.plot_beta_vs_itr(
        mode_numbers=['LP01', 'LP02', 'LP03'],
        itr_list=numpy.linspace(1.0, 0.5, 10)
    ).show()

    # fibermode_solver.plot_field_distribution(
    #     itr=0.5,
    #     resolution=300,
    #     mode_numbers=['LP01', 'LP02'],
    #     normalization='max'
    # ).show()

# -
