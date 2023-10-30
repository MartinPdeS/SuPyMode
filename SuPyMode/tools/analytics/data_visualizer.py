#!/usr/bin/env python
# -*- coding: utf-8 -*-

# Built-in imports
import numpy

# Third-party imports
from MPSPlots.render2D import SceneList

# Local imports
from SuPyMode.tools.analytics.superset import Superset


class DataVisualizer():
    def __init__(self, wavelength: float):
        self.wavelength = wavelength
        self.superset = Superset(wavelength=wavelength)

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
            data_set = self.superset.get_mode_field_vs_r(
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

            coupling_data_set = self.superset.get_normalized_coupling_vs_itr(
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
            data_set = self.superset.get_beta_vs_itr(
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
            y_label=r'Propagation constant [$\beta$]',
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
    fibermode_solver = DataVisualizer(wavelength=1550e-9)

    # fibermode_solver.plot_normalized_coupling(
    #     resolution=500,
    #     mode_couples=[('LP01', 'LP02'), ('LP01', 'LP03'), ('LP02', 'LP03')],
    #     itr_list=numpy.linspace(1.0, 0.1, 100)
    # ).show()

    fibermode_solver.plot_beta_vs_itr(
        mode_numbers=['LP01', 'LP02', 'LP03'],
        itr_list=numpy.linspace(1.0, 0.1, 10)
    ).show()

    # fibermode_solver.plot_field_distribution(
    #     itr=0.5,
    #     resolution=300,
    #     mode_numbers=['LP01', 'LP02'],
    #     normalization='max'
    # ).show()

# -
