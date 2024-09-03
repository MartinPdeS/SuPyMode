#!/usr/bin/env python
# -*- coding: utf-8 -*-

# Built-in imports
import numpy
from typing import NoReturn
from functools import wraps

# Local imports
from SuPyMode.supermode import SuperMode
from SuPyMode.utils import get_intersection, interpret_mode_of_interest, interpret_slice_number_and_itr, parse_filename
from SuPyMode.profiles import AlphaProfile
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
from SuPyMode.utils import parse_mode_of_interest, parse_combination
from MPSPlots.styles import gg_plot as plot_style

class SuperSetPlots(object):
    # EFFECTIVE INDEX -------------------------------------------------------
    @parse_mode_of_interest
    def _logic_index(self,
            ax: plt.Axes,
            mode_of_interest: list[SuperMode],
            show_crossings: bool = False) -> NoReturn:
        """
        Plot effective index for each mode as a function of itr.

        Args:
            show_crossings (bool): Whether to show crossings in the plot.
            mode_of_interest (str | list[SuperMode]): The mode of interest.

        Returns:
            NoReturn
        """
        for mode in mode_of_interest:
            mode.index.render_on_ax(ax=ax)
            mode.index._dress_ax(ax=ax)

        if show_crossings:
            self.add_crossings_to_ax(ax=ax, mode_of_interest=mode_of_interest, data_type='index')

        plt.legend()

    @wraps(_logic_index)
    def get_figure_index(self, *args, **kwargs):
        figure, ax = plt.subplots(1, 1)
        self._logic_index(ax=ax, *args, **kwargs)
        return figure

    @wraps(_logic_index)
    def plot_index(self, *args, **kwargs):
         with plt.style.context(plot_style):
            figure, ax = plt.subplots(1, 1)
            self._logic_index(ax=ax, *args, **kwargs)
            plt.show()

    # PROPAGATION CONSTANT -------------------------------------------------------
    @parse_mode_of_interest
    def _logic_beta(self,
            ax: plt.Axes,
            mode_of_interest: list[SuperMode],
            show_crossings: bool = False) -> NoReturn:
        """
        Plot propagation constant for each mode as a function of itr.

        Args:
            show_crossings (bool): Whether to show crossings in the plot.
            mode_of_interest (str | list[SuperMode]): The mode of interest.

        Returns:
            NoReturn
        """
        for mode in mode_of_interest:
            mode.beta.render_on_ax(ax=ax)
            mode.beta._dress_ax(ax=ax)

        if show_crossings:
            self.add_crossings_to_ax(ax=ax, mode_of_interest=mode_of_interest, data_type='beta')

        plt.legend()

    @wraps(_logic_beta)
    def get_figure_beta(self, *args, **kwargs):
        figure, ax = plt.subplots(1, 1)
        self._logic_beta(ax=ax, *args, **kwargs)
        return figure

    @wraps(_logic_beta)
    def plot_beta(self, *args, **kwargs):
         with plt.style.context(plot_style):
            figure, ax = plt.subplots(1, 1)
            self._logic_beta(ax=ax, *args, **kwargs)
            plt.show()

    # EIGEN-VALUE -------------------------------------------------------
    @parse_mode_of_interest
    def _logic_eigen_value(self,
            ax: plt.Axes,
            mode_of_interest: list[SuperMode],
            show_crossings: bool = False) -> NoReturn:
        """
        Plot propagation constant for each mode as a function of itr.

        Args:
            mode_of_interest (str | list[SuperMode]): The mode of interest.
            show_crossings (bool): Whether to show crossings in the plot.

        Returns:
            None
        """
        for mode in mode_of_interest:
            mode.eigen_value.render_on_ax(ax=ax)
            mode.eigen_value._dress_ax(ax=ax)

        if show_crossings:
            self.add_crossings_to_ax(ax=ax, mode_of_interest=mode_of_interest, data_type='eigen_value')

        plt.legend()

    @wraps(_logic_eigen_value)
    def get_figure_eigen_value(self, *args, **kwargs):
        figure, ax = plt.subplots(1, 1)
        self._logic_eigen_value(ax=ax, *args, **kwargs)
        return figure

    @wraps(_logic_eigen_value)
    def plot_eigen_value(self, *args, **kwargs):
         with plt.style.context(plot_style):
            figure, ax = plt.subplots(1, 1)
            self._logic_eigen_value(ax=ax, *args, **kwargs)
            plt.show()

    # BEATING LENGTH-------------------------------------------------------
    @parse_mode_of_interest
    @parse_combination
    def _logic_beating_length(self,
            ax: plt.Axes,
            mode_of_interest: list[SuperMode],
            combination: list,
            add_profile: list[AlphaProfile] = []) -> NoReturn:
        """
        Render a figure representing beating_length for each mode as a function of itr.

        Args:
            mode_of_interest (list[SuperMode]): The mode of interest.
            combination (list): The mode combinations.
            add_profile (list[AlphaProfile]): List of profiles to add to the plot.

        Returns:
            None
        """
        for mode_0, mode_1 in combination:
            mode_0.beating_length.render_on_ax(ax=ax, other_supermode=mode_1)
            mode_0.beating_length._dress_ax(ax=ax)

        plt.legend()

    @wraps(_logic_beating_length)
    def get_figure_normalized_coupling(self, *args, **kwargs):
        figure, ax = plt.subplots(1, 1)
        self._logic_beating_length(ax=ax, *args, **kwargs)
        return figure

    @wraps(_logic_beating_length)
    def plot_beating_length(self, *args, **kwargs):
         with plt.style.context(plot_style):
            figure, ax = plt.subplots(1, 1)
            self._logic_beating_length(ax=ax, *args, **kwargs)
            plt.show()


    # NORMALIZED COUPLING-------------------------------------------------------
    @parse_mode_of_interest
    @parse_combination
    def _logic_normalized_coupling(self,
            ax: plt.Axes,
            mode_of_interest: list[SuperMode],
            combination: list,
            add_profile: list[AlphaProfile] = []) -> NoReturn:
        """
        Render a figure representing normalized coupling for each mode as a function of itr.

        Args:
            mode_of_interest (list[SuperMode]): The mode of interest.
            combination (list): The mode combinations.
            add_profile (list[AlphaProfile]): List of profiles to add to the plot.

        Returns:
            None
        """
        for mode_0, mode_1 in combination:
            mode_0.normalized_coupling.render_on_ax(ax=ax, other_supermode=mode_1)
            mode_0.normalized_coupling._dress_ax(ax=ax)

        plt.legend()

    @wraps(_logic_normalized_coupling)
    def get_figure_normalized_coupling(self, *args, **kwargs):
        figure, ax = plt.subplots(1, 1)
        self._logic_normalized_coupling(ax=ax, *args, **kwargs)
        return figure

    @wraps(_logic_normalized_coupling)
    def plot_normalized_coupling(self, *args, **kwargs):
         with plt.style.context(plot_style):
            figure, ax = plt.subplots(1, 1)
            self._logic_normalized_coupling(ax=ax, *args, **kwargs)
            plt.show()


    # ADIABATIC CRITERION-------------------------------------------------------
    @parse_mode_of_interest
    @parse_combination
    def _logic_adiabatic(self,
            ax: plt.Axes,
            mode_of_interest: list[SuperMode],
            combination: list,
            add_profile: list[AlphaProfile] = []) -> NoReturn:
        """
        Render a figure representing adiabatic criterion for each mode as a function of itr.

        Args:
            mode_of_interest (list[SuperMode]): The mode of interest.
            combination (list): The mode combinations.
            add_profile (list[AlphaProfile]): List of profiles to add to the plot.

        Returns:
            None
        """
        for mode_0, mode_1 in combination:
            mode_0.adiabatic.render_on_ax(ax=ax, other_supermode=mode_1)
            mode_0.adiabatic._dress_ax(ax=ax)

        for profile in numpy.atleast_1d(add_profile):
            profile.render_adiabatic_factor_vs_itr_on_ax(ax=ax, line_style='--')

        plt.legend()

    @wraps(_logic_adiabatic)
    def get_figure_adiabatic(self, *args, **kwargs):
        figure, ax = plt.subplots(1, 1)
        self._logic_adiabatic(ax=ax, *args, **kwargs)
        return figure

    @wraps(_logic_adiabatic)
    def plot_adiabatic(self, *args, **kwargs):
        with plt.style.context(plot_style):
            figure, ax = plt.subplots(1, 1)
            self._logic_adiabatic(ax=ax, *args, **kwargs)
            plt.show()

    # FIELD -------------------------------------------------------
    @parse_mode_of_interest
    def _logic_field(
            self,
            mode_of_interest: list = 'all',
            itr_list: list[float] = None,
            slice_list: list[int] = None,
            show_mode_label: bool = True,
            show_itr: bool = True,
            show_slice: bool = True) -> plt.Figure:
        """
        Render the mode field for different ITR values or slice numbers.

        Args:
            mode_of_interest (list): List of modes to be plotted. Default is 'all'.
            itr_list (list): List of ITR values for plotting. Default is None.
            slice_list (list): List of slice numbers for plotting. Default is None.
            show_mode_label (bool): Flag to display mode labels. Default is True.
            show_itr (bool): Flag to display ITR values. Default is True.
            show_slice (bool): Flag to display slice numbers. Default is True.

        Returns:
            plt.Figure: The figure object containing the generated plots.
        """
        # Interpret input lists
        slice_list, itr_list = interpret_slice_number_and_itr(
            itr_baseline=self.model_parameters.itr_list,
            itr_list=itr_list,
            slice_list=slice_list
        )

        mode_of_interest = interpret_mode_of_interest(
            superset=self,
            mode_of_interest=mode_of_interest
        )

        # Determine the grid size for subplots
        grid_size = numpy.array([len(slice_list), len(mode_of_interest)])
        figure, axes = plt.subplots(*grid_size, figsize=3 * numpy.flip(grid_size), squeeze=False)

        # Plot each mode field on the grid
        for m, mode in enumerate(mode_of_interest):
            for n, slice_number in enumerate(slice_list):
                mode.field.render_on_ax(
                    ax=axes[n, m],
                    slice_number=slice_number,
                    show_mode_label=show_mode_label,
                    show_itr=show_itr,
                    show_slice=show_slice
                )

        figure.tight_layout()
        return figure

    @wraps(_logic_field)
    def get_figure_field(self, *args, **kwargs):
        figure = self._logic_field(*args, **kwargs)
        return figure

    @wraps(_logic_field)
    def plot_field(self, *args, **kwargs):
         with plt.style.context(plot_style):
            figure = self._logic_field(*args, **kwargs)
            plt.show()

    def plot(self, plot_type: str, **kwargs) -> NoReturn:
        """
        General plotting function to handle different types of supermode plots.

        Args:
            plot_type (str): The type of plot to generate. Options include 'index', 'beta', 'eigen-value', etc.
            **kwargs: Additional keyword arguments for specific plot configurations.

        Raises:
            ValueError: If an unrecognized plot type is specified.
        """
        match plot_type.lower():
            case 'index':
                return self.plot_index(**kwargs)
            case 'beta':
                return self.plot_beta(**kwargs)
            case 'eigen-value':
                return self.plot_eigen_value(**kwargs)
            case 'normalized-coupling':
                return self.plot_normalized_coupling(**kwargs)
            case 'overlap':
                return self.plot_overlap(**kwargs)
            case 'adiabatic':
                return self.plot_adiabatic(**kwargs)
            case 'field':
                return self.plot_field(**kwargs)
            case 'beating-length':
                return self.plot_beating_length(**kwargs)
            case 'normalized-adiabatic':
                return self.plot_normalized_adiabatic(**kwargs)
            case _:
                raise ValueError(f'Invalid plot type: {plot_type}. Options are: index, beta, eigen-value, adiabatic, normalized-adiabatic, normalized-coupling, field, beating-length')

    @parse_filename
    def generate_pdf_report(
            self,
            filename: str = "auto",
            directory: str = '.',
            itr_list: list[float] | None = None,
            slice_list: list[int] | None = None,
            dpi: int = 200,
            mode_of_interest: list = 'all',
            combination: str = 'specific') -> None:
        """
        Generate a full report of the coupler properties as a .pdf file.

        Args:
            filename (str): Name of the report file to be output.
            directory (str): Directory to save the report.
            itr_list (List[float]): List of ITR values to evaluate the mode field.
            slice_list (List[int]): List of slice values to evaluate the mode field.
            dpi (int): Pixel density for the images included in the report.
            mode_of_interest (List): List of modes to consider in the adiabatic criterion plotting.
            combination (str): Method for selecting mode combinations.

        Returns:
            None
        """

        figure_list = [
            self.geometry.render_plot(),
            self.get_figure_field(itr_list=itr_list, slice_list=slice_list, mode_of_interest=mode_of_interest),
            self.get_figure_index(mode_of_interest=mode_of_interest),
            self.get_figure_beta(mode_of_interest=mode_of_interest),
            self.get_figure_normalized_coupling(mode_of_interest=mode_of_interest, combination=combination),
            self.get_figure_adiabatic(mode_of_interest=mode_of_interest, combination=combination)
        ]

        pp = PdfPages(filename.with_suffix('.pdf'))

        for fig in figure_list:
            fig.savefig(pp, format='pdf')

        pp.close()

        plt.close()

    def add_crossings_to_ax(self, ax: plt.Axes, mode_of_interest: list, data_type: str) -> None:
        combination = self.interpret_combination(
            mode_of_interest=mode_of_interest,
            combination='pairs'
        )

        for mode_0, mode_1 in combination:
            x, y = get_intersection(
                x=self.model_parameters.itr_list,
                y0=getattr(mode_0, data_type).data,
                y1=getattr(mode_1, data_type).data,
                average=True
            )

            if x is not None:
                ax.scatter(x=x, y=y, marker='o', color='black', s=20, label='mode crossing')


# -
