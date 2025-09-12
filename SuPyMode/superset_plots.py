#!/usr/bin/env python
# -*- coding: utf-8 -*-

# Built-in imports
import matplotlib.pyplot as plt
import numpy
from matplotlib.backends.backend_pdf import PdfPages
from MPSPlots import helper

from SuPyMode.profiles import AlphaProfile
from SuPyMode.supermode import SuperMode
from SuPyMode.utils import (
    get_intersection,
    interpret_mode_of_interest,
    interpret_slice_number_and_itr,
    parse_filename,
)


class SuperSetPlots(object):
    @helper.pre_plot(nrows=1, ncols=1)
    def plot_index(
        self,
        axes: plt.Axes,
        mode_of_interest: list[SuperMode] = "all",
        show_crossings: bool = False,
    ) -> plt.Figure:
        """
        Plot the effective index for each mode as a function of inverse taper ratio (ITR).

        Parameters
        ----------
        axes : matplotlib.axes.Axes
            The axes object on which to plot.
        mode_of_interest : list of SuperMode
            List of modes to be plotted.
        show_crossings : bool, optional
            Whether to show crossings in the plot (default is False).

        Examples
        --------
        >>> fig, ax = plt.subplots()
        >>> superset_plots.plot_index(ax=ax, mode_of_interest=[mode1, mode2], show_crossings=True)
        >>> plt.show()

        Returns
        -------
        plt.Figure
            The figure object containing the generated plots. This can be customized further or saved to a file.
        """
        mode_of_interest = interpret_mode_of_interest(
            superset=self, mode_of_interest=mode_of_interest
        )

        for mode in mode_of_interest:
            mode.index.plot(ax=axes, show=False)

        if show_crossings:
            self.add_crossings_to_ax(
                ax=axes, mode_of_interest=mode_of_interest, data_type="index"
            )

        axes.legend()

    @helper.pre_plot(nrows=1, ncols=1)
    def plot_beta(
        self,
        axes: plt.Axes,
        mode_of_interest: list[SuperMode] | str = "all",
        show_crossings: bool = False,
    ) -> plt.Figure:
        """
        Plot the effective propagation constant for each mode as a function of inverse taper ratio (ITR).

        Parameters
        ----------
        axes : matplotlib.axes.Axes
            The axes object on which to plot.
        mode_of_interest : list of SuperMode
            List of modes to be plotted.
        show_crossings : bool, optional
            Whether to show crossings in the plot (default is False).

        Examples
        --------
        >>> fig, ax = plt.subplots()
        >>> superset_plots.plot_index(ax=ax, mode_of_interest=[mode1, mode2], show_crossings=True)
        >>> plt.show()

        Returns
        -------
        plt.Figure
            The figure object containing the generated plots. This can be customized further or saved to a file.
        """
        mode_of_interest = interpret_mode_of_interest(
            superset=self, mode_of_interest=mode_of_interest
        )

        for mode in mode_of_interest:
            mode.beta.plot(ax=axes, show=False)

        if show_crossings:
            self.add_crossings_to_ax(
                ax=axes, mode_of_interest=mode_of_interest, data_type="beta"
            )

        axes.legend()

    @helper.pre_plot(nrows=1, ncols=1)
    def plot_eigen_value(
        self,
        axes: plt.Axes,
        mode_of_interest: list[SuperMode] | str = "all",
        show_crossings: bool = False,
    ) -> plt.Figure:
        """
        Plot the computed eigen-values for each mode as a function of inverse taper ratio (ITR).

        Parameters
        ----------
        axes : matplotlib.axes.Axes
            The axes object on which to plot.
        mode_of_interest : list of SuperMode
            List of modes to be plotted.
        show_crossings : bool, optional
            Whether to show crossings in the plot (default is False).

        Examples
        --------
        >>> fig, ax = plt.subplots()
        >>> superset_plots.plot_index(ax=ax, mode_of_interest=[mode1, mode2], show_crossings=True)
        >>> plt.show()

        Returns
        -------
        plt.Figure
            The figure object containing the generated plots. This can be customized further or saved to a file.
        """
        mode_of_interest = interpret_mode_of_interest(
            superset=self, mode_of_interest=mode_of_interest
        )

        for mode in mode_of_interest:
            mode.eigen_value.plot(ax=axes, show=False)

        if show_crossings:
            self.add_crossings_to_ax(
                ax=axes, mode_of_interest=mode_of_interest, data_type="eigen_value"
            )

        axes.legend()

    @helper.pre_plot(nrows=1, ncols=1)
    def plot_beating_length(
        self,
        axes: plt.Axes,
        mode_of_interest: list[SuperMode] = "all",
        combination: list = "pairs",
    ) -> plt.Figure:
        """
        Render a figure representing the beating length for each mode combination as a function of inverse taper ratio (ITR).

        Parameters
        ----------
        axes : matplotlib.axes.Axes
            The axes object on which to plot.
        mode_of_interest : list of SuperMode
            List of modes to be plotted.
        combination : list
            List of mode combinations.

        Returns
        -------
        plt.Figure
            The figure object containing the generated plots. This can be customized further or saved to a file.
        """
        mode_of_interest = interpret_mode_of_interest(
            superset=self, mode_of_interest=mode_of_interest
        )

        combination = self.interpret_combination(
            mode_of_interest=mode_of_interest, combination=combination
        )

        for mode_0, mode_1 in combination:
            mode_0.beating_length.plot(ax=axes, other_supermode=mode_1)

        axes.legend()

    @helper.pre_plot(nrows=1, ncols=1)
    def plot_normalized_coupling(
        self,
        axes: plt.Axes,
        mode_of_interest: list[SuperMode] | str = "all",
        combination: list | str = "pairs",
    ) -> plt.Figure:
        """
        Render a figure representing the normalized coupling for each mode combination as a function of inverse taper ratio (ITR).

        Parameters
        ----------
        ax : matplotlib.axes.Axes
            The axes object on which to plot.
        mode_of_interest : list of SuperMode
            List of modes to be plotted.
        combination : list
            List of mode combinations.

        Returns
        -------
        plt.Figure
            The figure object containing the generated plots. This can be customized further or saved to a file.
        """
        mode_of_interest = interpret_mode_of_interest(
            superset=self, mode_of_interest=mode_of_interest
        )

        combination = self.interpret_combination(
            mode_of_interest=mode_of_interest, combination=combination
        )

        for mode_0, mode_1 in combination:
            mode_0.normalized_coupling.plot(ax=axes, other_supermode=mode_1, show=False)

    @helper.pre_plot(nrows=1, ncols=1)
    def plot_adiabatic(
        self,
        axes: plt.Axes,
        mode_of_interest: list[SuperMode] | str = "all",
        combination: list | str = "pairs",
        add_profile: list[AlphaProfile] = None,
    ) -> plt.Figure:
        """
        Render a figure representing the adiabatic criterion for each mode combination as a function of inverse taper ratio (ITR).

        Parameters
        ----------
        ax : matplotlib.axes.Axes
            The axes object on which to plot.
        mode_of_interest : list of SuperMode
            List of modes to be plotted.
        combination : list
            List of mode combinations.
        add_profile : list of AlphaProfile, optional
            List of profiles to add to the plot (default is an empty list).

        Returns
        -------
        plt.Figure
            The figure object containing the generated plots. This can be customized further or saved to a file.
        """
        mode_of_interest = interpret_mode_of_interest(
            superset=self, mode_of_interest=mode_of_interest
        )

        combination = self.interpret_combination(
            mode_of_interest=mode_of_interest, combination=combination
        )

        for mode_0, mode_1 in combination:
            mode_0.adiabatic.plot(ax=axes, other_supermode=mode_1, show=False)

        if add_profile is not None:
            for profile in numpy.atleast_1d(add_profile):
                profile.render_adiabatic_factor_vs_itr_on_ax(ax=axes, line_style="--")

    @helper.post_mpl_plot
    def plot_field(
        self,
        mode_of_interest: list = "all",
        itr_list: list[float] = None,
        slice_list: list[int] = None,
        show_mode_label: bool = True,
        show_itr: bool = True,
        show_slice: bool = True,
    ) -> plt.Figure:
        """
        Render the mode field for different ITR values or slice numbers.

        Parameters
        ----------
        mode_of_interest : list, optional
            List of modes to be plotted (default is 'all').
        itr_list : list of float, optional
            List of ITR values for plotting (default is None).
        slice_list : list of int, optional
            List of slice numbers for plotting (default is None).
        show_mode_label : bool, optional
            Flag to display mode labels (default is True).
        show_itr : bool, optional
            Flag to display ITR values (default is True).
        show_slice : bool, optional
            Flag to display slice numbers (default is True).

        Returns
        -------
        plt.Figure
            The figure object containing the generated plots. This can be customized further or saved to a file.
        """
        mode_of_interest = interpret_mode_of_interest(
            superset=self, mode_of_interest=mode_of_interest
        )

        # Interpret input lists
        slice_list, itr_list = interpret_slice_number_and_itr(
            itr_baseline=self.model_parameters.itr_list,
            itr_list=itr_list,
            slice_list=slice_list,
        )

        mode_of_interest = interpret_mode_of_interest(
            superset=self, mode_of_interest=mode_of_interest
        )

        n_mode = len(mode_of_interest)
        n_slice = len(slice_list)
        grid_size = numpy.array([n_slice, n_mode])

        figure, axes = plt.subplots(
            *grid_size,
            figsize=3 * numpy.flip(grid_size),
            squeeze=False,
            sharex=False,
            sharey=False,
        )

        # Plot each mode field on the grid
        for m, mode in enumerate(mode_of_interest):
            for n, slice_number in enumerate(slice_list):
                mode.field.render_on_ax(
                    ax=axes[n, m],
                    slice_number=slice_number,
                    show_mode_label=show_mode_label,
                    show_itr=show_itr,
                    show_slice=show_slice,
                )

        return figure

    def plot(self, plot_type: str, **kwargs) -> None:
        """
        General plotting function to handle different types of supermode plots.

        Parameters
        ----------
        plot_type : str
            The type of plot to generate. Options include 'index', 'beta', 'eigen-value', etc.
        **kwargs : dict
            Additional keyword arguments for specific plot configurations.

        Raises
        ------
        ValueError
            If an unrecognized plot type is specified.

        Examples
        --------
        >>> superset_plots.plot(plot_type='index', ax=ax)
        Generates an effective index plot.

        >>> superset_plots.plot(plot_type='invalid')
        ValueError: Invalid plot type: invalid. Options are: index, beta, eigen-value, adiabatic, normalized-adiabatic, normalized-coupling, field, beating-length.
        """
        match plot_type.lower().replace("_", "-"):
            case "index":
                return self.plot_index(**kwargs)
            case "beta":
                return self.plot_beta(**kwargs)
            case "eigen-value":
                return self.plot_eigen_value(**kwargs)
            case "normalized-coupling":
                return self.plot_normalized_coupling(**kwargs)
            case "overlap":
                return self.plot_overlap(**kwargs)
            case "adiabatic":
                return self.plot_adiabatic(**kwargs)
            case "field":
                return self.plot_field(**kwargs)
            case "beating-length":
                return self.plot_beating_length(**kwargs)
            case "normalized-adiabatic":
                return self.plot_normalized_adiabatic(**kwargs)
            case _:
                raise ValueError(
                    f"Invalid plot type: {plot_type}. Options are: index, beta, eigen-value, adiabatic, normalized-adiabatic, normalized-coupling, field, beating-length"
                )

    @parse_filename
    def generate_pdf_report(
        self,
        filename: str = "auto",
        itr_list: list[float] | None = None,
        slice_list: list[int] | None = None,
        mode_of_interest: list[SuperMode] | str = "all",
        combination: str = "pairs",
    ) -> None:
        """
        Generate a full report of the coupler properties as a PDF file.

        Parameters
        ----------
        filename : str, optional, default="auto"
            Name of the report file to be output. If "auto", a default name will be generated based on the current timestamp.
        directory : str, optional, default='.'
            Directory to save the report.
        itr_list : list of float, optional, default=None
            List of ITR values to evaluate the mode field. If None, all available ITR values will be used.
        slice_list : list of int, optional, default=None
            List of slice values to evaluate the mode field. If None, all available slices will be used.
        mode_of_interest : list, optional, default='all'
            List of modes to consider in the report. If 'all', all available modes will be included.
        combination : str, optional, default='specific'
            Method for selecting mode combinations ('specific' or 'pairs').

        Examples
        --------
        >>> superset_plots.generate_pdf_report(filename="coupler_report", itr_list=[0.1, 0.2, 0.5], mode_of_interest=['LP01', 'LP11'])
        This will generate a PDF report named 'coupler_report.pdf' containing plots for the specified modes at given ITR values.
        """
        kwargs = dict(show=False, mode_of_interest=mode_of_interest)

        figure_list = [
            self.geometry.plot(show=False),
            self.plot_field(itr_list=itr_list, slice_list=slice_list, **kwargs),
            self.plot_index(**kwargs),
            self.plot_beta(**kwargs),
            self.plot_normalized_coupling(**kwargs, combination=combination),
            self.plot_adiabatic(show=False, combination=combination),
        ]

        pp = PdfPages(filename.with_suffix(".pdf"))

        for fig in figure_list:
            fig.savefig(pp, format="pdf")

        pp.close()

        plt.close()

    def add_crossings_to_ax(
        self, ax: plt.Axes, mode_of_interest: list, data_type: str
    ) -> None:
        """
        Add mode crossings to the given axis.

        Parameters
        ----------
        ax : matplotlib.axes.Axes
            The axes object on which to plot.
        mode_of_interest : list of SuperMode
            List of modes of interest.
        data_type : str
            The type of data for which to find crossings (e.g., 'index', 'beta', etc.).

        """
        combination = self.interpret_combination(
            mode_of_interest=mode_of_interest, combination="pairs"
        )

        for mode_0, mode_1 in combination:
            x, y = get_intersection(
                x=self.model_parameters.itr_list,
                y0=getattr(mode_0, data_type).data,
                y1=getattr(mode_1, data_type).data,
                average=True,
            )

            if x is not None:
                ax.scatter(
                    x=x, y=y, marker="o", color="black", s=20, label="mode crossing"
                )
