from typing import Callable
from MPSPlots.styles import mps
import matplotlib.pyplot as plt
from SuPyMode.utils import interpret_mode_of_interest


def simple_plot_helper(function: Callable) -> Callable:
    """
    A decorator that helps in plotting by wrapping a plotting function with additional functionality
    such as handling axes creation, setting the figure style, managing legends, and saving figures.

    Parameters
    ----------
    function : Callable
        The plotting function that is decorated. It should accept `self`, `ax`, and `mode_of_interest`
        as parameters.

    Returns
    -------
    Callable
        A wrapper function that adds the specified plotting functionalities.

    Notes
    -----
    This decorator expects the decorated function to have the following signature:
    `function(self, ax=None, mode_of_interest='all', **kwargs)`.
    """
    def wrapper(self, ax: plt.Axes = None, show: bool = True, save_filename: str = None, **kwargs) -> plt.Figure:
        """
        A wrapped version of the plotting function that provides additional functionality for creating
        and managing plots.

        Parameters
        ----------
        self : object
            The instance of the class calling this method.
        ax : plt.Axes, optional
            A matplotlib Axes object to draw the plot on. If None, a new figure and axes are created.
            Default is None.
        show : bool, optional
            Whether to display the plot. If False, the plot will not be shown but can still be saved
            or returned. Default is True.
        mode_of_interest : str, optional
            Specifies the mode of interest for the plot. If 'all', all available modes will be plotted.
            This parameter is interpreted using the `interpret_mode_of_interest` function. Default is 'all'.
        save_filename : str, optional
            A file path to save the figure. If None, the figure will not be saved. Default is None.
        **kwargs : dict
            Additional keyword arguments passed to the decorated function.

        Returns
        -------
        plt.Figure
            The matplotlib Figure object created or used for the plot.

        Notes
        -----
        - If no `ax` is provided, a new figure and axes are created using the style context `mps`.
        - The legend is only added if there are labels to display.
        - If `save_filename` is specified, the figure is saved to the given path.
        - The plot is shown if `show` is set to True.
        """
        if ax is None:
            with plt.style.context(mps):
                figure, ax = plt.subplots(1, 1)

        else:
            figure = ax.get_figure()

        function(self, ax=ax, **kwargs)

        _, labels = ax.get_legend_handles_labels()

        # Only add a legend if there are labels
        if labels:
            ax.legend()

        if save_filename:
            figure.savefig(save_filename)

        if show:
            plt.show()

        return figure

    return wrapper


def singular_plot_helper(function: Callable) -> Callable:
    """
    A decorator that helps in plotting by wrapping a plotting function with additional functionality
    such as handling axes creation, setting the figure style, managing legends, and saving figures.

    Parameters
    ----------
    function : Callable
        The plotting function that is decorated. It should accept `self`, `ax`, and `mode_of_interest`
        as parameters.

    Returns
    -------
    Callable
        A wrapper function that adds the specified plotting functionalities.

    Notes
    -----
    This decorator expects the decorated function to have the following signature:
    `function(self, ax=None, mode_of_interest='all', **kwargs)`.
    """
    def wrapper(self, ax: plt.Axes = None, show: bool = True, mode_of_interest: str = 'all', save_filename: str = None, **kwargs) -> plt.Figure:
        """
        A wrapped version of the plotting function that provides additional functionality for creating
        and managing plots.

        Parameters
        ----------
        self : object
            The instance of the class calling this method.
        ax : plt.Axes, optional
            A matplotlib Axes object to draw the plot on. If None, a new figure and axes are created.
            Default is None.
        show : bool, optional
            Whether to display the plot. If False, the plot will not be shown but can still be saved
            or returned. Default is True.
        mode_of_interest : str, optional
            Specifies the mode of interest for the plot. If 'all', all available modes will be plotted.
            This parameter is interpreted using the `interpret_mode_of_interest` function. Default is 'all'.
        save_filename : str, optional
            A file path to save the figure. If None, the figure will not be saved. Default is None.
        **kwargs : dict
            Additional keyword arguments passed to the decorated function.

        Returns
        -------
        plt.Figure
            The matplotlib Figure object created or used for the plot.

        Notes
        -----
        - If no `ax` is provided, a new figure and axes are created using the style context `mps`.
        - The legend is only added if there are labels to display.
        - If `save_filename` is specified, the figure is saved to the given path.
        - The plot is shown if `show` is set to True.
        """
        if ax is None:
            with plt.style.context(mps):
                figure, ax = plt.subplots(1, 1)

        else:
            figure = ax.get_figure()

        mode_of_interest = interpret_mode_of_interest(superset=self, mode_of_interest=mode_of_interest)

        function(self, ax=ax, mode_of_interest=mode_of_interest, **kwargs)

        _, labels = ax.get_legend_handles_labels()

        # Only add a legend if there are labels
        if labels:
            ax.legend()

        if save_filename:
            figure.savefig(save_filename)

        if show:
            plt.show()

        return figure

    return wrapper


def combination_plot_helper(function: Callable) -> Callable:
    """
    A decorator that enhances a plotting function by adding functionality for handling axes creation,
    setting the figure style, interpreting combinations and modes, and managing figure display and saving.

    Parameters
    ----------
    function : Callable
        The plotting function that is being decorated. It should accept `self`, `ax`, `mode_of_interest`,
        and `combination` as parameters.

    Returns
    -------
    Callable
        A wrapper function that provides additional functionalities for the plotting function.

    Notes
    -----
    The decorated function should have the following signature:
    `function(self, ax=None, mode_of_interest='all', combination='pairs', **kwargs)`.
    """
    def wrapper(self, ax: plt.Axes = None, show: bool = True, mode_of_interest: str = 'all', combination: str = 'pairs', save_filename: str = None, **kwargs) -> plt.Figure:
        """
        A wrapped version of the plotting function that provides additional functionality for creating
        and managing plots with specific combinations of modes.

        Parameters
        ----------
        self : object
            The instance of the class calling this method.
        ax : plt.Axes, optional
            A matplotlib Axes object to draw the plot on. If None, a new figure and axes are created.
            Default is None.
        show : bool, optional
            Whether to display the plot. If False, the plot will not be shown but can still be saved
            or returned. Default is True.
        mode_of_interest : str, optional
            Specifies the mode of interest for the plot. If 'all', all available modes will be plotted.
            This parameter is interpreted using the `interpret_mode_of_interest` function. Default is 'all'.
        combination : str, optional
            Specifies the type of combination to plot. The value is interpreted using
            `self.interpret_combination` and can be 'pairs', 'triplets', etc. Default is 'pairs'.
        save_filename : str, optional
            A file path to save the figure. If None, the figure will not be saved. Default is None.
        **kwargs : dict
            Additional keyword arguments passed to the decorated function.

        Returns
        -------
        plt.Figure
            The matplotlib Figure object created or used for the plot.

        Notes
        -----
        - If no `ax` is provided, a new figure and axes are created using the style context `mps`.
        - The `mode_of_interest` and `combination` are interpreted based on the instance's methods.
        - The legend is only added if there are labels to display.
        - If `save_filename` is specified, the figure is saved to the given path.
        - The plot is shown if `show` is set to True.
        """
        if ax is None:
            with plt.style.context(mps):
                figure, ax = plt.subplots(1, 1)

        mode_of_interest = interpret_mode_of_interest(superset=self, mode_of_interest=mode_of_interest)

        combination = self.interpret_combination(mode_of_interest=mode_of_interest, combination=combination)

        function(self, ax=ax, mode_of_interest=mode_of_interest, combination=combination, **kwargs)

        _, labels = ax.get_legend_handles_labels()

        # Only add a legend if there are labels
        if labels:
            ax.legend()

        if save_filename:
            figure.savefig(save_filename)

        if show:
            plt.show()

        return figure

    return wrapper


def parse_mode_of_interest(plot_function: Callable) -> Callable:
    """
    A decorator that parses and interprets the `mode_of_interest` parameter for a given plotting function.

    Parameters
    ----------
    plot_function : Callable
        The plotting function to be decorated. It should accept `self`, `mode_of_interest`, and other
        optional arguments.

    Returns
    -------
    Callable
        A wrapper function that interprets the `mode_of_interest` parameter and passes it to the
        decorated plotting function.

    Notes
    -----
    The decorated function should have the following signature:
    `plot_function(self, *args, mode_of_interest='all', **kwargs)`.
    """

    def wrapper(self, *args, mode_of_interest='all', **kwargs):
        """
        A wrapped version of the plotting function that interprets the `mode_of_interest` parameter.

        Parameters
        ----------
        self : object
            The instance of the class calling this method.
        *args : tuple
            Positional arguments passed to the decorated function.
        mode_of_interest : str, optional
            Specifies the mode of interest for the plot. This parameter is interpreted using the
            `interpret_mode_of_interest` function. Default is 'all'.
        **kwargs : dict
            Additional keyword arguments passed to the decorated function.

        Returns
        -------
        Any
            The return value of the decorated plotting function.
        """
        mode_of_interest = interpret_mode_of_interest(
            superset=self,
            mode_of_interest=mode_of_interest
        )

        return plot_function(self, *args, mode_of_interest=mode_of_interest, **kwargs)

    return wrapper


def parse_combination(plot_function: Callable) -> Callable:
    """
    A decorator that parses and interprets the `combination` parameter for a given plotting function,
    in addition to the `mode_of_interest` parameter.

    Parameters
    ----------
    plot_function : Callable
        The plotting function to be decorated. It should accept `self`, `mode_of_interest`, `combination`,
        and other optional arguments.

    Returns
    -------
    Callable
        A wrapper function that interprets the `mode_of_interest` and `combination` parameters and
        passes them to the decorated plotting function.

    Notes
    -----
    The decorated function should have the following signature:
    `plot_function(self, *args, mode_of_interest='all', combination='pairs', **kwargs)`.
    """
    def wrapper(self, *args, mode_of_interest='all', combination: str = 'pairs', **kwargs):
        """
        A wrapped version of the plotting function that interprets the `mode_of_interest` and `combination`
        parameters.

        Parameters
        ----------
        self : object
            The instance of the class calling this method.
        *args : tuple
            Positional arguments passed to the decorated function.
        mode_of_interest : str, optional
            Specifies the mode of interest for the plot. This parameter is interpreted using the
            `interpret_mode_of_interest` function. Default is 'all'.
        combination : str, optional
            Specifies the type of combination to plot. The value is interpreted using
            `self.interpret_combination` and can be 'pairs', 'triplets', etc. Default is 'pairs'.
        **kwargs : dict
            Additional keyword arguments passed to the decorated function.

        Returns
        -------
        Any
            The return value of the decorated plotting function.
        """
        combination = self.interpret_combination(
            mode_of_interest=mode_of_interest,
            combination=combination
        )

        return plot_function(self, *args, mode_of_interest=mode_of_interest, combination=combination, **kwargs)

    return wrapper
