from typing import Callable
from MPSPlots.styles import mps
import matplotlib.pyplot as plt
from SuPyMode.utils import interpret_mode_of_interest


def singular_plot_helper(function: Callable) -> Callable:
    def wrapper(self, ax: plt.Axes = None, show: bool = True, mode_of_interest: str = 'all', **kwargs) -> plt.Figure:
        if ax is None:
            with plt.style.context(mps):
                figure, ax = plt.subplots(1, 1)

        mode_of_interest = interpret_mode_of_interest(superset=self, mode_of_interest=mode_of_interest)

        function(self, ax=ax, mode_of_interest=mode_of_interest, **kwargs)

        _, labels = ax.get_legend_handles_labels()

        # Only add a legend if there are labels
        if labels:
            ax.legend()

        if show:
            plt.show()

        return figure

    return wrapper


def combination_plot_helper(function: Callable) -> Callable:
    def wrapper(self, ax: plt.Axes = None, show: bool = True, mode_of_interest: str = 'all', combination: str = 'pairs', **kwargs) -> plt.Figure:
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

        if show:
            plt.show()

        return figure

    return wrapper


def parse_mode_of_interest(plot_function: Callable) -> Callable:
    def wrapper(self, *args, mode_of_interest='all', **kwargs):
        mode_of_interest = interpret_mode_of_interest(
            superset=self,
            mode_of_interest=mode_of_interest
        )

        return plot_function(self, *args, mode_of_interest=mode_of_interest, **kwargs)

    return wrapper


def parse_combination(plot_function: Callable) -> Callable:
    def wrapper(self, *args, mode_of_interest='all', combination: str = 'pairs', **kwargs):
        combination = self.interpret_combination(
            mode_of_interest=mode_of_interest,
            combination=combination
        )

        return plot_function(self, *args, mode_of_interest=mode_of_interest, combination=combination, **kwargs)

    return wrapper
