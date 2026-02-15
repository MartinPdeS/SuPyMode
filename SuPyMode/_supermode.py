# Built-in imports
import numpy
from dataclasses import dataclass
from scipy.interpolate import RectBivariateSpline

# Local imports
from SuPyMode import representation
from SuPyMode.binary.interface_model_parameters import ModelParameters  # type: ignore
from SuPyMode.utils import interpret_slice_number_and_itr
from SuPyMode.binary.interface_boundaries import Boundaries, BoundaryValue


@dataclass(kw_only=True)
class SuperMode:
    """
    Represents supermodes within fiber optic structures. This class serves as a Python
    counterpart to a C++ SuperMode class, facilitating integration and computation via
    the SuPySolver. Instances of this class belong to a SuperSet, and each supermode
    is uniquely identified within its symmetry set by a mode number.

    Parameters
    ----------
    parent_set : object
        The SuperSet instance associated with this supermode.
    binding : object
        The corresponding C++ bound supermode object.
    solver_number : int
        Identifier linking this supermode to a specific Python solver.
    mode_number : int
        Unique identifier for this mode within a symmetry set.
    boundaries : dict
        Specifications of the boundary conditions for the supermode.

    """

    parent_set: object
    binding: object

    def __post_init__(self):
        self.ID = self.binding.ID
        self.label = self.binding.label
        self.boundaries = self.binding.boundaries
        self.mode_number = self.binding.mode_number
        self.solver_number = self.binding.solver_number
        self.model_parameters = self.binding.model_parameters

        self.is_computation_compatible = self.binding.is_computation_compatible
        self.stylized_label = self.binding.stylized_label
        self.__repr__ = self.binding.__repr__
        self.field = self.binding.field
        self.beating_length = self.binding.beating_length
        self.index = self.binding.index
        self.beta = self.binding.beta
        self.eigenvalue = self.binding.eigenvalue
        self.normalized_coupling = self.binding.normalized_coupling

    @property
    def amplitudes(self) -> numpy.ndarray:
        """
        Computes the amplitude array for this supermode, setting its own mode number
        to 1 and all others to 0.

        Returns
        -------
        numpy.ndarray
            Array of complex numbers representing amplitudes.
        """
        n_mode = len(self.parent_set.supermodes)
        amplitudes = numpy.zeros(n_mode, dtype=complex)
        amplitudes[self.mode_number] = 1
        return amplitudes

    def get_field_interpolation(
        self, itr: float = None, slice_number: int = None
    ) -> RectBivariateSpline:
        """
        Computes the field interpolation for a given iteration or slice number.
        Requires exactly one of the parameters to be specified.

        Parameters
        ----------
        itr : float, optional
            The iteration number for which to compute the interpolation.
        slice_number : int, optional
            The slice number for which to compute the interpolation.

        Returns
        -------
        RectBivariateSpline
            Interpolated field values over a grid.

        Raises
        ------
        ValueError
            If neither or both parameters are specified.
        """
        if (itr is None) == (slice_number is None):
            raise ValueError("Exactly one of itr or slice_number must be provided.")

        if slice_number is None:
            slice_number = self.parent_set.itr_to_slice(itr_list=itr)

        if itr is None:
            slice_number, itr = interpret_slice_number_and_itr(
                itr_baseline=self.itr_list, slice_list=slice_number
            )

        field = self.field.get_field(slice_number=slice_number, add_symmetries=True).T

        x_axis, y_axis = self.get_axis(slice_number=slice_number, add_symmetries=True)

        field_interpolation = RectBivariateSpline(
            x=x_axis * itr,
            y=y_axis * itr,
            z=field,
        )

        return field_interpolation


# -
