# #!/usr/bin/env python
# # -*- coding: utf-8 -*-

# Built-in imports
import numpy
from dataclasses import dataclass, field as field_arg
from scipy.interpolate import RectBivariateSpline

# Local imports
from SuPyMode import representation
from SuPyMode.binary.ModelParameters import ModelParameters
from SuPyMode.utils import interpret_slice_number_and_itr, get_symmetrized_vector


@dataclass(kw_only=True)
class SuperMode():
    """
    Represents supermodes within fiber optic structures. This class serves as a Python
    counterpart to a C++ SuperMode class, facilitating integration and computation via
    the SuPySolver. Instances of this class belong to a SuperSet, and each supermode
    is uniquely identified within its symmetry set by a mode number.

    Attributes:
        parent_set (None): The SuperSet instance associated with this supermode.
        binded_supermode (None): The corresponding C++ bound supermode object.
        solver_number (int): Identifier linking this supermode to a specific Python solver.
        mode_number (int): Unique identifier for this mode within a symmetry set.
        boundaries (dict): Specifications of the boundary conditions for the supermode.
        label (str, optional): An arbitrary descriptive label for the supermode.
    """
    parent_set: object
    binded_supermode: object
    solver_number: int
    mode_number: int
    boundaries: dict
    label: str = None
    ID: list = field_arg(init=False)
    # Other representations
    field: representation.Field = field_arg(init=False)
    index: representation.Index = field_arg(init=False)
    beta: representation.Beta = field_arg(init=False)
    normalized_coupling: representation.NormalizedCoupling = field_arg(init=False)
    adiabatic: representation.Adiabatic = field_arg(init=False)
    eigen_value: representation.EigenValue = field_arg(init=False)
    beating_length: representation.BeatingLength = field_arg(init=False)

    def __post_init__(self):
        self.ID = [self.solver_number, self.binding_number]
        self.field = representation.Field(parent_supermode=self)
        self.index = representation.Index(parent_supermode=self)
        self.beta = representation.Beta(parent_supermode=self)
        self.normalized_coupling = representation.NormalizedCoupling(parent_supermode=self)
        self.adiabatic = representation.Adiabatic(parent_supermode=self)
        self.eigen_value = representation.EigenValue(parent_supermode=self)
        self.beating_length = representation.BeatingLength(parent_supermode=self)

    def __hash__(self):
        """
        Returns a hash based on the binded supermode object, allowing this class
        instance to be used in hash-based collections like sets and dictionaries.

        Returns:
            int: The hash value of the binded supermode object.
        """
        return hash(self.binded_supermode)

    @property
    def binding_number(self) -> int:
        """Retrieves the binding number specific to the linked C++ solver."""
        return self.binded_supermode.binding_number

    @property
    def geometry(self) -> object:
        """
        Provides access to the geometric configuration associated with the supermode.

        Returns:
            object: The geometry of the parent SuperSet.
        """
        return self.parent_set.geometry

    @property
    def coordinate_system(self) -> object:
        """
        Accesses the coordinate system used by the supermode.

        Returns:
            object: The coordinate system of the parent SuperSet.
        """
        return self.parent_set.coordinate_system

    @property
    def itr_list(self) -> numpy.ndarray:
        """Provides a list of iteration parameters from the model."""
        return self.binded_supermode.model_parameters.itr_list

    @property
    def model_parameters(self) -> ModelParameters:
        """
        Retrieves parameters defining the model's computational aspects.

        Returns:
            ModelParameters: Computational parameters from the binded supermode.
        """
        return self.binded_supermode.model_parameters

    @property
    def mesh_gradient(self) -> numpy.ndarray:
        """Accesses the gradient mesh associated with the supermode."""
        return self.binded_supermode.mesh_gradient

    @property
    def amplitudes(self) -> numpy.ndarray:
        """
        Computes the amplitude array for this supermode, setting its own mode number
        to 1 and all others to 0.

        Returns:
            numpy.ndarray: Array of complex numbers representing amplitudes.
        """
        n_mode = len(self.parent_set.supermodes)
        amplitudes = numpy.zeros(n_mode, dtype=complex)
        amplitudes[self.mode_number] = 1
        return amplitudes

    @property
    def stylized_label(self) -> str:
        """
        Provides a stylized label for the supermode. If no custom label is provided,
        it defaults to a generic label with the mode ID.

        Returns:
            str: The stylized or default label.
        """
        if self.label is None:
            return f"Mode: {self.ID}"
        else:
            return f"${self.label}$"

    def is_computation_compatible(self, other: 'SuperMode') -> bool:
        """
        Determines if another supermode is compatible for computation, based on unique
        identifiers and boundary conditions.

        Parameters:
            other (SuperMode): The other supermode to compare.

        Returns:
            bool: True if the supermodes are compatible for computation, False otherwise.
        """
        return self.binded_supermode.is_computation_compatible(other.binded_supermode)

    def is_symmetry_compatible(self, other: 'SuperMode') -> bool:
        """
        Determines whether the specified other supermode has same symmetry.

        :param      other:  The other supermode
        :type       other:  SuperMode

        :returns:   True if the specified other is symmetry compatible, False otherwise.
        :rtype:     bool
        """
        return self.boundaries == other.boundaries

    def get_field_interpolation(self, itr: float = None, slice_number: int = None) -> RectBivariateSpline:
        """
        Computes the field interpolation for a given iteration or slice number. Requires
        exactly one of the parameters to be specified.

        Parameters:
            itr (float, optional): The iteration number for which to compute the interpolation.
            slice_number (int, optional): The slice number for which to compute the interpolation.

        Returns:
            RectBivariateSpline: Interpolated field values over a grid.

        Raises:
            ValueError: If neither or both parameters are specified.
        """
        if (itr is None) == (slice_number is None):
            raise ValueError("Exactly one of itr or slice_number must be provided.")

        if slice_number is None:
            slice_number = self.parent_set.itr_to_slice(itr_list=itr)

        if itr is None:
            slice_number, itr = interpret_slice_number_and_itr(
                itr_baseline=self.itr_list,
                slice_list=slice_number
            )

        field = self.field.get_field(slice_number=slice_number, add_symmetries=True)

        x_axis, y_axis = self.get_axis(slice_number=slice_number, add_symmetries=True)

        field_interpolation = RectBivariateSpline(
            x=x_axis * itr,
            y=y_axis * itr,
            z=field,
        )

        return field_interpolation

    def _get_axis_vector(self, add_symmetries: bool = True) -> tuple:
        full_x_axis = self.coordinate_system.x_vector
        full_y_axis = self.coordinate_system.y_vector

        if not add_symmetries:
            return full_x_axis, full_y_axis

        if self.boundaries.right in ['symmetric', 'anti-symmetric']:
            full_x_axis = get_symmetrized_vector(full_x_axis, symmetry_type='last')
            full_x_axis.sort()

        if self.boundaries.left in ['symmetric', 'anti-symmetric']:
            full_x_axis = get_symmetrized_vector(full_x_axis, symmetry_type='first')
            full_x_axis.sort()

        if self.boundaries.top in ['symmetric', 'anti-symmetric']:
            full_y_axis = get_symmetrized_vector(full_y_axis, symmetry_type='last')
            full_y_axis.sort()

        if self.boundaries.bottom in ['symmetric', 'anti-symmetric']:
            full_y_axis = get_symmetrized_vector(full_y_axis, symmetry_type='first')
            full_y_axis.sort()

        return full_x_axis, full_y_axis

    def get_axis(self, slice_number: int, add_symmetries: bool = True) -> tuple:
        itr = self.model_parameters.itr_list[slice_number]

        x_axis, y_axis = self._get_axis_vector(add_symmetries=add_symmetries)

        return (x_axis * itr, y_axis * itr)

    def __repr__(self) -> str:
        return self.label

    def plot(self, plot_type: str, **kwargs):
        """
        Plots various properties of the supermode based on specified type.

        Parameters:
            plot_type (str): The type of plot to generate (e.g., 'field', 'beta').
            *args: Additional positional arguments for the plot function.
            **kwargs: Additional keyword arguments for the plot function.

        Returns:
            The result of the plotting function, typically a plot object.

        Raises:
            ValueError: If an invalid plot type is specified.
        """
        match plot_type.lower():
            case 'field':
                return self.field.plot(**kwargs)
            case 'beta':
                return self.beta.plot(**kwargs)
            case 'index':
                return self.index.plot(**kwargs)
            case 'eigen-value':
                return self.eigen_value.plot(**kwargs)
            case 'beating-length':
                return self.beating_length.plot(**kwargs)
            case 'adiabatic':
                return self.adiabatic.plot(**kwargs)
            case 'normalized-coupling':
                return self.normalized_coupling.plot(**kwargs)
            case _:
                raise ValueError(f'Invalid plot type: {plot_type}. Options are: index, beta, eigen-value, field, beating-length, adiabatic, normalized-coupling')


# -
