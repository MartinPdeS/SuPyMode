# #!/usr/bin/env python
# # -*- coding: utf-8 -*-

# Built-in imports
import numpy
from dataclasses import dataclass
from scipy.interpolate import RectBivariateSpline

# Local imports
from SuPyMode import representation
from SuPyMode.binary.ModelParameters import ModelParameters
from SuPyMode.tools.utils import interpret_slice_number_and_itr


class InheritFromSuperSet():
    """
    Property class for inherited attribute from SuperSet.

    """
    @property
    def geometry(self) -> object:
        return self.parent_set.geometry

    @property
    def coordinate_system(self) -> object:
        return self.parent_set.coordinate_system

    @property
    def itr_list(self) -> numpy.ndarray:
        return self.binded_supermode.model_parameters.itr_list


@dataclass
class SuperMode(InheritFromSuperSet):
    """
    This class is a representation of the fiber optic structures SuperModes.
    Those mode belongs to a SuperSet class and are constructed with the SuPySolver.
    It links to c++ SuperMode class.
    """
    parent_set: None
    """SuperSet to which is associated the computed this mode"""
    binded_supermode: None
    """C++ binded sueprmode"""
    solver_number: int
    """Number which bind this mode to a specific python solver"""
    mode_number: int
    """Unique number associated to this mode in a particular symmetry set"""
    boundaries: dict
    """Boundary conditions"""
    label: str = None
    """Arbitrary label given to the mode"""

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
        return hash(self.binded_supermode)

    @property
    def binding_number(self) -> int:
        """ Returns the mode number specific to one CppSolver """
        return self.binded_supermode.binding_number

    @property
    def model_parameters(self) -> ModelParameters:
        return self.binded_supermode.model_parameters

    @property
    def mesh_gradient(self) -> numpy.ndarray:
        return self.binded_supermode.mesh_gradient

    @property
    def amplitudes(self) -> numpy.ndarray:
        n_mode = len(self.parent_set.supermodes)
        amplitudes = numpy.zeros(n_mode, dtype=complex)
        amplitudes[self.mode_number] = 1
        return amplitudes

    @property
    def stylized_label(self) -> str:
        if self.label is None:
            return f"Mode: {self.ID}"
        else:
            return f"${self.label}$"

    def is_computation_compatible(self, other: 'SuperMode') -> bool:
        """
        Determines whether the specified other supermode is compatible
        for computation of the modal coupling and adiabatic criterion.
        It, basically return False only if the mode is the same or if the
        boundaries symmetries differ in some way.

        :param      other:  The other SuperMode to compare with
        :type       other:  SuperMode

        :returns:   True if the specified other is computation compatible, False otherwise.
        :rtype:     bool
        """
        if self.ID != other.ID and self.is_symmetry_compatible(other):
            return True
        else:
            return False

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
        Gets the mode field interpolation at a certain itr or slice number value.

        :param      itr:           The itr
        :type       itr:           float
        :param      slice_number:  The slice number
        :type       slice_number:  int

        :returns:   The mode field interpolation.
        :rtype:     RectBivariateSpline
        """
        if not (itr is None) ^ (slice_number is None):
            raise ValueError("Exactly one of the two keyword argument [itr, slice_number] has to be provided.")

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

    def _get_symmetrize_vector(self, vector: numpy.ndarray, symmetry_type: str = 'last') -> numpy.ndarray:
        """
        Take as input a vector and return a symmetric version of that vector.
        The symmetry can be setted mirror the last or first element.

        :param      vector:          The vector
        :type       vector:          numpy.ndarray
        :param      symmetry_type:   The symmetry type
        :type       symmetry_type:   str

        :returns:   The symmetrized vector
        :rtype:     numpy.ndarray

        :raises     AssertionError:  Verify that input vector is 1-dimensionnal.
        """
        assert vector.ndim == 1, f'Vector should be 1d, instead {vector.ndim} dimensional is provided.'

        size = len(vector)
        dx = abs(vector[0] - vector[1])

        match symmetry_type.lower():
            case 'last':
                start_value = vector[0]
                return numpy.arange(0, 2 * size - 1) * dx + start_value

            case 'first':
                start_value = vector[-1]
                return numpy.arange(0, 2 * size - 1) * -dx + start_value

    def _get_axis_vector(self, add_symmetries: bool = True) -> tuple:
        full_x_axis = self.coordinate_system.x_vector
        full_y_axis = self.coordinate_system.y_vector

        if not add_symmetries:
            return full_x_axis, full_y_axis

        if self.boundaries.right in ['symmetric', 'anti-symmetric']:
            full_x_axis = self._get_symmetrize_vector(full_x_axis, symmetry_type='last')
            full_x_axis.sort()

        if self.boundaries.left in ['symmetric', 'anti-symmetric']:
            full_x_axis = self._get_symmetrize_vector(full_x_axis, symmetry_type='first')
            full_x_axis.sort()

        if self.boundaries.top in ['symmetric', 'anti-symmetric']:
            full_y_axis = self._get_symmetrize_vector(full_y_axis, symmetry_type='last')
            full_y_axis.sort()

        if self.boundaries.bottom in ['symmetric', 'anti-symmetric']:
            full_y_axis = self._get_symmetrize_vector(full_y_axis, symmetry_type='first')
            full_y_axis.sort()

        return full_x_axis, full_y_axis

    def get_axis(self, slice_number: int, add_symmetries: bool = True) -> tuple:
        itr = self.model_parameters.itr_list[slice_number]

        x_axis, y_axis = self._get_axis_vector(add_symmetries=add_symmetries)

        return (x_axis * itr, y_axis * itr)

    def __repr__(self) -> str:
        return self.label

    def plot(self, plot_type: str, *args, **kwargs):
        match plot_type.lower():
            case 'field':
                return self.field.plot(*args, **kwargs)
            case 'beta':
                return self.beta.plot(*args, **kwargs)
            case 'index':
                return self.beta.plot(*args, **kwargs)
            case 'eigen-value':
                return self.eigen_value.plot(*args, **kwargs)
            case 'beating-length':
                return self.beating_length.plot(*args, **kwargs)

# -
