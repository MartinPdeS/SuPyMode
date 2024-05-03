#!/usr/bin/env python
# -*- coding: utf-8 -*-

# Standard imports
import numpy
from dataclasses import dataclass
from MPSPlots.CMAP import BKR
from MPSPlots.render2D import Scene2D, Axis, Line, Mesh, ColorBar


@dataclass
class ModeSolver():
    itr_list: numpy.ndarray
    """ List of value representing the inverse taper ratio """
    shape: tuple
    """ Tuple representing the shape of the mesh from which the eigen matrix is constructred """
    wavelength: float
    """ Value representing the wavelength of the simulation """
    label: str = ''
    """ Label of the mode """
    fitting_max_degree: int = 5
    """ Max value for the coefficient fitting algorithm for eigen value extrapolation """
    save_fittings: bool = True

    def __post_init__(self):
        self.wavenumber = 2 * numpy.pi / self.wavelength
        self.construct_initial_arrays()
        self.fitting_coefficients = None

    def construct_initial_arrays(self) -> None:
        """
        Construct the arrays where the eigen values and vectors will be stored.

        :returns:   No return
        :rtype:     None
        """
        eigen_vector_size = self.shape[0] * self.shape[1]
        self.eigen_vectors = numpy.full([self.itr_list.size, eigen_vector_size], numpy.nan)
        self.eigen_values = numpy.full(self.itr_list.size, numpy.nan)
        self.effective_index = numpy.full(self.itr_list.size, numpy.nan)
        self.fit_values = numpy.full(self.itr_list.size, numpy.nan)

        self.eigen_values = self.eigen_values.astype(complex)
        self.effective_index = self.effective_index
        self.eigen_vectors = self.eigen_vectors.astype(complex)
        self.fit_values = self.fit_values.astype(complex)

    def index_to_eigen_value(self, itr: float, index: float) -> float:
        """
        Return the eigen value associated to the given effective index

        :param      index:  The mode effective index
        :type       index:  float

        :returns:   The associated eigen value
        :rtype:     float
        """
        return (self.wavenumber * index * itr)**2

    def eigen_value_to_index(self, eigen_value: float, itr: float) -> float:
        """
        Return the effective index associated toe the given eigen value.

        :param      eigen_value:  The eigen value
        :type       eigen_value:  float
        :param      itr:          The itr
        :type       itr:          float

        :returns:   The associated effective index
        :rtype:     float
        """
        return numpy.sqrt(eigen_value) / self.wavenumber / itr

    def get_next_slice_idx(self) -> int:
        """
        Return the index of the next slice that will be computed.

        :returns:   The next slice index.
        :rtype:     int
        """
        non_nan_idx = numpy.argwhere(~numpy.isnan(self.eigen_values))
        return non_nan_idx.size

    def get_last_slice_idx(self) -> int:
        """
        Return the index of the last slice that was computed.

        :returns:   The last slice index.
        :rtype:     int
        """
        next_slice_idx = self.get_next_slice_idx()
        return next_slice_idx - 1

    def append_eigen_value(self, itr: float, eigen_value: float) -> None:
        """
        Append eigen_value to the mode. Stored in self.eigen_values.

        :param      itr:  The corresponding itr
        :type       itr:  float

        :param      eigen_value:  The eigen value
        :type       eigen_value:  float

        :returns:   No return
        :rtype:     None
        """
        idx = self.get_next_slice_idx()

        self.eigen_values[idx] = eigen_value
        self.effective_index[idx] = self.eigen_value_to_index(eigen_value=eigen_value, itr=itr).real

    def append_eigen_vector(self, eigen_vector: numpy.ndarray) -> None:
        """
        Append eigen_vector to the mode. Stored in self.eigen_values.

        :param      eigen_vector:  The eigen vector
        :type       eigen_vector:  numpy.ndarray

        :returns:   No return
        :rtype:     None
        """
        idx = self.get_next_slice_idx()
        self.eigen_vectors[idx] = eigen_vector

    def append_next_slice(self, eigen_value: float, eigen_vector: numpy.ndarray, itr: float) -> None:
        """
        Append eigen_value to the mode. Stored in self.eigen_values.

        :param      itr:  The corresponding itr
        :type       itr:  float

        :param      eigen_value:  The eigen value
        :type       eigen_value:  float
        :param      eigen_vector:  The eigen vector
        :type       eigen_vector:  numpy.ndarray

        :returns:   No return
        :rtype:     None
        """
        last_slice = self.get_last_slice_idx()
        if last_slice > 0:
            projection_with_previous = abs(eigen_vector.dot(self.last_eigen_vector))
            if projection_with_previous < 0.5:
                print(f"Error bad mode matching with {self.label} -> [{last_slice=}] [{projection_with_previous = :.3f}]")

        self.append_eigen_vector(eigen_vector=eigen_vector)
        self.append_eigen_value(eigen_value=eigen_value, itr=itr)

    @property
    def last_eigen_vector(self) -> numpy.ndarray:
        """
        Return the last eigen vector that was computed.

        :returns:   The eigen vector
        :rtype:     numpy.ndarray
        """
        idx = self.get_last_slice_idx()
        return self.eigen_vectors[idx]

    @property
    def last_eigen_value(self) -> float:
        """
        Return the last eigen value that was computed.

        :returns:   The eigen value
        :rtype:     float
        """
        idx = self.get_last_slice_idx()
        return self.eigen_values[idx]

    def compute_fitting_coefficients(self, n_points: int = 10, max_degree: int = None):
        """
        Compute the fitting coefficient for the n last points:

        :param      max_degree:  The maximum degree of the fitting, if None is the default one.
        :type       max_degree:  int

        :param      n_points:  The n last points to evaluate for the fitting
        :type       n_points:  int
        """
        if max_degree is None:
            max_degree = self.fitting_max_degree

        idx = self.get_next_slice_idx()
        x = self.itr_list[:idx]
        y = self.eigen_values[:idx]

        x = x[-n_points:]
        y = y[-n_points:]

        if max_degree > x.size - 2:
            max_degree = x.size - 2

        if max_degree < 1:
            max_degree = 1

        self.fitting_coefficients = numpy.polyfit(x, y, max_degree)

    def evaluate_to_fitting(self, itr_list: float) -> numpy.ndarray:
        """
        Return the eigen values the was evaluated ovec the polynomial fitting.

        :param      itr_list:  The itr list
        :type       itr_list:  numpy.ndarray

        :returns:   The evaluated eigen values
        :rtype:     numpy.ndarray
        """
        assert self.fitting_coefficients is not None, "Error trying the fit the data without having computed the coefficient prior."

        fit_value = numpy.polyval(
            self.fitting_coefficients,
            itr_list
        )

        return fit_value

    def extrapolate_eigen_value(self, itr: float) -> float:
        """
        Extrapolate the eigen value of the mode at the given itr value.

        :param      itr:  The itr
        :type       itr:  float

        :returns:   The evaluated eigen value
        :rtype:     float
        """
        self.compute_fitting_coefficients()

        extrapolated_value = self.evaluate_to_fitting(itr_list=itr)

        if self.save_fittings:
            slice_number = self.itr_to_slice(itr=itr)
            self.fit_values[slice_number] = extrapolated_value

        return extrapolated_value

    def extrapolate_eigen_vector(self, itr: float) -> numpy.ndarray:
        """
        Extrapolate the eigen vector of the mode at the given itr value.
        This function dont truly extrapolate it just return the eigen vector
        which is the closest in term of itr.

        :param      itr:  The itr
        :type       itr:  float

        :returns:   The evaluated eigen vector
        :rtype:     numpy.ndarray
        """
        slice_number = self.itr_to_slice(itr=itr)

        return self.eigen_vectors[slice_number]

    def itr_to_slice(self, itr: float) -> int:
        """
        Return the slice that correspond the closest to the itr value given.

        :param      itr:  The itr
        :type       itr:  float

        :returns:   The slice number
        :rtype:     int
        """
        closest_idx = (numpy.abs(self.itr_list - itr)).argmin()
        return closest_idx

    @property
    def fields(self) -> numpy.ndarray:
        """
        Return the fields of the mode, the first axis is the itr_slice

        :returns:   The fields array
        :rtype:     numpy.ndarray
        """
        return self.eigen_vectors.reshape([self.itr_list.size, *self.shape]).real

    @property
    def get_field(self, itr: float) -> numpy.ndarray:
        """
        Return the fields of the mode, the first axis is the itr_slice

        :returns:   The fields array
        :rtype:     numpy.ndarray
        """
        return self.eigen_vectors.reshape([self.itr_list.size, *self.shape]).real

    def render_field(self, itr: float, slice_number, ax: Axis) -> None:
        """
        Render the mode field evaluated at a certain itr and on a given ax.

        :param      ax:   The ax to which render the field
        :type       ax:   { type_description }
        :param      itr:  The itr
        :type       itr:  float

        :returns:   No return
        :rtype:     None
        """
        if itr is not None:
            slice_number = self.itr_to_slice(itr=itr)

        elif slice_number is not None:
            itr = self.itr_list[slice_number]

        artist = Mesh(
            scalar=self.fields[slice_number],
            colormap=BKR,
        )

        ax.add_artist(artist)

        ax.title = f"{itr = :.4f}, {slice_number = }"
        ax.equal = True

        ax.colorbar = ColorBar(symmetric=True, position='right')

    def render_effective_index(self, ax: Axis) -> None:
        """
        Render the mode field evaluated at a certain itr and on a given ax.

        :param      ax:   The ax to which render the field
        :type       ax:   { type_description }

        :returns:   No return
        :rtype:     None
        """
        artist = Line(
            x=self.itr_list,
            y=self.effective_index,
        )

        ax.add_artist(artist)

        ax.x_label = 'Inverse taper ratio'
        ax.y_label = 'Effective index'

    def render_eigen_value(self, ax: Axis) -> None:
        """
        Render the mode field evaluated at a certain itr and on a given ax.

        :param      ax:   The ax to which render the field
        :type       ax:   { type_description }

        :returns:   No return
        :rtype:     None
        """

        artist = Line(
            x=self.itr_list,
            y=self.eigen_values.real,
        )

        ax.add_artist(artist)

        ax.x_label = 'Inverse taper ratio'
        ax.y_label = 'Effective index'

        # if ax is None:
        #     figure, ax = plt.subplots(1, 1, figsize=(12, 8))

        # ax.plot(self.itr_list, self.eigen_values.real)
        # ax.set_xlabel('Inverse taper ratio')
        # ax.set_ylabel('Eigen values')
        # if self.save_fittings:
        #     ax.plot(self.itr_list, self.fit_values.real, 'x')

    def plot(self, itr: float = None, slice_number: float = None, show_field: bool = True, show_index: bool = False, show_eigen_value: bool = False) -> None:
        """
        Render the mode field at given itr plus the effective index associated eigen value.

        :param      itr:  The itr
        :type       itr:  float

        :returns:   { description_of_the_return_value }
        :rtype:     None
        """
        figure = Scene2D(unit_size=(3, 3))

        plot_number = 0
        if show_field:
            ax = Axis(row=0, col=plot_number)
            figure.add_axes(ax)
            self.render_field(ax=ax, itr=itr, slice_number=slice_number)
            plot_number += 1

        if show_index:
            ax = Axis(row=0, col=plot_number)
            figure.add_axes(ax)
            self.render_effective_index(ax=ax)
            plot_number += 1

        if show_eigen_value:
            ax = Axis(row=0, col=plot_number)
            figure.add_axes(ax)
            self.render_eigen_value(ax=ax)
            plot_number += 1

        return figure
