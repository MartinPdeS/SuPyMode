#include <pybind11/pybind11.h>
#include <pybind11/numpy.h>
#include <pybind11/stl.h>

#include "taper.h"

namespace py = pybind11;

static py::array make_vector_view_1d(const std::vector<double>& values, py::handle base_object) {
    return py::array(
        py::buffer_info(
            const_cast<double*>(values.data()),
            static_cast<py::ssize_t>(sizeof(double)),
            py::format_descriptor<double>::format(),
            1,
            { static_cast<py::ssize_t>(values.size()) },
            { static_cast<py::ssize_t>(sizeof(double)) }
        ),
        base_object
    );
}

PYBIND11_MODULE(interface_taper, module_handle) {
    module_handle.doc() = R"pdoc(
        C++ implementation of an optical fiber taper alpha profile.

        Notes
        -----
        Some properties return NumPy arrays as zero copy views into the internal C++ buffers.
        The returned arrays remain valid as long as the owning Python object is alive.
        Calling initialize again can invalidate previously returned views.
        )pdoc";

    py::class_<Interpolator1D>(module_handle, "Interpolator1D", R"pdoc(
        One dimensional linear interpolator.

        This object behaves like a callable in Python and can be used as a drop in
        replacement for SciPy interp1d in typical workflows.

        Notes
        -----
        This implementation copies x and y into the interpolator object for robust lifetime handling.
        )pdoc")
        .def(
            py::init<std::vector<double>, std::vector<double>, bool, double>(),
            py::arg("x"),
            py::arg("y"),
            py::arg("bounds_error") = false,
            py::arg("fill_value") = 0.0,
            R"pdoc(
                Construct an interpolator from x and y samples.

                Parameters
                ----------
                x : list[float] or numpy.ndarray
                    Sorted x sample positions.
                y : list[float] or numpy.ndarray
                    y values corresponding to x.
                bounds_error : bool, optional
                    If True, out of bounds evaluation raises.
                fill_value : float, optional
                    Value returned for out of bounds queries when bounds_error is False.
                )pdoc"
        )
        .def(
            "__call__",
            [](const Interpolator1D& self, py::object x_like) -> py::object {
                if (!py::isinstance<py::array>(x_like) && !py::isinstance<py::sequence>(x_like)) {
                    const double x = py::cast<double>(x_like);
                    return py::float_(self(x));
                }

                py::array x_array = py::array::ensure(x_like);
                if (!x_array) throw std::invalid_argument("x must be a scalar or array-like.");

                py::array_t<double, py::array::c_style | py::array::forcecast> x_double(x_array);
                auto x_buffer = x_double.request();

                py::array_t<double> out(static_cast<py::ssize_t>(x_buffer.size));
                auto out_buffer = out.request();

                const auto* x_ptr = static_cast<const double*>(x_buffer.ptr);
                auto* out_ptr = static_cast<double*>(out_buffer.ptr);

                for (py::ssize_t i = 0; i < x_buffer.size; ++i) {
                    out_ptr[i] = self(x_ptr[i]);
                }

                const auto input_buffer = x_double.request();
                std::vector<py::ssize_t> output_shape(
                    input_buffer.shape.begin(),
                    input_buffer.shape.end()
                );

                out.resize(output_shape);
                return out;
            },
            py::arg("x"),
            R"pdoc(
                Evaluate interpolation.

                Parameters
                ----------
                x : float or numpy.ndarray
                    Query points.

                Returns
                -------
                float or numpy.ndarray
                    Interpolated values.
                )pdoc"
        )
        .def_property_readonly(
            "bounds_error",
            &Interpolator1D::bounds_error,
            R"pdoc(
                Whether out of bounds evaluation raises.

                Parameters
                ----------
                bounds_error : bool
                    Flag indicating whether out of bounds evaluation raises.

                Returns
                -------
                bool
                    bounds_error flag.
                )pdoc"
        )
        .def_property_readonly(
            "fill_value",
            &Interpolator1D::fill_value,
            R"pdoc(
                Fill value used when bounds_error is False.

                Parameters
                ----------
                fill_value : float
                    Fill value.

                Returns
                -------
                float
                    Fill value.
                )pdoc"
        )
        .def_property_readonly(
            "domain",
            [](const Interpolator1D& self) {
                return py::make_tuple(self.x_min(), self.x_max());
            },
            R"pdoc(
                Domain of the interpolator.

                Parameters
                ----------
                domain : tuple[float, float]
                    (x_min, x_max) bounds.

                Returns
                -------
                tuple[float, float]
                    Domain bounds.
                )pdoc"
        );

    py::class_<TaperSection>(module_handle, "TaperSection",
            R"pdoc(
                A single taper section defined by sampled z positions and corresponding radii.

                Notes
                -----
                This binding currently exposes scalar convenience properties.
            )pdoc"
        )
        .def(py::init<std::vector<double>, std::vector<double>, double, double>(),
            py::arg("z_array"),
            py::arg("radius_array"),
            py::arg("heating_length_initial") = std::numeric_limits<double>::quiet_NaN(),
            py::arg("heating_length_final") = std::numeric_limits<double>::quiet_NaN(),
            R"pdoc(
                Construct a taper section.

                Parameters
                ----------
                z_array : List[float]
                    List of z coordinates defining the section.
                radius_array : List[float]
                    List of radius values corresponding to each z coordinate.
                heating_length_initial : float, optional
                    Heating length at the start of the section. May be nan if not set.
                heating_length_final : float, optional
                    Heating length at the end of the section. May be nan if not set.

                Returns
                -------
                TaperSection
                    Newly created taper section object.

                Raises
                ------
                ValueError
                    If input conditions are invalid (e.g., size mismatch, empty arrays, unsorted z_array).
                )pdoc"
        )
        .def_property_readonly(
            "z_initial",
            &TaperSection::z_initial,
            R"pdoc(
                First z coordinate of the section.

                Parameters
                ----------
                z_initial : float
                    First z coordinate of the section.

                Returns
                -------
                float
                    z_array[0].
                )pdoc"
        )
        .def_property_readonly(
            "z_final",
            &TaperSection::z_final,
            R"pdoc(
                Last z coordinate of the section.

                Parameters
                ----------
                z_final : float
                    Last z coordinate of the section.

                Returns
                -------
                float
                    z_array[-1].
                )pdoc"
        )
        .def_property_readonly(
            "radius_initial",
            &TaperSection::radius_initial,
            R"pdoc(
                Radius at the start of the section.

                Parameters
                ----------
                radius_initial : float
                    Radius at the start of the section.

                Returns
                -------
                float
                    radius_array[0].
                )pdoc"
        )
        .def_property_readonly(
            "radius_final",
            &TaperSection::radius_final,
            R"pdoc(
                Radius at the end of the section.

                Parameters
                ----------
                radius_final : float
                    Radius at the end of the section.

                Returns
                -------
                float
                    radius_array[-1].
                )pdoc"
        )
        .def_property_readonly(
            "is_constant",
            &TaperSection::is_constant,
            R"pdoc(
                Whether the section radius is constant.

                Parameters
                ----------
                is_constant : bool
                    True if radius_array[0] equals radius_array[-1].

                Returns
                -------
                bool
                    Constant radius flag.
                )pdoc"
        )
        .def_property_readonly(
            "heating_length_initial",
            &TaperSection::heating_length_initial,
            R"pdoc(
                Heating length at the beginning of the section.

                Parameters
                ----------
                heating_length_initial : float
                    Heating length at the beginning of the section.

                Returns
                -------
                float
                    Heating length value. May be nan when not set.
                )pdoc"
        )
        .def_property_readonly(
            "heating_length_final",
            &TaperSection::heating_length_final,
            R"pdoc(
                Heating length at the end of the section.

                Parameters
                ----------
                heating_length_final : float
                    Heating length at the end of the section.

                Returns
                -------
                float
                    Heating length value. May be nan when not set.
                )pdoc"
        );

    py::class_<AlphaProfile>(module_handle, "AlphaProfile",
        R"pdoc(
            Taper profile built from multiple taper sections and sampled onto a uniform grid.

            The workflow is:
            1) Add one or more sections using add_taper_segment and or add_constant_segment.
            2) Call initialize to assemble the full profile arrays.
            3) Access scalar values and numpy views of the assembled arrays.

            Notes
            -----
            The numpy properties return zero copy views. They are invalidated if you call initialize again.
        )pdoc")
        .def(
            py::init<double, std::size_t, bool, std::string, bool>(),
            py::arg("initial_radius") = 1.0,
            py::arg("n_point") = 200,
            py::arg("symmetric") = false,
            py::arg("label") = "profile",
            py::arg("add_end_of_taper_section") = true,
            R"pdoc(
                Construct an alpha profile.

                Parameters
                ----------
                initial_radius : float, optional
                    Reference radius used to define the inverse taper ratio.
                n_point : int, optional
                    Number of points used when sampling the assembled profile in initialize.
                symmetric : bool, optional
                    If True, mirror the profile around the waist after initialization.
                label : str, optional
                    Label stored for identification.
                add_end_of_taper_section : bool, optional
                    If True, append a constant section at the end when the last taper section is not constant.

                Returns
                -------
                AlphaProfile
                    Newly created profile object.
            )pdoc"
        )

        .def(
            "add_taper_segment",
            &AlphaProfile::add_taper_segment,
            py::arg("alpha"),
            py::arg("initial_heating_length"),
            py::arg("stretching_length"),
            py::arg("n_point") = 100,
            R"pdoc(
                Append a tapered section following the last defined section.

                Parameters
                ----------
                alpha : float
                    Alpha parameter controlling how heating length evolves with stretching.
                initial_heating_length : float
                    Heating length at the start of the segment.
                stretching_length : float
                    Stretching length applied for this segment.
                n_point : int, optional
                    Number of discretization points used for this segment.

                Returns
                -------
                None

                Raises
                ------
                ValueError
                    If input conditions are invalid.
            )pdoc"
        )

        .def(
            "add_constant_segment",
            &AlphaProfile::add_constant_segment,
            py::arg("length"),
            py::arg("n_point") = 100,
            R"pdoc(
                Append a constant radius segment to the profile.

                Parameters
                ----------
                length : float
                    Length of the constant segment.
                n_point : int, optional
                    Number of discretization points for the constant segment.

                Returns
                -------
                None
            )pdoc"
        )

        .def(
            "add_end_of_taper_segment",
            &AlphaProfile::add_end_of_taper_segment,
            py::arg("n_point") = 100,
            R"pdoc(
                Append a constant segment at the end if the last section is not constant.

                Parameters
                ----------
                n_point : int, optional
                    Number of discretization points for the constant segment.

                Returns
                -------
                None

                Raises
                ------
                RuntimeError
                    If the last section has no finite heating_length_final value.
            )pdoc"
        )
        .def(
            "initialize",
            &AlphaProfile::initialize,
            R"pdoc(
                Assemble the full profile arrays and compute derived quantities.

                This samples the piecewise section definition onto a uniform grid of n_point points,
                then computes:
                - radius
                - inverse taper ratio (ITR)
                - adiabatic factor
                - taper angle

                If symmetric is True, all arrays are mirrored without repeating the waist point.

                Returns
                -------
                None

                Raises
                ------
                RuntimeError
                    If the profile has no length or has not been properly defined.
            )pdoc"
        )
        .def_property_readonly(
            "n_point",
            &AlphaProfile::n_point,
            R"pdoc(
                Number of sampling points used in initialize.

                Parameters
                ----------
                n_point : int
                    Number of sampling points used in :meth:`initialize`.

                Returns
                -------
                int
                    Sampling count.
            )pdoc"
        )
        .def_property_readonly(
            "initial_radius",
            &AlphaProfile::initial_radius,
            R"pdoc(
                Reference radius used to define the inverse taper ratio.

                Parameters
                ----------
                initial_radius : float
                    Reference radius used to define ITR.

                Returns
                -------
                float
                    Reference radius.
            )pdoc"
        )
        .def_property_readonly(
            "symmetric",
            &AlphaProfile::symmetric,
            R"pdoc(
                Whether the assembled profile is mirrored.

                Parameters
                ----------
                symmetric : bool
                    Mirroring flag.

                Returns
                -------
                bool
                    Mirroring flag.
            )pdoc"
        )
        .def_property_readonly(
            "label",
            &AlphaProfile::label,
            R"pdoc(
                Label stored for identification.

                Parameters
                ----------
                label : str
                    Profile label.

                Returns
                -------
                str
                    Profile label.
            )pdoc"
        )
        .def_property_readonly(
            "last_z",
            &AlphaProfile::last_z,
            R"pdoc(
                End coordinate of the piecewise section definition.

                Parameters
                ----------
                last_z : float
                    End coordinate of the piecewise section definition.

                Returns
                -------
                float
                    Last z value, or 0 if no sections exist.
            )pdoc"
        )
        .def_property_readonly(
            "total_length",
            &AlphaProfile::total_length,
            R"pdoc(
                Total length of the piecewise section definition.

                Parameters
                ----------
                total_length : float
                    Total length of the piecewise section definition.

                Returns
                -------
                float
                    Total length, or 0 if no sections exist.
            )pdoc"
        )
        .def_property_readonly(
            "last_radius",
            &AlphaProfile::last_radius,
            R"pdoc(
                Radius at the end of the piecewise section definition.

                Parameters
                ----------
                last_radius : float
                    Radius at the end of the piecewise section definition.

                Returns
                -------
                float
                    Last radius, or initial_radius if no sections exist.
            )pdoc"
        )
        .def_property_readonly(
            "smallest_itr",
            &AlphaProfile::smallest_itr,
            R"pdoc(
                Smallest inverse taper ratio after initialization.

                Parameters
                ----------
                smallest_itr : float
                    Smallest inverse taper ratio after initialization.

                Returns
                -------
                float
                    Minimum value of itr_list.

                Raises
                ------
                RuntimeError
                    If the profile has not been initialized.
            )pdoc"
        )

        .def_property_readonly(
            "distance",
            [](py::object self_object) {
                const auto& self_reference = self_object.cast<const AlphaProfile&>();
                return make_vector_view_1d(self_reference.distance(), self_object);
            },
            R"pdoc(
                Sampled distance coordinate array as a zero copy numpy view.

                Parameters
                ----------
                distance_numpy : numpy.ndarray
                    Sampled distance coordinate array.

                Returns
                -------
                numpy.ndarray
                    1D array of shape (n_point,) if symmetric is False,
                    or (2*n_point-1,) if symmetric is True.

                Raises
                ------
                RuntimeError
                    If the profile has not been initialized.

                Notes
                -----
                This is a zero copy view into internal C++ memory.
            )pdoc"
        )
        .def_property_readonly(
            "radius",
            [](py::object self_object) {
                const auto& self_reference = self_object.cast<const AlphaProfile&>();
                return make_vector_view_1d(self_reference.radius(), self_object);
            },
            R"pdoc(
                Sampled radius array as a zero copy numpy view.

                Returns
                -------
                numpy.ndarray
                    1D array with the same shape as distance.

                Raises
                ------
                RuntimeError
                    If the profile has not been initialized.

                Notes
                -----
                This is a zero copy view into internal C++ memory.
            )pdoc"
        )
        .def_property_readonly(
            "itr_list",
            [](py::object self_object) {
                const auto& self_reference = self_object.cast<const AlphaProfile&>();
                return make_vector_view_1d(self_reference.itr_list(), self_object);
            },
            R"pdoc(
                Inverse taper ratio array as a zero copy numpy view.

                Returns
                -------
                numpy.ndarray
                    1D array with the same shape as distance_numpy.

                Raises
                ------
                RuntimeError
                    If the profile has not been initialized.

                Notes
                -----
                This is a zero copy view into internal C++ memory.
            )pdoc"
        )
        .def_property_readonly(
            "adiabatic",
            [](py::object self_object) {
                const auto& self_reference = self_object.cast<const AlphaProfile&>();
                return make_vector_view_1d(self_reference.adiabatic(), self_object);
            },
            R"pdoc(
                Adiabatic factor along the taper as a zero copy numpy view.

                Returns
                -------
                numpy.ndarray
                    1D array with the same shape as distance_numpy.

                Raises
                ------
                RuntimeError
                    If the profile has not been initialized.

                Notes
                -----
                This is a zero copy view into internal C++ memory.
            )pdoc"
        )
        .def_property_readonly(
            "taper_angle",
            [](py::object self_object) {
                const auto& self_reference = self_object.cast<const AlphaProfile&>();
                return make_vector_view_1d(self_reference.taper_angle(), self_object);
            },
            R"pdoc(
                Taper angle proxy along the taper as a zero copy numpy view.

                Returns
                -------
                numpy.ndarray
                    1D array with the same shape as distance_numpy.

                Raises
                ------
                RuntimeError
                    If the profile has not been initialized.

                Notes
                -----
                This is a zero copy view into internal C++ memory.
            )pdoc"
        )
        .def(
            "evaluate_adiabatic_factor",
            [](AlphaProfile& profile, py::array_t<double, py::array::c_style | py::array::forcecast> itr_array) {
                const auto buffer = itr_array.request();
                const auto* data_ptr = static_cast<const double*>(buffer.ptr);

                std::vector<double> itr_values(static_cast<std::size_t>(buffer.size));
                std::copy(data_ptr, data_ptr + buffer.size, itr_values.begin());

                auto result_values = profile.evaluate_adiabatic_factor(itr_values);

                py::array_t<double> result(static_cast<py::ssize_t>(result_values.size()));
                auto result_buffer = result.request();
                auto* out_ptr = static_cast<double*>(result_buffer.ptr);
                std::copy(result_values.begin(), result_values.end(), out_ptr);
                return result;
            },
            py::arg("itr"),
            R"pdoc(
                Evaluate the adiabatic factor for a given inverse taper ratio array.

                Returns
                -------
                numpy.ndarray
                    1D array of adiabatic factor values. Values outside the interpolation
                    domain return nan.

                Raises
                ------
                RuntimeError
                    If the profile has not been initialized.
            )pdoc"
        )
        .def(
            "evaluate_distance_vs_itr",
            [](AlphaProfile& profile, py::array_t<double, py::array::c_style | py::array::forcecast> distance_array) {
                const auto buffer = distance_array.request();
                const auto* data_ptr = static_cast<const double*>(buffer.ptr);

                std::vector<double> distance_values(static_cast<std::size_t>(buffer.size));
                std::copy(data_ptr, data_ptr + buffer.size, distance_values.begin());

                auto result_values = profile.evaluate_distance_vs_itr(distance_values);

                py::array_t<double> result(static_cast<py::ssize_t>(result_values.size()));
                auto result_buffer = result.request();
                auto* out_ptr = static_cast<double*>(result_buffer.ptr);
                std::copy(result_values.begin(), result_values.end(), out_ptr);
                return result;
            },
            py::arg("distance"),
            R"pdoc(
                Evaluate the inverse taper ratio for a given distance array.

                Returns
                -------
                numpy.ndarray
                    1D array of inverse taper ratio values.

                Raises
                ------
                RuntimeError
                    If the profile has not been initialized.
                IndexError
                    If any distance value is outside the interpolation range.
                )pdoc"
        )
        .def(
            "get_itr_vs_distance_interpolation",
            [](const AlphaProfile& profile) {
                return profile.get_itr_vs_distance_interpolation_object();
            },
            R"pdoc(
                Return an interpolator mapping distance to inverse taper ratio.

                Parameters
                ----------
                None

                Returns
                -------
                Interpolator1D
                    Callable interpolator f(z) -> itr.
            )pdoc"
        )
        .def(
            "get_distance_vs_itr_interpolation",
            [](const AlphaProfile& profile) {
                return profile.get_distance_vs_itr_interpolation_object();
            },
            R"pdoc(
                Return an interpolator mapping inverse taper ratio to distance.

                Parameters
                ----------
                None

                Returns
                -------
                Interpolator1D
                    Callable interpolator f(itr) -> distance.
            )pdoc"
        )
        ;
}
