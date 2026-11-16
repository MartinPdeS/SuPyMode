#pragma once
#include <pybind11/pybind11.h>
#include <string>

/** Validate through TypedUnit, then extract canonical numerical magnitudes.
 * Bare dimensional values retain the legacy SI convention. rank=-1 accepts
 * any shape; rank=0,1,2 selects Scalar, Vector, or Matrix respectively.
 */
pybind11::object numerical_input(pybind11::handle value, const std::string &unit,
                                 const std::string &name, int rank = -1, bool positive = false,
                                 bool finite = true);
