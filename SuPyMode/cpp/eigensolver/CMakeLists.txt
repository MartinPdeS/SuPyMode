# cpp/eigensolver/CMakeLists.txt

set(NAME "eigensolver")

# Create a shared library for functionality.
add_library("${NAME}" STATIC "${NAME}.cpp")

target_link_libraries("${NAME}" PUBLIC pybind11::module supermode)

target_include_directories("${NAME}" PUBLIC "${CMAKE_SOURCE_DIR}/extern/spectra/include")

# Create a Python module, if needed.
pybind11_add_module("interface_${NAME}" MODULE interface.cpp)
set_target_properties("interface_${NAME}" PROPERTIES OUTPUT_NAME "interface_${NAME}")
target_link_libraries("interface_${NAME}" PRIVATE "${NAME}")

set_property(GLOBAL APPEND PROPERTY SUPYMODE_TARGETS "${NAME}" "interface_${NAME}")
