# cpp/boundaries/CMakeLists.txt

set(NAME "boundaries")

# Create a shared library for functionality.
add_library("${NAME}_lib" STATIC "${NAME}.cpp")

target_link_libraries("${NAME}_lib" PRIVATE pybind11::module)

# Create a Python module, if needed.
pybind11_add_module("${NAME}" MODULE interface.cpp)
target_link_libraries("${NAME}" PRIVATE "${NAME}_lib")

set_property(GLOBAL APPEND PROPERTY SUPYMODE_TARGETS "${NAME}")
