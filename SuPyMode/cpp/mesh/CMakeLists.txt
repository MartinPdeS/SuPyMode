# cpp/mesh/CMakeLists.txt

set(NAME "mesh")

# Create a shared library for functionality.
add_library("${NAME}" STATIC "${NAME}.cpp")

target_include_directories("${NAME}" PUBLIC "${CMAKE_SOURCE_DIR}/extern/eigen")

target_link_libraries("${NAME}" PRIVATE pybind11::module)

# Create a Python module, if needed.
pybind11_add_module("interface_${NAME}" MODULE interface.cpp)
set_target_properties("interface_${NAME}" PROPERTIES OUTPUT_NAME "interface_${NAME}")
target_link_libraries("interface_${NAME}" PRIVATE "${NAME}")

set_property(GLOBAL APPEND PROPERTY SUPYMODE_TARGETS "${NAME}" "interface_${NAME}")
