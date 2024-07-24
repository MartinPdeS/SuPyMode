cmake_minimum_required(VERSION 3.16)

project(SuPyMode VERSION "${SUPYMODE_VERSION}" LANGUAGES CXX)

set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_VERBOSE_MAKEFILE OFF)
set(SUPYMODE_CXX_DIR "SuPyMode/cpp")
set(SUPYMODE_BIN_DIR "${CMAKE_SOURCE_DIR}/SuPyMode/binary")
set(CMAKE_BUILD_TYPE "Release" CACHE STRING "Type of build" FORCE)

# Ensure -fPIC is used globally
set(CMAKE_POSITION_INDEPENDENT_CODE ON)
add_compile_options(-fPIC)

#FIND_PACKAGE------------------------
find_package(Python REQUIRED COMPONENTS Interpreter Development.Module)
find_package(pybind11 CONFIG REQUIRED)

include_directories("${CMAKE_SOURCE_DIR}/SuPyMode/cpp/includes")
include_directories("${CMAKE_SOURCE_DIR}/SuPyMode/cpp/headers")

add_compile_options("-w" "-Wdeprecated-declarations")

# Platform-specific settings for static linking
if (WIN32)
    message("Windows OS detected")
    set(STATIC_LINK_OPTIONS "-static")
elseif (CMAKE_CXX_COMPILER_ID STREQUAL "AppleClang")
    message("AppleClang compiler detected")
elseif (CMAKE_CXX_COMPILER_ID STREQUAL "GNU")
    message("GNU compiler detected")
endif()

message(STATUS "SuPyMode  version is: " ${SUPYMODE_VERSION})
message(STATUS "SuPyMode includes directory is: " ${SUPYMODE_CXX_DIR})
message(STATUS "Python version to being compiled against: " ${PYBIND11_PYTHON_VERSION})
message(STATUS "Binary will be installed in location: " ${SUPYMODE_BIN_DIR})

#Building EigenSolver----------------------------------------------
pybind11_add_module(EigenSolver MODULE "${SUPYMODE_CXX_DIR}/eigensolver_interface.cpp")
set_target_properties(EigenSolver PROPERTIES OUTPUT_NAME "CppSolver")
target_link_libraries(EigenSolver PRIVATE "${Eigen_LIBRARIES}")
target_link_options(EigenSolver PUBLIC ${STATIC_LINK_OPTIONS})
target_include_directories(EigenSolver PRIVATE "${CMAKE_CURRENT_SOURCE_DIR}/extern/eigen")
target_include_directories(EigenSolver PRIVATE "${CMAKE_CURRENT_SOURCE_DIR}/extern/spectra/include")

#Building SuperModes----------------------------------------------
pybind11_add_module(SuperMode MODULE "${SUPYMODE_CXX_DIR}/supermode_interface.cpp")
set_target_properties(SuperMode PROPERTIES OUTPUT_NAME "SuperMode")
target_link_options(SuperMode PUBLIC ${STATIC_LINK_OPTIONS})
target_include_directories(SuperMode PRIVATE "${CMAKE_CURRENT_SOURCE_DIR}/extern/eigen")

#Building ModelParameters----------------------------------------------
pybind11_add_module(ModelParameters MODULE "${SUPYMODE_CXX_DIR}/model_parameters_interface.cpp")
set_target_properties(ModelParameters PROPERTIES OUTPUT_NAME "ModelParameters")
target_link_options(ModelParameters PUBLIC ${STATIC_LINK_OPTIONS})
target_include_directories(ModelParameters PRIVATE "${CMAKE_CURRENT_SOURCE_DIR}/extern/eigen")

#Building test example----------------------------------------------
pybind11_add_module(Example MODULE "${SUPYMODE_CXX_DIR}/includes/mesh.cpp")
set_target_properties(Example PROPERTIES OUTPUT_NAME "Example")
target_link_options(Example PUBLIC ${STATIC_LINK_OPTIONS})
target_include_directories(Example PRIVATE "${CMAKE_CURRENT_SOURCE_DIR}/extern/eigen")

install(
    TARGETS EigenSolver SuperMode ModelParameters
    DESTINATION "${SUPYMODE_BIN_DIR}"
    # RUNTIME_DEPENDENCIES # do not uncomment this as it will mixed dependencies that may not be manylinux-compatible
)

