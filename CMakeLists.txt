cmake_minimum_required(VERSION 3.20)
set(CMAKE_OSX_DEPLOYMENT_TARGET "15.0" CACHE STRING "Minimum OS X deployment version")
project(SuPyMode VERSION "${SUPYMODE_VERSION}" LANGUAGES CXX)

# CMake settings
set(CMAKE_VERBOSE_MAKEFILE ON)
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_POSITION_INDEPENDENT_CODE ON)
set(CMAKE_BUILD_TYPE "Release" CACHE STRING "Type of build")

# Set the output directory for libraries
set(LOCAL_CXX_DIR "${PROJECT_NAME}/cpp")
set(LOCAL_BIN_DIR "${CMAKE_SOURCE_DIR}/${PROJECT_NAME}/binary")

# Find dependencies
if(APPLE AND CMAKE_CXX_COMPILER_ID STREQUAL "AppleClang")
    set(OpenMP_CXX_FLAGS "-Xpreprocessor -fopenmp" CACHE STRING "")
    set(OpenMP_CXX_LIB_NAMES "omp" CACHE STRING "")
    find_library(OpenMP_omp_LIBRARY omp HINTS $ENV{LIBOMP_PREFIX}/lib /opt/homebrew/opt/libomp/lib /usr/local/opt/libomp/lib)
    if(OpenMP_omp_LIBRARY)
        set(OpenMP_CXX_LIBRARIES "${OpenMP_omp_LIBRARY}" CACHE STRING "")
    endif()
endif()

find_package(OpenMP REQUIRED)
find_package(Python REQUIRED COMPONENTS Interpreter Development.Module)
find_package(pybind11 CONFIG REQUIRED)

# Platform-specific settings for static linking
if (WIN32 AND CMAKE_CXX_COMPILER_ID STREQUAL "GNU")
    message("MinGW detected on Windows")
    set(STATIC_LINK_OPTIONS "-static")
    add_compile_options(-fopenmp)
    add_link_options(-static -fopenmp -Wl,--whole-archive -lgomp -Wl,--no-whole-archive)
endif()

# Compiler and linker options
add_compile_options(-Wall -Wextra -pedantic-errors)

# Print some messages
message(STATUS "OPENMP flags: ${OpenMP_CXX_FLAGS}")
message(STATUS "SuPyMode version is: ${SUPYMODE_VERSION}")
message(STATUS "SuPyMode includes directory is: ${SUPYMODE_CXX_DIR}")
message(STATUS "Python version to be compiled against: ${PYBIND11_PYTHON_VERSION}")
message(STATUS "Binary will be installed in location: ${LOCAL_BIN_DIR}")

add_subdirectory(SuPyMode/cpp/mesh)                  # mesh
add_subdirectory(SuPyMode/cpp/model_parameters)      # model_parameters
add_subdirectory(SuPyMode/cpp/supermode)             # supermode
add_subdirectory(SuPyMode/cpp/eigensolver)           # eigen_solver
