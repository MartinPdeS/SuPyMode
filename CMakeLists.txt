cmake_minimum_required(VERSION 3.0)

project(SuPyModes LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 17)

set(CMAKE_CXX_STANDARD_REQUIRED ON)

find_package(Python 3 COMPONENTS Interpreter Development EXACT)

add_subdirectory(extern/pybind11)

find_package (Eigen3 REQUIRED NO_MODULE)
include_directories(${EIGEN_INCLUDE_DIR})

FIND_PACKAGE(PythonInterp 3 REQUIRED)

INCLUDE_DIRECTORIES( ${PYTHON_INCLUDE_DIRS} ${PYBIND11_INCLUDE_DIRS} )


add_custom_command(
  OUTPUT ${CMAKE_CURRENT_SOURCE_DIR}/Eigen.txt
  COMMAND rm -rf extern/eigen/build && mkdir extern/eigen/build && cd extern/eigen/build && cmake .. && make install && cd ../..
  COMMENT "Installing Eigen3 c++ library")

add_custom_target(installEigen ALL DEPENDS ${CMAKE_CURRENT_SOURCE_DIR}/Eigen.txt)


add_custom_command(
  OUTPUT ${CMAKE_CURRENT_SOURCE_DIR}/Spectra.txt
  COMMAND rm -rf extern/spectra/build && mkdir extern/spectra/build && cd extern/spectra/build && cmake .. && make install && cd ../..
  COMMENT "Installing Spectra c++ library")

add_custom_target(installSpectra ALL DEPENDS ${CMAKE_CURRENT_SOURCE_DIR}/Spectra.txt)


add_custom_command(
  OUTPUT ${CMAKE_CURRENT_SOURCE_DIR}/Dependencies.txt
  COMMAND make installEigen && make installSpectra
  COMMENT "Installing dependencies")

add_custom_target(installDependence ALL DEPENDS ${CMAKE_CURRENT_SOURCE_DIR}/Dependencies.txt)


add_custom_command(
  OUTPUT ${CMAKE_CURRENT_SOURCE_DIR}/Upload.txt
  COMMAND python3.8 -m twine upload --password ${PipPassword} --username ${PipUsername} --repository pypi ./dist/*
  COMMENT "Upload on Pypi")

add_custom_target(Upload ALL DEPENDS ${CMAKE_CURRENT_SOURCE_DIR}/Upload.txt)

FIND_PACKAGE(Boost REQUIRED COMPONENTS system iostreams)

set(CMAKE_LIBRARY_OUTPUT_DIRECTORY ${PROJECT_BINARY_DIR}/SuPyModes/includes)
pybind11_add_module(EigenSolver MODULE SuPyModes/includes/interface.cpp )

include_directories(${Boost_INCLUDE_DIRS})
target_link_libraries(EigenSolver PRIVATE Eigen3::Eigen ${Boost_LIBRARIES})
target_compile_options (EigenSolver PRIVATE -std=c++17 -O3)
