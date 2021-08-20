

cmake_minimum_required(VERSION 3.0)

project(SuPyModes LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 17)

set(CMAKE_CXX_STANDARD_REQUIRED ON)

find_package(Python 3 COMPONENTS Interpreter Development EXACT)

add_subdirectory(extern/pybind11)


set(Boost_USE_STATIC_LIBS ON)
find_package(Boost COMPONENTS filesystem REQUIRED)
include_directories(${Boost_INCLUDE_DIR})



find_package(Eigen3 REQUIRED NO_MODULE)
set(EIGEN_INCLUDE_DIR ../extern/eigen)
include_directories(${EIGEN_INCLUDE_DIR})

find_package(PythonInterp 3 REQUIRED)

include_directories( ${PYTHON_INCLUDE_DIRS} ${PYBIND11_INCLUDE_DIRS} )

FIND_PACKAGE(Boost REQUIRED COMPONENTS system iostreams)


set(BIN_PATH ${PROJECT_BINARY_DIR}/../SuPyModes/includes/)
#set(CMAKE_LIBRARY_OUTPUT_DIRECTORY ${PROJECT_BINARY_DIR}/../SuPyModes/includes)
pybind11_add_module(EigenSolver MODULE ../SuPyModes/includes/interface.cpp )

include_directories(${Boost_INCLUDE_DIRS})
target_link_libraries(EigenSolver PRIVATE Eigen3::Eigen ${Boost_LIBRARIES})
target_compile_options (EigenSolver PRIVATE -std=c++17 -O3)


message(STATUS "Installing file in ${BIN_PATH}" )
install(TARGETS EigenSolver DESTINATION ${BIN_PATH} CONFIGURATIONS Release)
install(TARGETS EigenSolver DESTINATION ${BIN_PATH} CONFIGURATIONS Debug)
