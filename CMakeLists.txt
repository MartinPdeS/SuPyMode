cmake_minimum_required(VERSION 3.0)

macro(RemoveDebugCXXFlag flag)
    string(REPLACE "${flag}" "" CMAKE_CXX_FLAGS_DEBUG "${CMAKE_CXX_FLAGS_DEBUG}")
endmacro()

macro(RemoveReleaseCXXFlag flag)
    string(REPLACE "${flag}" "" CMAKE_CXX_FLAGS_RELEASE "${CMAKE_CXX_FLAGS_RELEASE}")
endmacro()

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
pybind11_add_module(EigenSolver MODULE ../SuPyModes/includes/interface.cpp )

include_directories(${Boost_INCLUDE_DIRS})
target_link_libraries(EigenSolver PRIVATE Eigen3::Eigen ${Boost_LIBRARIES})
target_compile_options (EigenSolver PRIVATE -O2)

message("${CMAKE_CXX_FLAGS_RELEASE} \n\n\n")
RemoveReleaseCXXFlag("/RTC1")
message("${CMAKE_CXX_FLAGS_RELEASE} \n\n\n")


message("${CMAKE_CXX_FLAGS_DEBUG} \n\n\n")
RemoveDebugCXXFlag("/RTC1")
message("${CMAKE_CXX_FLAGS_DEBUG} \n\n\n")

message(STATUS "Installing file in ${BIN_PATH}" )

install(TARGETS EigenSolver DESTINATION ${BIN_PATH} CONFIGURATIONS Release)
install(TARGETS EigenSolver DESTINATION ${BIN_PATH} CONFIGURATIONS Debug)


#Upload Pypi package----------------------------------------------
add_custom_command(
  OUTPUT ${CMAKE_CURRENT_SOURCE_DIR}/UploadPypi.txt
  COMMAND python3.8 -m twine upload --password ${PipPassword} --username ${PipUsername} --repository pypi ./dist/*
  COMMENT "Upload on Pypi")

add_custom_target(UploadPypi DEPENDS ${CMAKE_CURRENT_SOURCE_DIR}/UploadPypi.txt)


#Cleaning command----------------------------------------------
add_custom_command(
  OUTPUT ${CMAKE_CURRENT_SOURCE_DIR}/Clean.txt
  COMMAND rm -rf CMakeCache.txt cmake_install.cmake CMakeFiles
  COMMENT "Cleaning cmake output files")

add_custom_target(Clean DEPENDS ${CMAKE_CURRENT_SOURCE_DIR}/Clean.txt)