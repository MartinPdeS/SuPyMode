cmake_minimum_required(VERSION 3.16)

file(STRINGS SuPyMode/VERSION SUPYMODE_VERSION)

project(SuPyMode VERSION "${SUPYMODE_VERSION}" LANGUAGES CXX)

set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_STANDARD 17)

set(CMAKE_VERBOSE_MAKEFILE OFF)
set(SUPYMODE_CXX_DIR "SuPyMode/cpp")
set(SUPYMODE_BIN_DIR "${CMAKE_SOURCE_DIR}/SuPyMode/binary")
add_compile_options("-w" "-Wdeprecated-declarations")

set(CMAKE_BUILD_TYPE "Release" CACHE STRING "Type of build" FORCE)

#FIND_PACKAGE------------------------
find_package(Python ${PYBIND11_PYTHON_VERSION} EXACT COMPONENTS Interpreter Development.Module)
add_subdirectory("extern/pybind11")
include_directories(${CMAKE_CURRENT_SOURCE_DIR}/extern/eigen)
include_directories("${CMAKE_SOURCE_DIR}/SuPyMode/cpp/includes")

if (WIN32)
    message("Windows OS detected")
    set(STATIC_STDLIB "-w -lm -static-libgcc -static-libstdc++ -Wl,-Bstatic,--whole-archive -lwinpthread -Wl,--no-whole-archive") # Need to keep -lwinpthread static because of MinGW dll directory not the same as windows' MSVC
elseif (CMAKE_CXX_COMPILER_ID STREQUAL "AppleClang")
    message("AppleClang compiler detected")
    set(STATIC_STDLIB "-w -lm -static-libstdc++")
elseif (CMAKE_CXX_COMPILER_ID STREQUAL "GNU")
    message("GNU compiler detected")
    # set(STATIC_STDLIB "-w -lm -static-libstdc++ -static-libgcc")
    set(STATIC_STDLIB "-w -lm") # <iostream> crash on github action when static-linking c++
endif()

message(STATUS "SuPyMode  version is: " ${SUPYMODE_VERSION})
message(STATUS "SuPyMode includes directory is: " ${SUPYMODE_CXX_DIR})
message(STATUS "Python version to being compiled against: " ${PYBIND11_PYTHON_VERSION})
message(STATUS "Binary will be installed in location: " ${SUPYMODE_BIN_DIR})
message(STATUS "Custom flags passed to linker: " ${STATIC_STDLIB})

#Building EigenSolver----------------------------------------------
pybind11_add_module(EigenSolver MODULE "${SUPYMODE_CXX_DIR}/eigensolver_interface.cpp")
set_target_properties(EigenSolver PROPERTIES OUTPUT_NAME "CppSolver")
target_link_libraries(EigenSolver PRIVATE "${STATIC_STDLIB}" "${Eigen_LIBRARIES}")

#Building SuperModes----------------------------------------------
pybind11_add_module(SuperMode MODULE "${SUPYMODE_CXX_DIR}/supermode_interface.cpp")
set_target_properties(SuperMode PROPERTIES OUTPUT_NAME "SuperMode")
target_link_libraries(SuperMode PRIVATE "${STATIC_STDLIB}")


#Building ModelParameters----------------------------------------------
pybind11_add_module(ModelParameters MODULE "${SUPYMODE_CXX_DIR}/model_parameters_interface.cpp")
set_target_properties(ModelParameters PROPERTIES OUTPUT_NAME "ModelParameters")
target_link_libraries(ModelParameters PRIVATE "${STATIC_STDLIB}")

set_target_properties(
    EigenSolver SuperMode ModelParameters
    PROPERTIES
    ARCHIVE_OUTPUT_DIRECTORY "${SUPYMODE_BIN_DIR}"
    LIBRARY_OUTPUT_DIRECTORY "${SUPYMODE_BIN_DIR}"
    RUNTIME_OUTPUT_DIRECTORY "${SUPYMODE_BIN_DIR}"
    POSITION_INDEPENDENT_CODE TRUE
)

install(
    TARGETS EigenSolver SuperMode ModelParameters
    DESTINATION "${SUPYMODE_BIN_DIR}"
    # RUNTIME_DEPENDENCIES # do not uncomment this as it will mixed dependencies that may not be manylinux-compatible
)


