cmake_minimum_required(VERSION 3.0)

message("PYTHON VERSION=${VERSION}")

set(Token    "$ENV{PyPiUsername}")
set(Password "$ENV{PyPiPassword}")


macro(RemoveDebugCXXFlag flag)
    string(REPLACE "${flag}" "" CMAKE_CXX_FLAGS_DEBUG "${CMAKE_CXX_FLAGS_DEBUG}")
endmacro()

macro(RemoveReleaseCXXFlag flag)
    string(REPLACE "${flag}" "" CMAKE_CXX_FLAGS_RELEASE "${CMAKE_CXX_FLAGS_RELEASE}")
endmacro()

if (WIN32)

    # if windows is used
    set(CMAKE_SYSTEM_NAME Windows)
    set(TOOLCHAIN_PREFIX clang)

    # cross compilers to use for C, C++ and Fortran
    set(CMAKE_C_COMPILER ${TOOLCHAIN_PREFIX}-gcc)
    set(CMAKE_CXX_COMPILER ${TOOLCHAIN_PREFIX}-g++)
    set(CMAKE_Fortran_COMPILER ${TOOLCHAIN_PREFIX}-gfortran)
    set(CMAKE_RC_COMPILER ${TOOLCHAIN_PREFIX}-windres)

    # target environment on the build host system
    set(CMAKE_FIND_ROOT_PATH /usr/${TOOLCHAIN_PREFIX})

    # modify default behavior of FIND_XXX() commands
    set(CMAKE_FIND_ROOT_PATH_MODE_PROGRAM NEVER)
    set(CMAKE_FIND_ROOT_PATH_MODE_LIBRARY ONLY)
    set(CMAKE_FIND_ROOT_PATH_MODE_INCLUDE ONLY)

    RemoveReleaseCXXFlag("/RTC1")
    RemoveDebugCXXFlag("/RTC1")

endif (WIN32)



project(SuPyModes LANGUAGES CXX)

set(CMAKE_CXX_STANDARD          17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

set(CMAKE_LIBRARY_OUTPUT_DIRECTORY ${PROJECT_BINARY_DIR}/../SuPyMode/bin)


add_subdirectory(extern/pybind11)
add_subdirectory(extern/eigen)

find_package( Python ${VERSION} COMPONENTS Interpreter Development EXACT )

find_package( Eigen3 REQUIRED NO_MODULE )

set(EIGEN_INCLUDE_DIR ../extern/eigen)


include_directories( ${EIGEN_INCLUDE_DIR} )
include_directories( ${PYTHON_INCLUDE_DIRS} ${PYBIND11_INCLUDE_DIRS} )

pybind11_add_module(EigenSolver MODULE ../SuPyMode/includes/interface.cpp  )
target_link_libraries(EigenSolver PRIVATE Eigen3::Eigen )
target_compile_options (EigenSolver PRIVATE -O2)



install(TARGETS EigenSolver DESTINATION ${CMAKE_LIBRARY_OUTPUT_DIRECTORY} CONFIGURATIONS Release)
install(TARGETS EigenSolver DESTINATION ${CMAKE_LIBRARY_OUTPUT_DIRECTORY} CONFIGURATIONS Debug)


#Upload Pypi package------------------------------------------------------------
add_custom_command(
  OUTPUT ${CMAKE_CURRENT_SOURCE_DIR}/UploadPypi.txt
  COMMAND python3.8 -m twine upload --password ${Password} --username ${Token} --repository pypi /Project/SuPyMode/output/*
  COMMENT "Upload on Pypi")

add_custom_target(UploadPypi DEPENDS ${CMAKE_CURRENT_SOURCE_DIR}/UploadPypi.txt)


#Cleaning command---------------------------------------------------------------
add_custom_command(
  OUTPUT ${CMAKE_CURRENT_SOURCE_DIR}/Clean.txt
  COMMAND rm -rf CMakeCache.txt cmake_install.cmake CMakeFiles
  COMMENT "Cleaning cmake output files")

add_custom_target(Clean DEPENDS ${CMAKE_CURRENT_SOURCE_DIR}/Clean.txt)
