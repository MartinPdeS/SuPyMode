# Install script for directory: /home/martth/Desktop/git_project/SuPyModes/extern/eigen/unsupported/Eigen

# Set the install prefix
if(NOT DEFINED CMAKE_INSTALL_PREFIX)
  set(CMAKE_INSTALL_PREFIX "/usr/local")
endif()
string(REGEX REPLACE "/$" "" CMAKE_INSTALL_PREFIX "${CMAKE_INSTALL_PREFIX}")

# Set the install configuration name.
if(NOT DEFINED CMAKE_INSTALL_CONFIG_NAME)
  if(BUILD_TYPE)
    string(REGEX REPLACE "^[^A-Za-z0-9_]+" ""
           CMAKE_INSTALL_CONFIG_NAME "${BUILD_TYPE}")
  else()
    set(CMAKE_INSTALL_CONFIG_NAME "Release")
  endif()
  message(STATUS "Install configuration: \"${CMAKE_INSTALL_CONFIG_NAME}\"")
endif()

# Set the component getting installed.
if(NOT CMAKE_INSTALL_COMPONENT)
  if(COMPONENT)
    message(STATUS "Install component: \"${COMPONENT}\"")
    set(CMAKE_INSTALL_COMPONENT "${COMPONENT}")
  else()
    set(CMAKE_INSTALL_COMPONENT)
  endif()
endif()

# Install shared libraries without execute permission?
if(NOT DEFINED CMAKE_INSTALL_SO_NO_EXE)
  set(CMAKE_INSTALL_SO_NO_EXE "1")
endif()

# Is this installation the result of a crosscompile?
if(NOT DEFINED CMAKE_CROSSCOMPILING)
  set(CMAKE_CROSSCOMPILING "FALSE")
endif()

if("x${CMAKE_INSTALL_COMPONENT}x" STREQUAL "xDevelx" OR NOT CMAKE_INSTALL_COMPONENT)
  file(INSTALL DESTINATION "${CMAKE_INSTALL_PREFIX}/include/eigen3/unsupported/Eigen" TYPE FILE FILES
    "/home/martth/Desktop/git_project/SuPyModes/extern/eigen/unsupported/Eigen/AdolcForward"
    "/home/martth/Desktop/git_project/SuPyModes/extern/eigen/unsupported/Eigen/AlignedVector3"
    "/home/martth/Desktop/git_project/SuPyModes/extern/eigen/unsupported/Eigen/ArpackSupport"
    "/home/martth/Desktop/git_project/SuPyModes/extern/eigen/unsupported/Eigen/AutoDiff"
    "/home/martth/Desktop/git_project/SuPyModes/extern/eigen/unsupported/Eigen/BVH"
    "/home/martth/Desktop/git_project/SuPyModes/extern/eigen/unsupported/Eigen/EulerAngles"
    "/home/martth/Desktop/git_project/SuPyModes/extern/eigen/unsupported/Eigen/FFT"
    "/home/martth/Desktop/git_project/SuPyModes/extern/eigen/unsupported/Eigen/IterativeSolvers"
    "/home/martth/Desktop/git_project/SuPyModes/extern/eigen/unsupported/Eigen/KroneckerProduct"
    "/home/martth/Desktop/git_project/SuPyModes/extern/eigen/unsupported/Eigen/LevenbergMarquardt"
    "/home/martth/Desktop/git_project/SuPyModes/extern/eigen/unsupported/Eigen/MatrixFunctions"
    "/home/martth/Desktop/git_project/SuPyModes/extern/eigen/unsupported/Eigen/MoreVectorization"
    "/home/martth/Desktop/git_project/SuPyModes/extern/eigen/unsupported/Eigen/MPRealSupport"
    "/home/martth/Desktop/git_project/SuPyModes/extern/eigen/unsupported/Eigen/NonLinearOptimization"
    "/home/martth/Desktop/git_project/SuPyModes/extern/eigen/unsupported/Eigen/NumericalDiff"
    "/home/martth/Desktop/git_project/SuPyModes/extern/eigen/unsupported/Eigen/OpenGLSupport"
    "/home/martth/Desktop/git_project/SuPyModes/extern/eigen/unsupported/Eigen/Polynomials"
    "/home/martth/Desktop/git_project/SuPyModes/extern/eigen/unsupported/Eigen/Skyline"
    "/home/martth/Desktop/git_project/SuPyModes/extern/eigen/unsupported/Eigen/SparseExtra"
    "/home/martth/Desktop/git_project/SuPyModes/extern/eigen/unsupported/Eigen/SpecialFunctions"
    "/home/martth/Desktop/git_project/SuPyModes/extern/eigen/unsupported/Eigen/Splines"
    )
endif()

if("x${CMAKE_INSTALL_COMPONENT}x" STREQUAL "xDevelx" OR NOT CMAKE_INSTALL_COMPONENT)
  file(INSTALL DESTINATION "${CMAKE_INSTALL_PREFIX}/include/eigen3/unsupported/Eigen" TYPE DIRECTORY FILES "/home/martth/Desktop/git_project/SuPyModes/extern/eigen/unsupported/Eigen/src" FILES_MATCHING REGEX "/[^/]*\\.h$")
endif()

if(NOT CMAKE_INSTALL_LOCAL_ONLY)
  # Include the install script for each subdirectory.
  include("/home/martth/Desktop/git_project/SuPyModes/build/extern/eigen/unsupported/Eigen/CXX11/cmake_install.cmake")

endif()

