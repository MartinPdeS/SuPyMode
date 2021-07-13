# Install script for directory: /home/martth/Desktop/git_project/SuPyModes/extern/eigen/Eigen

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
  file(INSTALL DESTINATION "${CMAKE_INSTALL_PREFIX}/include/eigen3/Eigen" TYPE FILE FILES
    "/home/martth/Desktop/git_project/SuPyModes/extern/eigen/Eigen/Cholesky"
    "/home/martth/Desktop/git_project/SuPyModes/extern/eigen/Eigen/CholmodSupport"
    "/home/martth/Desktop/git_project/SuPyModes/extern/eigen/Eigen/Core"
    "/home/martth/Desktop/git_project/SuPyModes/extern/eigen/Eigen/Dense"
    "/home/martth/Desktop/git_project/SuPyModes/extern/eigen/Eigen/Eigen"
    "/home/martth/Desktop/git_project/SuPyModes/extern/eigen/Eigen/Eigenvalues"
    "/home/martth/Desktop/git_project/SuPyModes/extern/eigen/Eigen/Geometry"
    "/home/martth/Desktop/git_project/SuPyModes/extern/eigen/Eigen/Householder"
    "/home/martth/Desktop/git_project/SuPyModes/extern/eigen/Eigen/IterativeLinearSolvers"
    "/home/martth/Desktop/git_project/SuPyModes/extern/eigen/Eigen/Jacobi"
    "/home/martth/Desktop/git_project/SuPyModes/extern/eigen/Eigen/LU"
    "/home/martth/Desktop/git_project/SuPyModes/extern/eigen/Eigen/MetisSupport"
    "/home/martth/Desktop/git_project/SuPyModes/extern/eigen/Eigen/OrderingMethods"
    "/home/martth/Desktop/git_project/SuPyModes/extern/eigen/Eigen/PaStiXSupport"
    "/home/martth/Desktop/git_project/SuPyModes/extern/eigen/Eigen/PardisoSupport"
    "/home/martth/Desktop/git_project/SuPyModes/extern/eigen/Eigen/QR"
    "/home/martth/Desktop/git_project/SuPyModes/extern/eigen/Eigen/QtAlignedMalloc"
    "/home/martth/Desktop/git_project/SuPyModes/extern/eigen/Eigen/SPQRSupport"
    "/home/martth/Desktop/git_project/SuPyModes/extern/eigen/Eigen/SVD"
    "/home/martth/Desktop/git_project/SuPyModes/extern/eigen/Eigen/Sparse"
    "/home/martth/Desktop/git_project/SuPyModes/extern/eigen/Eigen/SparseCholesky"
    "/home/martth/Desktop/git_project/SuPyModes/extern/eigen/Eigen/SparseCore"
    "/home/martth/Desktop/git_project/SuPyModes/extern/eigen/Eigen/SparseLU"
    "/home/martth/Desktop/git_project/SuPyModes/extern/eigen/Eigen/SparseQR"
    "/home/martth/Desktop/git_project/SuPyModes/extern/eigen/Eigen/StdDeque"
    "/home/martth/Desktop/git_project/SuPyModes/extern/eigen/Eigen/StdList"
    "/home/martth/Desktop/git_project/SuPyModes/extern/eigen/Eigen/StdVector"
    "/home/martth/Desktop/git_project/SuPyModes/extern/eigen/Eigen/SuperLUSupport"
    "/home/martth/Desktop/git_project/SuPyModes/extern/eigen/Eigen/UmfPackSupport"
    )
endif()

if("x${CMAKE_INSTALL_COMPONENT}x" STREQUAL "xDevelx" OR NOT CMAKE_INSTALL_COMPONENT)
  file(INSTALL DESTINATION "${CMAKE_INSTALL_PREFIX}/include/eigen3/Eigen" TYPE DIRECTORY FILES "/home/martth/Desktop/git_project/SuPyModes/extern/eigen/Eigen/src" FILES_MATCHING REGEX "/[^/]*\\.h$")
endif()

