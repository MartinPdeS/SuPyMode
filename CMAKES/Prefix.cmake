

include(${CMAKE_CURRENT_SOURCE_DIR}/CMAKES/Utils.cmake)

IF (WIN32)
    message("_______________________Setup compilation for Windows OS_______________________")
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

ELSEIF(LINUX)
    message("_______________________Setup compilation for Linux OS  _______________________")

ELSEIF(APPLE)
    message("_______________________Setup compilation for Mac OS    _______________________")


ENDIF()



project(PyMieSim LANGUAGES Fortran CXX)


set(Token    "$ENV{PyPiUsername}")
set(Password "$ENV{PyPiPassword}")
