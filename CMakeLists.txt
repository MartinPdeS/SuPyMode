cmake_minimum_required(VERSION 3.0)


#set(CMAKE_SYSTEM_PROCESSOR x86_64)
#set(CMAKE_CROSSCOMPILING 1)
#set(CMAKE_CXX_COMPILER_WORKS 1)
#set(TOOLCHAIN_PREFIX /usr/bin/x86_64-w64-mingw32)
#set(CMAKE_CXX_COMPILER ${TOOLCHAIN_PREFIX}-g++)

macro(RemoveDebugCXXFlag flag)
    string(REPLACE "${flag}" "" CMAKE_CXX_FLAGS_DEBUG "${CMAKE_CXX_FLAGS_DEBUG}")
endmacro()

macro(RemoveReleaseCXXFlag flag)
    string(REPLACE "${flag}" "" CMAKE_CXX_FLAGS_RELEASE "${CMAKE_CXX_FLAGS_RELEASE}")
endmacro()

project(SuPyModes LANGUAGES CXX)

set(CMAKE_CXX_STANDARD          17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(Token    "$ENV{PipUsername}")
set(Password "$ENV{PipPassword}")
set(Boost_USE_STATIC_LIBS ON)
set(CMAKE_LIBRARY_OUTPUT_DIRECTORY ${PROJECT_BINARY_DIR}/../SuPyMode/bin)


add_subdirectory(extern/pybind11)

find_package(Python 3 COMPONENTS Interpreter Development EXACT)
find_package( Boost REQUIRED )
find_package(Eigen3 REQUIRED NO_MODULE)
find_package(PythonInterp 3 REQUIRED)

set(EIGEN_INCLUDE_DIR ../extern/eigen)

include_directories(${EIGEN_INCLUDE_DIR})
include_directories(${Boost_INCLUDE_DIR})
include_directories( ${PYTHON_INCLUDE_DIRS} ${PYBIND11_INCLUDE_DIRS} )

pybind11_add_module(EigenSolver MODULE ../SuPyMode/includes/interface.cpp  )
target_link_libraries(EigenSolver PRIVATE Eigen3::Eigen ${Boost_LIBRARIES})
target_compile_options (EigenSolver PRIVATE -O2)

RemoveReleaseCXXFlag("/RTC1")
RemoveDebugCXXFlag("/RTC1")

install(TARGETS EigenSolver DESTINATION ${CMAKE_LIBRARY_OUTPUT_DIRECTORY} CONFIGURATIONS Release)
install(TARGETS EigenSolver DESTINATION ${CMAKE_LIBRARY_OUTPUT_DIRECTORY} CONFIGURATIONS Debug)


#Upload Pypi package------------------------------------------------------------
add_custom_command(
  OUTPUT ${CMAKE_CURRENT_SOURCE_DIR}/UploadPypi.txt
  COMMAND python3.8 -m twine upload --password ${Password} --username ${Token} --repository pypi ./dist/*
  COMMENT "Upload on Pypi")

add_custom_target(UploadPypi DEPENDS ${CMAKE_CURRENT_SOURCE_DIR}/UploadPypi.txt)


#Cleaning command---------------------------------------------------------------
add_custom_command(
  OUTPUT ${CMAKE_CURRENT_SOURCE_DIR}/Clean.txt
  COMMAND rm -rf CMakeCache.txt cmake_install.cmake CMakeFiles
  COMMENT "Cleaning cmake output files")

add_custom_target(Clean DEPENDS ${CMAKE_CURRENT_SOURCE_DIR}/Clean.txt)
