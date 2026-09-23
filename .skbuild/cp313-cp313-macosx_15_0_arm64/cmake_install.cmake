# Install script for directory: /Users/m.poinsinetdesivry-houle/Desktop/GitRepositories/SuPyMode

# Set the install prefix
if(NOT DEFINED CMAKE_INSTALL_PREFIX)
  set(CMAKE_INSTALL_PREFIX "/var/folders/nx/6l9ysx0j7yb5szphhl8jskyh0000gn/T/tmpk4luksx8/wheel/platlib")
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

# Is this installation the result of a crosscompile?
if(NOT DEFINED CMAKE_CROSSCOMPILING)
  set(CMAKE_CROSSCOMPILING "FALSE")
endif()

# Set path to fallback-tool for dependency-resolution.
if(NOT DEFINED CMAKE_OBJDUMP)
  set(CMAKE_OBJDUMP "/Library/Developer/CommandLineTools/usr/bin/objdump")
endif()

if(NOT CMAKE_INSTALL_LOCAL_ONLY)
  # Include the install script for the subdirectory.
  include("/Users/m.poinsinetdesivry-houle/Desktop/GitRepositories/SuPyMode/.skbuild/cp313-cp313-macosx_15_0_arm64/SuPyMode/cpp/mesh/cmake_install.cmake")
endif()

if(NOT CMAKE_INSTALL_LOCAL_ONLY)
  # Include the install script for the subdirectory.
  include("/Users/m.poinsinetdesivry-houle/Desktop/GitRepositories/SuPyMode/.skbuild/cp313-cp313-macosx_15_0_arm64/SuPyMode/cpp/model_parameters/cmake_install.cmake")
endif()

if(NOT CMAKE_INSTALL_LOCAL_ONLY)
  # Include the install script for the subdirectory.
  include("/Users/m.poinsinetdesivry-houle/Desktop/GitRepositories/SuPyMode/.skbuild/cp313-cp313-macosx_15_0_arm64/SuPyMode/cpp/supermode/cmake_install.cmake")
endif()

if(NOT CMAKE_INSTALL_LOCAL_ONLY)
  # Include the install script for the subdirectory.
  include("/Users/m.poinsinetdesivry-houle/Desktop/GitRepositories/SuPyMode/.skbuild/cp313-cp313-macosx_15_0_arm64/SuPyMode/cpp/eigensolver/cmake_install.cmake")
endif()

if(CMAKE_INSTALL_COMPONENT STREQUAL "Unspecified" OR NOT CMAKE_INSTALL_COMPONENT)
  list(APPEND CMAKE_ABSOLUTE_DESTINATION_FILES
   "/Users/m.poinsinetdesivry-houle/Desktop/GitRepositories/SuPyMode/SuPyMode/binary/libmesh.a")
  if(CMAKE_WARN_ON_ABSOLUTE_INSTALL_DESTINATION)
    message(WARNING "ABSOLUTE path INSTALL DESTINATION : ${CMAKE_ABSOLUTE_DESTINATION_FILES}")
  endif()
  if(CMAKE_ERROR_ON_ABSOLUTE_INSTALL_DESTINATION)
    message(FATAL_ERROR "ABSOLUTE path INSTALL DESTINATION forbidden (by caller): ${CMAKE_ABSOLUTE_DESTINATION_FILES}")
  endif()
  file(INSTALL DESTINATION "/Users/m.poinsinetdesivry-houle/Desktop/GitRepositories/SuPyMode/SuPyMode/binary" TYPE STATIC_LIBRARY FILES "/Users/m.poinsinetdesivry-houle/Desktop/GitRepositories/SuPyMode/.skbuild/cp313-cp313-macosx_15_0_arm64/SuPyMode/cpp/mesh/libmesh.a")
  if(EXISTS "$ENV{DESTDIR}/Users/m.poinsinetdesivry-houle/Desktop/GitRepositories/SuPyMode/SuPyMode/binary/libmesh.a" AND
     NOT IS_SYMLINK "$ENV{DESTDIR}/Users/m.poinsinetdesivry-houle/Desktop/GitRepositories/SuPyMode/SuPyMode/binary/libmesh.a")
    execute_process(COMMAND "/Library/Developer/CommandLineTools/usr/bin/ranlib" "$ENV{DESTDIR}/Users/m.poinsinetdesivry-houle/Desktop/GitRepositories/SuPyMode/SuPyMode/binary/libmesh.a")
  endif()
endif()

if(CMAKE_INSTALL_COMPONENT STREQUAL "Unspecified" OR NOT CMAKE_INSTALL_COMPONENT)
  list(APPEND CMAKE_ABSOLUTE_DESTINATION_FILES
   "/Users/m.poinsinetdesivry-houle/Desktop/GitRepositories/SuPyMode/SuPyMode/binary/interface_mesh.cpython-313-darwin.so")
  if(CMAKE_WARN_ON_ABSOLUTE_INSTALL_DESTINATION)
    message(WARNING "ABSOLUTE path INSTALL DESTINATION : ${CMAKE_ABSOLUTE_DESTINATION_FILES}")
  endif()
  if(CMAKE_ERROR_ON_ABSOLUTE_INSTALL_DESTINATION)
    message(FATAL_ERROR "ABSOLUTE path INSTALL DESTINATION forbidden (by caller): ${CMAKE_ABSOLUTE_DESTINATION_FILES}")
  endif()
  file(INSTALL DESTINATION "/Users/m.poinsinetdesivry-houle/Desktop/GitRepositories/SuPyMode/SuPyMode/binary" TYPE MODULE FILES "/Users/m.poinsinetdesivry-houle/Desktop/GitRepositories/SuPyMode/.skbuild/cp313-cp313-macosx_15_0_arm64/SuPyMode/cpp/mesh/interface_mesh.cpython-313-darwin.so")
  if(EXISTS "$ENV{DESTDIR}/Users/m.poinsinetdesivry-houle/Desktop/GitRepositories/SuPyMode/SuPyMode/binary/interface_mesh.cpython-313-darwin.so" AND
     NOT IS_SYMLINK "$ENV{DESTDIR}/Users/m.poinsinetdesivry-houle/Desktop/GitRepositories/SuPyMode/SuPyMode/binary/interface_mesh.cpython-313-darwin.so")
    if(CMAKE_INSTALL_DO_STRIP)
      execute_process(COMMAND "/Library/Developer/CommandLineTools/usr/bin/strip" -x "$ENV{DESTDIR}/Users/m.poinsinetdesivry-houle/Desktop/GitRepositories/SuPyMode/SuPyMode/binary/interface_mesh.cpython-313-darwin.so")
    endif()
  endif()
endif()

if(CMAKE_INSTALL_COMPONENT STREQUAL "Unspecified" OR NOT CMAKE_INSTALL_COMPONENT)
  list(APPEND CMAKE_ABSOLUTE_DESTINATION_FILES
   "/Users/m.poinsinetdesivry-houle/Desktop/GitRepositories/SuPyMode/SuPyMode/binary/libmodel_parameters.a")
  if(CMAKE_WARN_ON_ABSOLUTE_INSTALL_DESTINATION)
    message(WARNING "ABSOLUTE path INSTALL DESTINATION : ${CMAKE_ABSOLUTE_DESTINATION_FILES}")
  endif()
  if(CMAKE_ERROR_ON_ABSOLUTE_INSTALL_DESTINATION)
    message(FATAL_ERROR "ABSOLUTE path INSTALL DESTINATION forbidden (by caller): ${CMAKE_ABSOLUTE_DESTINATION_FILES}")
  endif()
  file(INSTALL DESTINATION "/Users/m.poinsinetdesivry-houle/Desktop/GitRepositories/SuPyMode/SuPyMode/binary" TYPE STATIC_LIBRARY FILES "/Users/m.poinsinetdesivry-houle/Desktop/GitRepositories/SuPyMode/.skbuild/cp313-cp313-macosx_15_0_arm64/SuPyMode/cpp/model_parameters/libmodel_parameters.a")
  if(EXISTS "$ENV{DESTDIR}/Users/m.poinsinetdesivry-houle/Desktop/GitRepositories/SuPyMode/SuPyMode/binary/libmodel_parameters.a" AND
     NOT IS_SYMLINK "$ENV{DESTDIR}/Users/m.poinsinetdesivry-houle/Desktop/GitRepositories/SuPyMode/SuPyMode/binary/libmodel_parameters.a")
    execute_process(COMMAND "/Library/Developer/CommandLineTools/usr/bin/ranlib" "$ENV{DESTDIR}/Users/m.poinsinetdesivry-houle/Desktop/GitRepositories/SuPyMode/SuPyMode/binary/libmodel_parameters.a")
  endif()
endif()

if(CMAKE_INSTALL_COMPONENT STREQUAL "Unspecified" OR NOT CMAKE_INSTALL_COMPONENT)
  list(APPEND CMAKE_ABSOLUTE_DESTINATION_FILES
   "/Users/m.poinsinetdesivry-houle/Desktop/GitRepositories/SuPyMode/SuPyMode/binary/interface_model_parameters.cpython-313-darwin.so")
  if(CMAKE_WARN_ON_ABSOLUTE_INSTALL_DESTINATION)
    message(WARNING "ABSOLUTE path INSTALL DESTINATION : ${CMAKE_ABSOLUTE_DESTINATION_FILES}")
  endif()
  if(CMAKE_ERROR_ON_ABSOLUTE_INSTALL_DESTINATION)
    message(FATAL_ERROR "ABSOLUTE path INSTALL DESTINATION forbidden (by caller): ${CMAKE_ABSOLUTE_DESTINATION_FILES}")
  endif()
  file(INSTALL DESTINATION "/Users/m.poinsinetdesivry-houle/Desktop/GitRepositories/SuPyMode/SuPyMode/binary" TYPE MODULE FILES "/Users/m.poinsinetdesivry-houle/Desktop/GitRepositories/SuPyMode/.skbuild/cp313-cp313-macosx_15_0_arm64/SuPyMode/cpp/model_parameters/interface_model_parameters.cpython-313-darwin.so")
  if(EXISTS "$ENV{DESTDIR}/Users/m.poinsinetdesivry-houle/Desktop/GitRepositories/SuPyMode/SuPyMode/binary/interface_model_parameters.cpython-313-darwin.so" AND
     NOT IS_SYMLINK "$ENV{DESTDIR}/Users/m.poinsinetdesivry-houle/Desktop/GitRepositories/SuPyMode/SuPyMode/binary/interface_model_parameters.cpython-313-darwin.so")
    if(CMAKE_INSTALL_DO_STRIP)
      execute_process(COMMAND "/Library/Developer/CommandLineTools/usr/bin/strip" -x "$ENV{DESTDIR}/Users/m.poinsinetdesivry-houle/Desktop/GitRepositories/SuPyMode/SuPyMode/binary/interface_model_parameters.cpython-313-darwin.so")
    endif()
  endif()
endif()

if(CMAKE_INSTALL_COMPONENT STREQUAL "Unspecified" OR NOT CMAKE_INSTALL_COMPONENT)
  list(APPEND CMAKE_ABSOLUTE_DESTINATION_FILES
   "/Users/m.poinsinetdesivry-houle/Desktop/GitRepositories/SuPyMode/SuPyMode/binary/libsupermode.a")
  if(CMAKE_WARN_ON_ABSOLUTE_INSTALL_DESTINATION)
    message(WARNING "ABSOLUTE path INSTALL DESTINATION : ${CMAKE_ABSOLUTE_DESTINATION_FILES}")
  endif()
  if(CMAKE_ERROR_ON_ABSOLUTE_INSTALL_DESTINATION)
    message(FATAL_ERROR "ABSOLUTE path INSTALL DESTINATION forbidden (by caller): ${CMAKE_ABSOLUTE_DESTINATION_FILES}")
  endif()
  file(INSTALL DESTINATION "/Users/m.poinsinetdesivry-houle/Desktop/GitRepositories/SuPyMode/SuPyMode/binary" TYPE STATIC_LIBRARY FILES "/Users/m.poinsinetdesivry-houle/Desktop/GitRepositories/SuPyMode/.skbuild/cp313-cp313-macosx_15_0_arm64/SuPyMode/cpp/supermode/libsupermode.a")
  if(EXISTS "$ENV{DESTDIR}/Users/m.poinsinetdesivry-houle/Desktop/GitRepositories/SuPyMode/SuPyMode/binary/libsupermode.a" AND
     NOT IS_SYMLINK "$ENV{DESTDIR}/Users/m.poinsinetdesivry-houle/Desktop/GitRepositories/SuPyMode/SuPyMode/binary/libsupermode.a")
    execute_process(COMMAND "/Library/Developer/CommandLineTools/usr/bin/ranlib" "$ENV{DESTDIR}/Users/m.poinsinetdesivry-houle/Desktop/GitRepositories/SuPyMode/SuPyMode/binary/libsupermode.a")
  endif()
endif()

if(CMAKE_INSTALL_COMPONENT STREQUAL "Unspecified" OR NOT CMAKE_INSTALL_COMPONENT)
  list(APPEND CMAKE_ABSOLUTE_DESTINATION_FILES
   "/Users/m.poinsinetdesivry-houle/Desktop/GitRepositories/SuPyMode/SuPyMode/binary/interface_supermode.cpython-313-darwin.so")
  if(CMAKE_WARN_ON_ABSOLUTE_INSTALL_DESTINATION)
    message(WARNING "ABSOLUTE path INSTALL DESTINATION : ${CMAKE_ABSOLUTE_DESTINATION_FILES}")
  endif()
  if(CMAKE_ERROR_ON_ABSOLUTE_INSTALL_DESTINATION)
    message(FATAL_ERROR "ABSOLUTE path INSTALL DESTINATION forbidden (by caller): ${CMAKE_ABSOLUTE_DESTINATION_FILES}")
  endif()
  file(INSTALL DESTINATION "/Users/m.poinsinetdesivry-houle/Desktop/GitRepositories/SuPyMode/SuPyMode/binary" TYPE MODULE FILES "/Users/m.poinsinetdesivry-houle/Desktop/GitRepositories/SuPyMode/.skbuild/cp313-cp313-macosx_15_0_arm64/SuPyMode/cpp/supermode/interface_supermode.cpython-313-darwin.so")
  if(EXISTS "$ENV{DESTDIR}/Users/m.poinsinetdesivry-houle/Desktop/GitRepositories/SuPyMode/SuPyMode/binary/interface_supermode.cpython-313-darwin.so" AND
     NOT IS_SYMLINK "$ENV{DESTDIR}/Users/m.poinsinetdesivry-houle/Desktop/GitRepositories/SuPyMode/SuPyMode/binary/interface_supermode.cpython-313-darwin.so")
    if(CMAKE_INSTALL_DO_STRIP)
      execute_process(COMMAND "/Library/Developer/CommandLineTools/usr/bin/strip" -x "$ENV{DESTDIR}/Users/m.poinsinetdesivry-houle/Desktop/GitRepositories/SuPyMode/SuPyMode/binary/interface_supermode.cpython-313-darwin.so")
    endif()
  endif()
endif()

if(CMAKE_INSTALL_COMPONENT STREQUAL "Unspecified" OR NOT CMAKE_INSTALL_COMPONENT)
  list(APPEND CMAKE_ABSOLUTE_DESTINATION_FILES
   "/Users/m.poinsinetdesivry-houle/Desktop/GitRepositories/SuPyMode/SuPyMode/binary/libeigensolver.a")
  if(CMAKE_WARN_ON_ABSOLUTE_INSTALL_DESTINATION)
    message(WARNING "ABSOLUTE path INSTALL DESTINATION : ${CMAKE_ABSOLUTE_DESTINATION_FILES}")
  endif()
  if(CMAKE_ERROR_ON_ABSOLUTE_INSTALL_DESTINATION)
    message(FATAL_ERROR "ABSOLUTE path INSTALL DESTINATION forbidden (by caller): ${CMAKE_ABSOLUTE_DESTINATION_FILES}")
  endif()
  file(INSTALL DESTINATION "/Users/m.poinsinetdesivry-houle/Desktop/GitRepositories/SuPyMode/SuPyMode/binary" TYPE STATIC_LIBRARY FILES "/Users/m.poinsinetdesivry-houle/Desktop/GitRepositories/SuPyMode/.skbuild/cp313-cp313-macosx_15_0_arm64/SuPyMode/cpp/eigensolver/libeigensolver.a")
  if(EXISTS "$ENV{DESTDIR}/Users/m.poinsinetdesivry-houle/Desktop/GitRepositories/SuPyMode/SuPyMode/binary/libeigensolver.a" AND
     NOT IS_SYMLINK "$ENV{DESTDIR}/Users/m.poinsinetdesivry-houle/Desktop/GitRepositories/SuPyMode/SuPyMode/binary/libeigensolver.a")
    execute_process(COMMAND "/Library/Developer/CommandLineTools/usr/bin/ranlib" "$ENV{DESTDIR}/Users/m.poinsinetdesivry-houle/Desktop/GitRepositories/SuPyMode/SuPyMode/binary/libeigensolver.a")
  endif()
endif()

if(CMAKE_INSTALL_COMPONENT STREQUAL "Unspecified" OR NOT CMAKE_INSTALL_COMPONENT)
  list(APPEND CMAKE_ABSOLUTE_DESTINATION_FILES
   "/Users/m.poinsinetdesivry-houle/Desktop/GitRepositories/SuPyMode/SuPyMode/binary/interface_eigensolver.cpython-313-darwin.so")
  if(CMAKE_WARN_ON_ABSOLUTE_INSTALL_DESTINATION)
    message(WARNING "ABSOLUTE path INSTALL DESTINATION : ${CMAKE_ABSOLUTE_DESTINATION_FILES}")
  endif()
  if(CMAKE_ERROR_ON_ABSOLUTE_INSTALL_DESTINATION)
    message(FATAL_ERROR "ABSOLUTE path INSTALL DESTINATION forbidden (by caller): ${CMAKE_ABSOLUTE_DESTINATION_FILES}")
  endif()
  file(INSTALL DESTINATION "/Users/m.poinsinetdesivry-houle/Desktop/GitRepositories/SuPyMode/SuPyMode/binary" TYPE MODULE FILES "/Users/m.poinsinetdesivry-houle/Desktop/GitRepositories/SuPyMode/.skbuild/cp313-cp313-macosx_15_0_arm64/SuPyMode/cpp/eigensolver/interface_eigensolver.cpython-313-darwin.so")
  if(EXISTS "$ENV{DESTDIR}/Users/m.poinsinetdesivry-houle/Desktop/GitRepositories/SuPyMode/SuPyMode/binary/interface_eigensolver.cpython-313-darwin.so" AND
     NOT IS_SYMLINK "$ENV{DESTDIR}/Users/m.poinsinetdesivry-houle/Desktop/GitRepositories/SuPyMode/SuPyMode/binary/interface_eigensolver.cpython-313-darwin.so")
    if(CMAKE_INSTALL_DO_STRIP)
      execute_process(COMMAND "/Library/Developer/CommandLineTools/usr/bin/strip" -x "$ENV{DESTDIR}/Users/m.poinsinetdesivry-houle/Desktop/GitRepositories/SuPyMode/SuPyMode/binary/interface_eigensolver.cpython-313-darwin.so")
    endif()
  endif()
endif()

string(REPLACE ";" "\n" CMAKE_INSTALL_MANIFEST_CONTENT
       "${CMAKE_INSTALL_MANIFEST_FILES}")
if(CMAKE_INSTALL_LOCAL_ONLY)
  file(WRITE "/Users/m.poinsinetdesivry-houle/Desktop/GitRepositories/SuPyMode/.skbuild/cp313-cp313-macosx_15_0_arm64/install_local_manifest.txt"
     "${CMAKE_INSTALL_MANIFEST_CONTENT}")
endif()
if(CMAKE_INSTALL_COMPONENT)
  if(CMAKE_INSTALL_COMPONENT MATCHES "^[a-zA-Z0-9_.+-]+$")
    set(CMAKE_INSTALL_MANIFEST "install_manifest_${CMAKE_INSTALL_COMPONENT}.txt")
  else()
    string(MD5 CMAKE_INST_COMP_HASH "${CMAKE_INSTALL_COMPONENT}")
    set(CMAKE_INSTALL_MANIFEST "install_manifest_${CMAKE_INST_COMP_HASH}.txt")
    unset(CMAKE_INST_COMP_HASH)
  endif()
else()
  set(CMAKE_INSTALL_MANIFEST "install_manifest.txt")
endif()

if(NOT CMAKE_INSTALL_LOCAL_ONLY)
  file(WRITE "/Users/m.poinsinetdesivry-houle/Desktop/GitRepositories/SuPyMode/.skbuild/cp313-cp313-macosx_15_0_arm64/${CMAKE_INSTALL_MANIFEST}"
     "${CMAKE_INSTALL_MANIFEST_CONTENT}")
endif()
