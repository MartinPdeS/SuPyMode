# CMAKE generated file: DO NOT EDIT!
# Generated by "Unix Makefiles" Generator, CMake Version 3.16

# Delete rule output on recipe failure.
.DELETE_ON_ERROR:


#=============================================================================
# Special targets provided by cmake.

# Disable implicit rules so canonical targets will work.
.SUFFIXES:


# Remove some rules from gmake that .SUFFIXES does not remove.
SUFFIXES =

.SUFFIXES: .hpux_make_needs_suffix_list


# Suppress display of executed commands.
$(VERBOSE).SILENT:


# A target that is always out of date.
cmake_force:

.PHONY : cmake_force

#=============================================================================
# Set environment variables for the build.

# The shell in which to execute make rules.
SHELL = /bin/sh

# The CMake executable.
CMAKE_COMMAND = /usr/bin/cmake

# The command to remove a file.
RM = /usr/bin/cmake -E remove -f

# Escaping for special characters.
EQUALS = =

# The top-level source directory on which CMake was run.
CMAKE_SOURCE_DIR = /home/martth/Desktop/git_project/SuPyModes

# The top-level build directory on which CMake was run.
CMAKE_BINARY_DIR = /home/martth/Desktop/git_project/SuPyModes/temp

# Utility rule file for sparse_ref.

# Include the progress variables for this target.
include extern/eigen/test/CMakeFiles/sparse_ref.dir/progress.make

sparse_ref: extern/eigen/test/CMakeFiles/sparse_ref.dir/build.make

.PHONY : sparse_ref

# Rule to build all files generated by this target.
extern/eigen/test/CMakeFiles/sparse_ref.dir/build: sparse_ref

.PHONY : extern/eigen/test/CMakeFiles/sparse_ref.dir/build

extern/eigen/test/CMakeFiles/sparse_ref.dir/clean:
	cd /home/martth/Desktop/git_project/SuPyModes/temp/extern/eigen/test && $(CMAKE_COMMAND) -P CMakeFiles/sparse_ref.dir/cmake_clean.cmake
.PHONY : extern/eigen/test/CMakeFiles/sparse_ref.dir/clean

extern/eigen/test/CMakeFiles/sparse_ref.dir/depend:
	cd /home/martth/Desktop/git_project/SuPyModes/temp && $(CMAKE_COMMAND) -E cmake_depends "Unix Makefiles" /home/martth/Desktop/git_project/SuPyModes /home/martth/Desktop/git_project/SuPyModes/extern/eigen/test /home/martth/Desktop/git_project/SuPyModes/temp /home/martth/Desktop/git_project/SuPyModes/temp/extern/eigen/test /home/martth/Desktop/git_project/SuPyModes/temp/extern/eigen/test/CMakeFiles/sparse_ref.dir/DependInfo.cmake --color=$(COLOR)
.PHONY : extern/eigen/test/CMakeFiles/sparse_ref.dir/depend

