# CMAKE generated file: DO NOT EDIT!
# Generated by "Unix Makefiles" Generator, CMake Version 2.4

#=============================================================================
# Special targets provided by cmake.

# Disable implicit rules so canoncical targets will work.
.SUFFIXES:

.SUFFIXES: .hpux_make_needs_suffix_list

# Suppress display of executed commands.
$(VERBOSE).SILENT:

# A target that is always out of date.
cmake_force:

#=============================================================================
# Set environment variables for the build.

# The shell in which to execute make rules.
SHELL = /bin/sh

# The CMake executable.
CMAKE_COMMAND = /usr/bin/cmake

# The command to remove a file.
RM = /usr/bin/cmake -E remove -f

# The program to use to edit the cache.
CMAKE_EDIT_COMMAND = /usr/bin/ccmake

# The top-level source directory on which CMake was run.
CMAKE_SOURCE_DIR = /Users/bhattgopal/dev/cpp/qgis_1.0.2

# The top-level build directory on which CMake was run.
CMAKE_BINARY_DIR = /Users/bhattgopal/dev/cpp/qgis_1.0.2/build

# Utility rule file for Info.plist.

src/mac/Contents/CMakeFiles/Info.plist: ../.svn/entries
	cd /Users/bhattgopal/dev/cpp/qgis_1.0.2/build/src/mac/Contents && /usr/bin/cmake -D CURRENT_SOURCE_DIR=/Users/bhattgopal/dev/cpp/qgis_1.0.2/src/mac/Contents -D SOURCE_DIR=/Users/bhattgopal/dev/cpp/qgis_1.0.2 -D VERSION=1.0.2-Kore -P /Users/bhattgopal/dev/cpp/qgis_1.0.2/src/mac/Contents/Info.plist.cmake

Info.plist: src/mac/Contents/CMakeFiles/Info.plist.dir/build.make

# Rule to build all files generated by this target.
src/mac/Contents/CMakeFiles/Info.plist.dir/build: Info.plist
src/mac/Contents/CMakeFiles/Info.plist.dir/build: src/mac/Contents/CMakeFiles/Info.plist

src/mac/Contents/CMakeFiles/Info.plist.dir/clean:
	cd /Users/bhattgopal/dev/cpp/qgis_1.0.2/build/src/mac/Contents && $(CMAKE_COMMAND) -P CMakeFiles/Info.plist.dir/cmake_clean.cmake
