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

# Include any dependencies generated for this target.
include src/providers/wms/CMakeFiles/wmsprovider.dir/depend.make

# Include the progress variables for this target.
include src/providers/wms/CMakeFiles/wmsprovider.dir/progress.make

# Include the compile flags for this target's objects.
include src/providers/wms/CMakeFiles/wmsprovider.dir/flags.make

src/providers/wms/CMakeFiles/wmsprovider.dir/depend.make.mark: src/providers/wms/CMakeFiles/wmsprovider.dir/flags.make
src/providers/wms/CMakeFiles/wmsprovider.dir/depend.make.mark: ../src/providers/wms/qgswmsprovider.cpp

src/providers/wms/CMakeFiles/wmsprovider.dir/qgswmsprovider.o: src/providers/wms/CMakeFiles/wmsprovider.dir/flags.make
src/providers/wms/CMakeFiles/wmsprovider.dir/qgswmsprovider.o: ../src/providers/wms/qgswmsprovider.cpp
	$(CMAKE_COMMAND) -E cmake_progress_report /Users/bhattgopal/dev/cpp/qgis_1.0.2/build/CMakeFiles $(CMAKE_PROGRESS_1)
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Building CXX object src/providers/wms/CMakeFiles/wmsprovider.dir/qgswmsprovider.o"
	/usr/bin/c++   $(CXX_FLAGS) -o src/providers/wms/CMakeFiles/wmsprovider.dir/qgswmsprovider.o -c /Users/bhattgopal/dev/cpp/qgis_1.0.2/src/providers/wms/qgswmsprovider.cpp

src/providers/wms/CMakeFiles/wmsprovider.dir/qgswmsprovider.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to src/providers/wms/CMakeFiles/wmsprovider.dir/qgswmsprovider.i"
	/usr/bin/c++  $(CXX_FLAGS) -E /Users/bhattgopal/dev/cpp/qgis_1.0.2/src/providers/wms/qgswmsprovider.cpp > src/providers/wms/CMakeFiles/wmsprovider.dir/qgswmsprovider.i

src/providers/wms/CMakeFiles/wmsprovider.dir/qgswmsprovider.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly src/providers/wms/CMakeFiles/wmsprovider.dir/qgswmsprovider.s"
	/usr/bin/c++  $(CXX_FLAGS) -S /Users/bhattgopal/dev/cpp/qgis_1.0.2/src/providers/wms/qgswmsprovider.cpp -o src/providers/wms/CMakeFiles/wmsprovider.dir/qgswmsprovider.s

src/providers/wms/CMakeFiles/wmsprovider.dir/qgswmsprovider.o.requires:

src/providers/wms/CMakeFiles/wmsprovider.dir/qgswmsprovider.o.provides: src/providers/wms/CMakeFiles/wmsprovider.dir/qgswmsprovider.o.requires
	$(MAKE) -f src/providers/wms/CMakeFiles/wmsprovider.dir/build.make src/providers/wms/CMakeFiles/wmsprovider.dir/qgswmsprovider.o.provides.build

src/providers/wms/CMakeFiles/wmsprovider.dir/qgswmsprovider.o.provides.build: src/providers/wms/CMakeFiles/wmsprovider.dir/qgswmsprovider.o

src/providers/wms/CMakeFiles/wmsprovider.dir/depend.make.mark: src/providers/wms/CMakeFiles/wmsprovider.dir/flags.make
src/providers/wms/CMakeFiles/wmsprovider.dir/depend.make.mark: src/providers/wms/moc_qgswmsprovider.cxx

src/providers/wms/CMakeFiles/wmsprovider.dir/moc_qgswmsprovider.o: src/providers/wms/CMakeFiles/wmsprovider.dir/flags.make
src/providers/wms/CMakeFiles/wmsprovider.dir/moc_qgswmsprovider.o: src/providers/wms/moc_qgswmsprovider.cxx
	$(CMAKE_COMMAND) -E cmake_progress_report /Users/bhattgopal/dev/cpp/qgis_1.0.2/build/CMakeFiles $(CMAKE_PROGRESS_2)
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Building CXX object src/providers/wms/CMakeFiles/wmsprovider.dir/moc_qgswmsprovider.o"
	/usr/bin/c++   $(CXX_FLAGS) -o src/providers/wms/CMakeFiles/wmsprovider.dir/moc_qgswmsprovider.o -c /Users/bhattgopal/dev/cpp/qgis_1.0.2/build/src/providers/wms/moc_qgswmsprovider.cxx

src/providers/wms/CMakeFiles/wmsprovider.dir/moc_qgswmsprovider.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to src/providers/wms/CMakeFiles/wmsprovider.dir/moc_qgswmsprovider.i"
	/usr/bin/c++  $(CXX_FLAGS) -E /Users/bhattgopal/dev/cpp/qgis_1.0.2/build/src/providers/wms/moc_qgswmsprovider.cxx > src/providers/wms/CMakeFiles/wmsprovider.dir/moc_qgswmsprovider.i

src/providers/wms/CMakeFiles/wmsprovider.dir/moc_qgswmsprovider.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly src/providers/wms/CMakeFiles/wmsprovider.dir/moc_qgswmsprovider.s"
	/usr/bin/c++  $(CXX_FLAGS) -S /Users/bhattgopal/dev/cpp/qgis_1.0.2/build/src/providers/wms/moc_qgswmsprovider.cxx -o src/providers/wms/CMakeFiles/wmsprovider.dir/moc_qgswmsprovider.s

src/providers/wms/CMakeFiles/wmsprovider.dir/moc_qgswmsprovider.o.requires:

src/providers/wms/CMakeFiles/wmsprovider.dir/moc_qgswmsprovider.o.provides: src/providers/wms/CMakeFiles/wmsprovider.dir/moc_qgswmsprovider.o.requires
	$(MAKE) -f src/providers/wms/CMakeFiles/wmsprovider.dir/build.make src/providers/wms/CMakeFiles/wmsprovider.dir/moc_qgswmsprovider.o.provides.build

src/providers/wms/CMakeFiles/wmsprovider.dir/moc_qgswmsprovider.o.provides.build: src/providers/wms/CMakeFiles/wmsprovider.dir/moc_qgswmsprovider.o

src/providers/wms/moc_qgswmsprovider.cxx: ../src/providers/wms/qgswmsprovider.h
	$(CMAKE_COMMAND) -E cmake_progress_report /Users/bhattgopal/dev/cpp/qgis_1.0.2/build/CMakeFiles $(CMAKE_PROGRESS_3)
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --blue --bold "Generating moc_qgswmsprovider.cxx"
	cd /Users/bhattgopal/dev/cpp/qgis_1.0.2/build/src/providers/wms && /Developer/Tools/Qt/moc -I /usr/include -I /Library/Frameworks/Qt3Support.framework/Headers -I /Library/Frameworks/QtSvg.framework/Headers -I /Library/Frameworks/QtGui.framework/Headers -I /Library/Frameworks/QtXml.framework/Headers -I /Library/Frameworks/QtSql.framework/Headers -I /Library/Frameworks/QtNetwork.framework/Headers -I /Library/Frameworks/QtCore.framework/Headers -I /Users/bhattgopal/dev/cpp/qgis_1.0.2/build -o /Users/bhattgopal/dev/cpp/qgis_1.0.2/build/src/providers/wms/moc_qgswmsprovider.cxx /Users/bhattgopal/dev/cpp/qgis_1.0.2/src/providers/wms/qgswmsprovider.h

src/providers/wms/CMakeFiles/wmsprovider.dir/depend: src/providers/wms/CMakeFiles/wmsprovider.dir/depend.make.mark

src/providers/wms/CMakeFiles/wmsprovider.dir/depend.make.mark: src/providers/wms/moc_qgswmsprovider.cxx
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --magenta --bold "Scanning dependencies of target wmsprovider"
	cd /Users/bhattgopal/dev/cpp/qgis_1.0.2/build && $(CMAKE_COMMAND) -E cmake_depends "Unix Makefiles" /Users/bhattgopal/dev/cpp/qgis_1.0.2 /Users/bhattgopal/dev/cpp/qgis_1.0.2/src/providers/wms /Users/bhattgopal/dev/cpp/qgis_1.0.2/build /Users/bhattgopal/dev/cpp/qgis_1.0.2/build/src/providers/wms /Users/bhattgopal/dev/cpp/qgis_1.0.2/build/src/providers/wms/CMakeFiles/wmsprovider.dir/DependInfo.cmake

# Object files for target wmsprovider
wmsprovider_OBJECTS = \
"CMakeFiles/wmsprovider.dir/qgswmsprovider.o" \
"CMakeFiles/wmsprovider.dir/moc_qgswmsprovider.o"

# External object files for target wmsprovider
wmsprovider_EXTERNAL_OBJECTS =

src/providers/wms/libwmsprovider.so: src/providers/wms/CMakeFiles/wmsprovider.dir/qgswmsprovider.o
src/providers/wms/libwmsprovider.so: src/providers/wms/CMakeFiles/wmsprovider.dir/moc_qgswmsprovider.o
src/providers/wms/libwmsprovider.so: /Library/Frameworks/QtCore.framework
src/providers/wms/libwmsprovider.so: /Library/Frameworks/QtXml.framework
src/providers/wms/libwmsprovider.so: src/core/libqgis_core.dylib
src/providers/wms/libwmsprovider.so: /Library/Frameworks/QtGui.framework
src/providers/wms/libwmsprovider.so: /Library/Frameworks/QtSvg.framework
src/providers/wms/libwmsprovider.so: /Library/Frameworks/QtNetwork.framework
src/providers/wms/libwmsprovider.so: /Library/Frameworks/proj.framework
src/providers/wms/libwmsprovider.so: /Library/Frameworks/GEOS.framework/unix/lib/libgeos_c.dylib
src/providers/wms/libwmsprovider.so: /Library/Frameworks/GDAL.framework/Versions/1.5/unix/lib/libgdal.dylib
src/providers/wms/libwmsprovider.so: /Library/Frameworks/sqlite3.framework
src/providers/wms/libwmsprovider.so: src/providers/wms/CMakeFiles/wmsprovider.dir/build.make
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --red --bold "Linking CXX shared module libwmsprovider.so"
	cd /Users/bhattgopal/dev/cpp/qgis_1.0.2/build/src/providers/wms && $(CMAKE_COMMAND) -P CMakeFiles/wmsprovider.dir/cmake_clean_target.cmake
	cd /Users/bhattgopal/dev/cpp/qgis_1.0.2/build/src/providers/wms && $(CMAKE_COMMAND) -E cmake_link_script CMakeFiles/wmsprovider.dir/link.txt --verbose=$(VERBOSE)

# Rule to build all files generated by this target.
src/providers/wms/CMakeFiles/wmsprovider.dir/build: src/providers/wms/libwmsprovider.so

src/providers/wms/CMakeFiles/wmsprovider.dir/requires: src/providers/wms/CMakeFiles/wmsprovider.dir/qgswmsprovider.o.requires
src/providers/wms/CMakeFiles/wmsprovider.dir/requires: src/providers/wms/CMakeFiles/wmsprovider.dir/moc_qgswmsprovider.o.requires

src/providers/wms/CMakeFiles/wmsprovider.dir/clean:
	cd /Users/bhattgopal/dev/cpp/qgis_1.0.2/build/src/providers/wms && $(CMAKE_COMMAND) -P CMakeFiles/wmsprovider.dir/cmake_clean.cmake
