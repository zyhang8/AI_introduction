# CMAKE generated file: DO NOT EDIT!
# Generated by "MinGW Makefiles" Generator, CMake Version 3.15

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

SHELL = cmd.exe

# The CMake executable.
CMAKE_COMMAND = "L:\software\clion\CLion 2019.2.5\bin\cmake\win\bin\cmake.exe"

# The command to remove a file.
RM = "L:\software\clion\CLion 2019.2.5\bin\cmake\win\bin\cmake.exe" -E remove -f

# Escaping for special characters.
EQUALS = =

# The top-level source directory on which CMake was run.
CMAKE_SOURCE_DIR = L:\github\workspace\AI_introduction

# The top-level build directory on which CMake was run.
CMAKE_BINARY_DIR = L:\github\workspace\AI_introduction\cmake-build-debug

# Include any dependencies generated for this target.
include CMakeFiles/ruadhan.dir/depend.make

# Include the progress variables for this target.
include CMakeFiles/ruadhan.dir/progress.make

# Include the compile flags for this target's objects.
include CMakeFiles/ruadhan.dir/flags.make

CMakeFiles/ruadhan.dir/astar.cpp.obj: CMakeFiles/ruadhan.dir/flags.make
CMakeFiles/ruadhan.dir/astar.cpp.obj: ../astar.cpp
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=L:\github\workspace\AI_introduction\cmake-build-debug\CMakeFiles --progress-num=$(CMAKE_PROGRESS_1) "Building CXX object CMakeFiles/ruadhan.dir/astar.cpp.obj"
	L:\software\mingw\bin\g++.exe  $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -o CMakeFiles\ruadhan.dir\astar.cpp.obj -c L:\github\workspace\AI_introduction\astar.cpp

CMakeFiles/ruadhan.dir/astar.cpp.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/ruadhan.dir/astar.cpp.i"
	L:\software\mingw\bin\g++.exe $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E L:\github\workspace\AI_introduction\astar.cpp > CMakeFiles\ruadhan.dir\astar.cpp.i

CMakeFiles/ruadhan.dir/astar.cpp.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/ruadhan.dir/astar.cpp.s"
	L:\software\mingw\bin\g++.exe $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S L:\github\workspace\AI_introduction\astar.cpp -o CMakeFiles\ruadhan.dir\astar.cpp.s

# Object files for target ruadhan
ruadhan_OBJECTS = \
"CMakeFiles/ruadhan.dir/astar.cpp.obj"

# External object files for target ruadhan
ruadhan_EXTERNAL_OBJECTS =

ruadhan.exe: CMakeFiles/ruadhan.dir/astar.cpp.obj
ruadhan.exe: CMakeFiles/ruadhan.dir/build.make
ruadhan.exe: CMakeFiles/ruadhan.dir/linklibs.rsp
ruadhan.exe: CMakeFiles/ruadhan.dir/objects1.rsp
ruadhan.exe: CMakeFiles/ruadhan.dir/link.txt
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --bold --progress-dir=L:\github\workspace\AI_introduction\cmake-build-debug\CMakeFiles --progress-num=$(CMAKE_PROGRESS_2) "Linking CXX executable ruadhan.exe"
	$(CMAKE_COMMAND) -E cmake_link_script CMakeFiles\ruadhan.dir\link.txt --verbose=$(VERBOSE)

# Rule to build all files generated by this target.
CMakeFiles/ruadhan.dir/build: ruadhan.exe

.PHONY : CMakeFiles/ruadhan.dir/build

CMakeFiles/ruadhan.dir/clean:
	$(CMAKE_COMMAND) -P CMakeFiles\ruadhan.dir\cmake_clean.cmake
.PHONY : CMakeFiles/ruadhan.dir/clean

CMakeFiles/ruadhan.dir/depend:
	$(CMAKE_COMMAND) -E cmake_depends "MinGW Makefiles" L:\github\workspace\AI_introduction L:\github\workspace\AI_introduction L:\github\workspace\AI_introduction\cmake-build-debug L:\github\workspace\AI_introduction\cmake-build-debug L:\github\workspace\AI_introduction\cmake-build-debug\CMakeFiles\ruadhan.dir\DependInfo.cmake --color=$(COLOR)
.PHONY : CMakeFiles/ruadhan.dir/depend

