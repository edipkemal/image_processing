# CMAKE generated file: DO NOT EDIT!
# Generated by "Unix Makefiles" Generator, CMake Version 3.20

# Delete rule output on recipe failure.
.DELETE_ON_ERROR:

#=============================================================================
# Special targets provided by cmake.

# Disable implicit rules so canonical targets will work.
.SUFFIXES:

# Disable VCS-based implicit rules.
% : %,v

# Disable VCS-based implicit rules.
% : RCS/%

# Disable VCS-based implicit rules.
% : RCS/%,v

# Disable VCS-based implicit rules.
% : SCCS/s.%

# Disable VCS-based implicit rules.
% : s.%

.SUFFIXES: .hpux_make_needs_suffix_list

# Command-line flag to silence nested $(MAKE).
$(VERBOSE)MAKESILENT = -s

#Suppress display of executed commands.
$(VERBOSE).SILENT:

# A target that is always out of date.
cmake_force:
.PHONY : cmake_force

#=============================================================================
# Set environment variables for the build.

# The shell in which to execute make rules.
SHELL = /bin/sh

# The CMake executable.
CMAKE_COMMAND = /usr/local/bin/cmake

# The command to remove a file.
RM = /usr/local/bin/cmake -E rm -f

# Escaping for special characters.
EQUALS = =

# The top-level source directory on which CMake was run.
CMAKE_SOURCE_DIR = /home/edip/hw2/v5/libceng391

# The top-level build directory on which CMake was run.
CMAKE_BINARY_DIR = /home/edip/hw2/v5/libceng391/build

# Include any dependencies generated for this target.
include app/CMakeFiles/image-box-filter.dir/depend.make
# Include any dependencies generated by the compiler for this target.
include app/CMakeFiles/image-box-filter.dir/compiler_depend.make

# Include the progress variables for this target.
include app/CMakeFiles/image-box-filter.dir/progress.make

# Include the compile flags for this target's objects.
include app/CMakeFiles/image-box-filter.dir/flags.make

app/CMakeFiles/image-box-filter.dir/image_box_filter.cpp.o: app/CMakeFiles/image-box-filter.dir/flags.make
app/CMakeFiles/image-box-filter.dir/image_box_filter.cpp.o: ../app/image_box_filter.cpp
app/CMakeFiles/image-box-filter.dir/image_box_filter.cpp.o: app/CMakeFiles/image-box-filter.dir/compiler_depend.ts
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/home/edip/hw2/v5/libceng391/build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_1) "Building CXX object app/CMakeFiles/image-box-filter.dir/image_box_filter.cpp.o"
	cd /home/edip/hw2/v5/libceng391/build/app && /usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -MD -MT app/CMakeFiles/image-box-filter.dir/image_box_filter.cpp.o -MF CMakeFiles/image-box-filter.dir/image_box_filter.cpp.o.d -o CMakeFiles/image-box-filter.dir/image_box_filter.cpp.o -c /home/edip/hw2/v5/libceng391/app/image_box_filter.cpp

app/CMakeFiles/image-box-filter.dir/image_box_filter.cpp.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/image-box-filter.dir/image_box_filter.cpp.i"
	cd /home/edip/hw2/v5/libceng391/build/app && /usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E /home/edip/hw2/v5/libceng391/app/image_box_filter.cpp > CMakeFiles/image-box-filter.dir/image_box_filter.cpp.i

app/CMakeFiles/image-box-filter.dir/image_box_filter.cpp.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/image-box-filter.dir/image_box_filter.cpp.s"
	cd /home/edip/hw2/v5/libceng391/build/app && /usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S /home/edip/hw2/v5/libceng391/app/image_box_filter.cpp -o CMakeFiles/image-box-filter.dir/image_box_filter.cpp.s

# Object files for target image-box-filter
image__box__filter_OBJECTS = \
"CMakeFiles/image-box-filter.dir/image_box_filter.cpp.o"

# External object files for target image-box-filter
image__box__filter_EXTERNAL_OBJECTS =

app/image-box-filter: app/CMakeFiles/image-box-filter.dir/image_box_filter.cpp.o
app/image-box-filter: app/CMakeFiles/image-box-filter.dir/build.make
app/image-box-filter: src/libceng391.a
app/image-box-filter: /usr/lib/x86_64-linux-gnu/libpng.so
app/image-box-filter: /usr/lib/x86_64-linux-gnu/libz.so
app/image-box-filter: app/CMakeFiles/image-box-filter.dir/link.txt
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --bold --progress-dir=/home/edip/hw2/v5/libceng391/build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_2) "Linking CXX executable image-box-filter"
	cd /home/edip/hw2/v5/libceng391/build/app && $(CMAKE_COMMAND) -E cmake_link_script CMakeFiles/image-box-filter.dir/link.txt --verbose=$(VERBOSE)

# Rule to build all files generated by this target.
app/CMakeFiles/image-box-filter.dir/build: app/image-box-filter
.PHONY : app/CMakeFiles/image-box-filter.dir/build

app/CMakeFiles/image-box-filter.dir/clean:
	cd /home/edip/hw2/v5/libceng391/build/app && $(CMAKE_COMMAND) -P CMakeFiles/image-box-filter.dir/cmake_clean.cmake
.PHONY : app/CMakeFiles/image-box-filter.dir/clean

app/CMakeFiles/image-box-filter.dir/depend:
	cd /home/edip/hw2/v5/libceng391/build && $(CMAKE_COMMAND) -E cmake_depends "Unix Makefiles" /home/edip/hw2/v5/libceng391 /home/edip/hw2/v5/libceng391/app /home/edip/hw2/v5/libceng391/build /home/edip/hw2/v5/libceng391/build/app /home/edip/hw2/v5/libceng391/build/app/CMakeFiles/image-box-filter.dir/DependInfo.cmake --color=$(COLOR)
.PHONY : app/CMakeFiles/image-box-filter.dir/depend
