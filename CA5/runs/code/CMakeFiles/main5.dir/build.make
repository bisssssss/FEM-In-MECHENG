# CMAKE generated file: DO NOT EDIT!
# Generated by "Unix Makefiles" Generator, CMake Version 2.8

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

# The program to use to edit the cache.
CMAKE_EDIT_COMMAND = /usr/bin/ccmake

# The top-level source directory on which CMake was run.
CMAKE_SOURCE_DIR = /home/dianawww/ME505/CA5/runs/code

# The top-level build directory on which CMake was run.
CMAKE_BINARY_DIR = /home/dianawww/ME505/CA5/runs/code

# Include any dependencies generated for this target.
include CMakeFiles/main5.dir/depend.make

# Include the progress variables for this target.
include CMakeFiles/main5.dir/progress.make

# Include the compile flags for this target's objects.
include CMakeFiles/main5.dir/flags.make

CMakeFiles/main5.dir/main5.cc.o: CMakeFiles/main5.dir/flags.make
CMakeFiles/main5.dir/main5.cc.o: main5.cc
	$(CMAKE_COMMAND) -E cmake_progress_report /home/dianawww/ME505/CA5/runs/code/CMakeFiles $(CMAKE_PROGRESS_1)
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Building CXX object CMakeFiles/main5.dir/main5.cc.o"
	/usr/caen/intel-2018/bin/icpc   $(CXX_DEFINES) $(CXX_FLAGS) -o CMakeFiles/main5.dir/main5.cc.o -c /home/dianawww/ME505/CA5/runs/code/main5.cc

CMakeFiles/main5.dir/main5.cc.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/main5.dir/main5.cc.i"
	/usr/caen/intel-2018/bin/icpc  $(CXX_DEFINES) $(CXX_FLAGS) -E /home/dianawww/ME505/CA5/runs/code/main5.cc > CMakeFiles/main5.dir/main5.cc.i

CMakeFiles/main5.dir/main5.cc.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/main5.dir/main5.cc.s"
	/usr/caen/intel-2018/bin/icpc  $(CXX_DEFINES) $(CXX_FLAGS) -S /home/dianawww/ME505/CA5/runs/code/main5.cc -o CMakeFiles/main5.dir/main5.cc.s

CMakeFiles/main5.dir/main5.cc.o.requires:
.PHONY : CMakeFiles/main5.dir/main5.cc.o.requires

CMakeFiles/main5.dir/main5.cc.o.provides: CMakeFiles/main5.dir/main5.cc.o.requires
	$(MAKE) -f CMakeFiles/main5.dir/build.make CMakeFiles/main5.dir/main5.cc.o.provides.build
.PHONY : CMakeFiles/main5.dir/main5.cc.o.provides

CMakeFiles/main5.dir/main5.cc.o.provides.build: CMakeFiles/main5.dir/main5.cc.o

CMakeFiles/main5.dir/FEM5.cc.o: CMakeFiles/main5.dir/flags.make
CMakeFiles/main5.dir/FEM5.cc.o: FEM5.cc
	$(CMAKE_COMMAND) -E cmake_progress_report /home/dianawww/ME505/CA5/runs/code/CMakeFiles $(CMAKE_PROGRESS_2)
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Building CXX object CMakeFiles/main5.dir/FEM5.cc.o"
	/usr/caen/intel-2018/bin/icpc   $(CXX_DEFINES) $(CXX_FLAGS) -o CMakeFiles/main5.dir/FEM5.cc.o -c /home/dianawww/ME505/CA5/runs/code/FEM5.cc

CMakeFiles/main5.dir/FEM5.cc.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/main5.dir/FEM5.cc.i"
	/usr/caen/intel-2018/bin/icpc  $(CXX_DEFINES) $(CXX_FLAGS) -E /home/dianawww/ME505/CA5/runs/code/FEM5.cc > CMakeFiles/main5.dir/FEM5.cc.i

CMakeFiles/main5.dir/FEM5.cc.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/main5.dir/FEM5.cc.s"
	/usr/caen/intel-2018/bin/icpc  $(CXX_DEFINES) $(CXX_FLAGS) -S /home/dianawww/ME505/CA5/runs/code/FEM5.cc -o CMakeFiles/main5.dir/FEM5.cc.s

CMakeFiles/main5.dir/FEM5.cc.o.requires:
.PHONY : CMakeFiles/main5.dir/FEM5.cc.o.requires

CMakeFiles/main5.dir/FEM5.cc.o.provides: CMakeFiles/main5.dir/FEM5.cc.o.requires
	$(MAKE) -f CMakeFiles/main5.dir/build.make CMakeFiles/main5.dir/FEM5.cc.o.provides.build
.PHONY : CMakeFiles/main5.dir/FEM5.cc.o.provides

CMakeFiles/main5.dir/FEM5.cc.o.provides.build: CMakeFiles/main5.dir/FEM5.cc.o

# Object files for target main5
main5_OBJECTS = \
"CMakeFiles/main5.dir/main5.cc.o" \
"CMakeFiles/main5.dir/FEM5.cc.o"

# External object files for target main5
main5_EXTERNAL_OBJECTS =

main5: CMakeFiles/main5.dir/main5.cc.o
main5: CMakeFiles/main5.dir/FEM5.cc.o
main5: CMakeFiles/main5.dir/build.make
main5: /usr/um/dealii-9.0.0/lib/libdeal_II.g.so.9.0.0
main5: /usr/caen/intel-2018/mkl/lib/intel64/libmkl_intel_lp64.a
main5: /usr/caen/intel-2018/mkl/lib/intel64/libmkl_sequential.a
main5: /usr/caen/intel-2018/mkl/lib/intel64/libmkl_core.a
main5: /usr/lib64/libz.so
main5: /usr/lib64/libgsl.so
main5: /usr/lib64/libgslcblas.so
main5: /usr/lib64/libhdf5_hl.so
main5: /usr/lib64/libhdf5.so
main5: CMakeFiles/main5.dir/link.txt
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --red --bold "Linking CXX executable main5"
	$(CMAKE_COMMAND) -E cmake_link_script CMakeFiles/main5.dir/link.txt --verbose=$(VERBOSE)

# Rule to build all files generated by this target.
CMakeFiles/main5.dir/build: main5
.PHONY : CMakeFiles/main5.dir/build

CMakeFiles/main5.dir/requires: CMakeFiles/main5.dir/main5.cc.o.requires
CMakeFiles/main5.dir/requires: CMakeFiles/main5.dir/FEM5.cc.o.requires
.PHONY : CMakeFiles/main5.dir/requires

CMakeFiles/main5.dir/clean:
	$(CMAKE_COMMAND) -P CMakeFiles/main5.dir/cmake_clean.cmake
.PHONY : CMakeFiles/main5.dir/clean

CMakeFiles/main5.dir/depend:
	cd /home/dianawww/ME505/CA5/runs/code && $(CMAKE_COMMAND) -E cmake_depends "Unix Makefiles" /home/dianawww/ME505/CA5/runs/code /home/dianawww/ME505/CA5/runs/code /home/dianawww/ME505/CA5/runs/code /home/dianawww/ME505/CA5/runs/code /home/dianawww/ME505/CA5/runs/code/CMakeFiles/main5.dir/DependInfo.cmake --color=$(COLOR)
.PHONY : CMakeFiles/main5.dir/depend

