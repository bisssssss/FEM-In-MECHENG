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
CMAKE_SOURCE_DIR = /home/dianawww/ME505/CA2/runs/CA2a/code

# The top-level build directory on which CMake was run.
CMAKE_BINARY_DIR = /home/dianawww/ME505/CA2/runs/CA2a/code

# Include any dependencies generated for this target.
include CMakeFiles/main2a.dir/depend.make

# Include the progress variables for this target.
include CMakeFiles/main2a.dir/progress.make

# Include the compile flags for this target's objects.
include CMakeFiles/main2a.dir/flags.make

CMakeFiles/main2a.dir/main2a.cc.o: CMakeFiles/main2a.dir/flags.make
CMakeFiles/main2a.dir/main2a.cc.o: main2a.cc
	$(CMAKE_COMMAND) -E cmake_progress_report /home/dianawww/ME505/CA2/runs/CA2a/code/CMakeFiles $(CMAKE_PROGRESS_1)
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Building CXX object CMakeFiles/main2a.dir/main2a.cc.o"
	/usr/caen/intel-2018/bin/icpc   $(CXX_DEFINES) $(CXX_FLAGS) -o CMakeFiles/main2a.dir/main2a.cc.o -c /home/dianawww/ME505/CA2/runs/CA2a/code/main2a.cc

CMakeFiles/main2a.dir/main2a.cc.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/main2a.dir/main2a.cc.i"
	/usr/caen/intel-2018/bin/icpc  $(CXX_DEFINES) $(CXX_FLAGS) -E /home/dianawww/ME505/CA2/runs/CA2a/code/main2a.cc > CMakeFiles/main2a.dir/main2a.cc.i

CMakeFiles/main2a.dir/main2a.cc.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/main2a.dir/main2a.cc.s"
	/usr/caen/intel-2018/bin/icpc  $(CXX_DEFINES) $(CXX_FLAGS) -S /home/dianawww/ME505/CA2/runs/CA2a/code/main2a.cc -o CMakeFiles/main2a.dir/main2a.cc.s

CMakeFiles/main2a.dir/main2a.cc.o.requires:
.PHONY : CMakeFiles/main2a.dir/main2a.cc.o.requires

CMakeFiles/main2a.dir/main2a.cc.o.provides: CMakeFiles/main2a.dir/main2a.cc.o.requires
	$(MAKE) -f CMakeFiles/main2a.dir/build.make CMakeFiles/main2a.dir/main2a.cc.o.provides.build
.PHONY : CMakeFiles/main2a.dir/main2a.cc.o.provides

CMakeFiles/main2a.dir/main2a.cc.o.provides.build: CMakeFiles/main2a.dir/main2a.cc.o

CMakeFiles/main2a.dir/FEM2a.cc.o: CMakeFiles/main2a.dir/flags.make
CMakeFiles/main2a.dir/FEM2a.cc.o: FEM2a.cc
	$(CMAKE_COMMAND) -E cmake_progress_report /home/dianawww/ME505/CA2/runs/CA2a/code/CMakeFiles $(CMAKE_PROGRESS_2)
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Building CXX object CMakeFiles/main2a.dir/FEM2a.cc.o"
	/usr/caen/intel-2018/bin/icpc   $(CXX_DEFINES) $(CXX_FLAGS) -o CMakeFiles/main2a.dir/FEM2a.cc.o -c /home/dianawww/ME505/CA2/runs/CA2a/code/FEM2a.cc

CMakeFiles/main2a.dir/FEM2a.cc.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/main2a.dir/FEM2a.cc.i"
	/usr/caen/intel-2018/bin/icpc  $(CXX_DEFINES) $(CXX_FLAGS) -E /home/dianawww/ME505/CA2/runs/CA2a/code/FEM2a.cc > CMakeFiles/main2a.dir/FEM2a.cc.i

CMakeFiles/main2a.dir/FEM2a.cc.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/main2a.dir/FEM2a.cc.s"
	/usr/caen/intel-2018/bin/icpc  $(CXX_DEFINES) $(CXX_FLAGS) -S /home/dianawww/ME505/CA2/runs/CA2a/code/FEM2a.cc -o CMakeFiles/main2a.dir/FEM2a.cc.s

CMakeFiles/main2a.dir/FEM2a.cc.o.requires:
.PHONY : CMakeFiles/main2a.dir/FEM2a.cc.o.requires

CMakeFiles/main2a.dir/FEM2a.cc.o.provides: CMakeFiles/main2a.dir/FEM2a.cc.o.requires
	$(MAKE) -f CMakeFiles/main2a.dir/build.make CMakeFiles/main2a.dir/FEM2a.cc.o.provides.build
.PHONY : CMakeFiles/main2a.dir/FEM2a.cc.o.provides

CMakeFiles/main2a.dir/FEM2a.cc.o.provides.build: CMakeFiles/main2a.dir/FEM2a.cc.o

# Object files for target main2a
main2a_OBJECTS = \
"CMakeFiles/main2a.dir/main2a.cc.o" \
"CMakeFiles/main2a.dir/FEM2a.cc.o"

# External object files for target main2a
main2a_EXTERNAL_OBJECTS =

main2a: CMakeFiles/main2a.dir/main2a.cc.o
main2a: CMakeFiles/main2a.dir/FEM2a.cc.o
main2a: CMakeFiles/main2a.dir/build.make
main2a: /usr/um/dealii-9.0.0/lib/libdeal_II.so.9.0.0
main2a: /usr/caen/intel-2018/mkl/lib/intel64/libmkl_intel_lp64.a
main2a: /usr/caen/intel-2018/mkl/lib/intel64/libmkl_sequential.a
main2a: /usr/caen/intel-2018/mkl/lib/intel64/libmkl_core.a
main2a: /usr/lib64/libz.so
main2a: /usr/lib64/libgsl.so
main2a: /usr/lib64/libgslcblas.so
main2a: /usr/lib64/libhdf5_hl.so
main2a: /usr/lib64/libhdf5.so
main2a: CMakeFiles/main2a.dir/link.txt
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --red --bold "Linking CXX executable main2a"
	$(CMAKE_COMMAND) -E cmake_link_script CMakeFiles/main2a.dir/link.txt --verbose=$(VERBOSE)

# Rule to build all files generated by this target.
CMakeFiles/main2a.dir/build: main2a
.PHONY : CMakeFiles/main2a.dir/build

CMakeFiles/main2a.dir/requires: CMakeFiles/main2a.dir/main2a.cc.o.requires
CMakeFiles/main2a.dir/requires: CMakeFiles/main2a.dir/FEM2a.cc.o.requires
.PHONY : CMakeFiles/main2a.dir/requires

CMakeFiles/main2a.dir/clean:
	$(CMAKE_COMMAND) -P CMakeFiles/main2a.dir/cmake_clean.cmake
.PHONY : CMakeFiles/main2a.dir/clean

CMakeFiles/main2a.dir/depend:
	cd /home/dianawww/ME505/CA2/runs/CA2a/code && $(CMAKE_COMMAND) -E cmake_depends "Unix Makefiles" /home/dianawww/ME505/CA2/runs/CA2a/code /home/dianawww/ME505/CA2/runs/CA2a/code /home/dianawww/ME505/CA2/runs/CA2a/code /home/dianawww/ME505/CA2/runs/CA2a/code /home/dianawww/ME505/CA2/runs/CA2a/code/CMakeFiles/main2a.dir/DependInfo.cmake --color=$(COLOR)
.PHONY : CMakeFiles/main2a.dir/depend
