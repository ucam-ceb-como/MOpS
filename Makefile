# CMAKE generated file: DO NOT EDIT!
# Generated by "Unix Makefiles" Generator, CMake Version 2.8

# Default target executed when no arguments are given to make.
default_target: all
.PHONY : default_target

#=============================================================================
# Special targets provided by cmake.

# Disable implicit rules so canoncical targets will work.
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

# The top-level source directory on which CMake was run.
CMAKE_SOURCE_DIR = /home/anthony/workspace1/LocalMops

# The top-level build directory on which CMake was run.
CMAKE_BINARY_DIR = /home/anthony/workspace1/LocalMops

#=============================================================================
# Targets provided globally by CMake.

# Special rule for the target edit_cache
edit_cache:
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --cyan "Running interactive CMake command-line interface..."
	/usr/bin/cmake -i .
.PHONY : edit_cache

# Special rule for the target edit_cache
edit_cache/fast: edit_cache
.PHONY : edit_cache/fast

# Special rule for the target rebuild_cache
rebuild_cache:
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --cyan "Running CMake to regenerate build system..."
	/usr/bin/cmake -H$(CMAKE_SOURCE_DIR) -B$(CMAKE_BINARY_DIR)
.PHONY : rebuild_cache

# Special rule for the target rebuild_cache
rebuild_cache/fast: rebuild_cache
.PHONY : rebuild_cache/fast

# Special rule for the target test
test:
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --cyan "Running tests..."
	/usr/bin/ctest --force-new-ctest-process $(ARGS)
.PHONY : test

# Special rule for the target test
test/fast: test
.PHONY : test/fast

# The main all target
all: cmake_check_build_system
	$(CMAKE_COMMAND) -E cmake_progress_start /home/anthony/workspace1/LocalMops/CMakeFiles /home/anthony/workspace1/LocalMops/CMakeFiles/progress.marks
	$(MAKE) -f CMakeFiles/Makefile2 all
	$(CMAKE_COMMAND) -E cmake_progress_start /home/anthony/workspace1/LocalMops/CMakeFiles 0
.PHONY : all

# The main clean target
clean:
	$(MAKE) -f CMakeFiles/Makefile2 clean
.PHONY : clean

# The main clean target
clean/fast: clean
.PHONY : clean/fast

# Prepare targets for installation.
preinstall: all
	$(MAKE) -f CMakeFiles/Makefile2 preinstall
.PHONY : preinstall

# Prepare targets for installation.
preinstall/fast:
	$(MAKE) -f CMakeFiles/Makefile2 preinstall
.PHONY : preinstall/fast

# clear depends
depend:
	$(CMAKE_COMMAND) -H$(CMAKE_SOURCE_DIR) -B$(CMAKE_BINARY_DIR) --check-build-system CMakeFiles/Makefile.cmake 1
.PHONY : depend

#=============================================================================
# Target rules for targets named brush-app

# Build rule for target.
brush-app: cmake_check_build_system
	$(MAKE) -f CMakeFiles/Makefile2 brush-app
.PHONY : brush-app

# fast build rule for target.
brush-app/fast:
	$(MAKE) -f CMakeFiles/brush-app.dir/build.make CMakeFiles/brush-app.dir/build
.PHONY : brush-app/fast

#=============================================================================
# Target rules for targets named brush-moonmd-test

# Build rule for target.
brush-moonmd-test: cmake_check_build_system
	$(MAKE) -f CMakeFiles/Makefile2 brush-moonmd-test
.PHONY : brush-moonmd-test

# fast build rule for target.
brush-moonmd-test/fast:
	$(MAKE) -f CMakeFiles/brush-moonmd-test.dir/build.make CMakeFiles/brush-moonmd-test.dir/build
.PHONY : brush-moonmd-test/fast

#=============================================================================
# Target rules for targets named camflow-app

# Build rule for target.
camflow-app: cmake_check_build_system
	$(MAKE) -f CMakeFiles/Makefile2 camflow-app
.PHONY : camflow-app

# fast build rule for target.
camflow-app/fast:
	$(MAKE) -f CMakeFiles/camflow-app.dir/build.make CMakeFiles/camflow-app.dir/build
.PHONY : camflow-app/fast

#=============================================================================
# Target rules for targets named chemkinReader-test

# Build rule for target.
chemkinReader-test: cmake_check_build_system
	$(MAKE) -f CMakeFiles/Makefile2 chemkinReader-test
.PHONY : chemkinReader-test

# fast build rule for target.
chemkinReader-test/fast:
	$(MAKE) -f CMakeFiles/chemkinReader-test.dir/build.make CMakeFiles/chemkinReader-test.dir/build
.PHONY : chemkinReader-test/fast

#=============================================================================
# Target rules for targets named mops-app

# Build rule for target.
mops-app: cmake_check_build_system
	$(MAKE) -f CMakeFiles/Makefile2 mops-app
.PHONY : mops-app

# fast build rule for target.
mops-app/fast:
	$(MAKE) -f CMakeFiles/mops-app.dir/build.make CMakeFiles/mops-app.dir/build
.PHONY : mops-app/fast

#=============================================================================
# Target rules for targets named sprogc-test

# Build rule for target.
sprogc-test: cmake_check_build_system
	$(MAKE) -f CMakeFiles/Makefile2 sprogc-test
.PHONY : sprogc-test

# fast build rule for target.
sprogc-test/fast:
	$(MAKE) -f CMakeFiles/sprogc-test.dir/build.make CMakeFiles/sprogc-test.dir/build
.PHONY : sprogc-test/fast

#=============================================================================
# Target rules for targets named chemkinReader

# Build rule for target.
chemkinReader: cmake_check_build_system
	$(MAKE) -f CMakeFiles/Makefile2 chemkinReader
.PHONY : chemkinReader

# fast build rule for target.
chemkinReader/fast:
	$(MAKE) -f src/io/chemkinReader/CMakeFiles/chemkinReader.dir/build.make src/io/chemkinReader/CMakeFiles/chemkinReader.dir/build
.PHONY : chemkinReader/fast

#=============================================================================
# Target rules for targets named comostrings

# Build rule for target.
comostrings: cmake_check_build_system
	$(MAKE) -f CMakeFiles/Makefile2 comostrings
.PHONY : comostrings

# fast build rule for target.
comostrings/fast:
	$(MAKE) -f src/io/comostrings/CMakeFiles/comostrings.dir/build.make src/io/comostrings/CMakeFiles/comostrings.dir/build
.PHONY : comostrings/fast

#=============================================================================
# Target rules for targets named camxml

# Build rule for target.
camxml: cmake_check_build_system
	$(MAKE) -f CMakeFiles/Makefile2 camxml
.PHONY : camxml

# fast build rule for target.
camxml/fast:
	$(MAKE) -f src/io/camxml/CMakeFiles/camxml.dir/build.make src/io/camxml/CMakeFiles/camxml.dir/build
.PHONY : camxml/fast

#=============================================================================
# Target rules for targets named geometry

# Build rule for target.
geometry: cmake_check_build_system
	$(MAKE) -f CMakeFiles/Makefile2 geometry
.PHONY : geometry

# fast build rule for target.
geometry/fast:
	$(MAKE) -f src/geometry/CMakeFiles/geometry.dir/build.make src/geometry/CMakeFiles/geometry.dir/build
.PHONY : geometry/fast

#=============================================================================
# Target rules for targets named sprog

# Build rule for target.
sprog: cmake_check_build_system
	$(MAKE) -f CMakeFiles/Makefile2 sprog
.PHONY : sprog

# fast build rule for target.
sprog/fast:
	$(MAKE) -f src/sprogc/CMakeFiles/sprog.dir/build.make src/sprogc/CMakeFiles/sprog.dir/build
.PHONY : sprog/fast

#=============================================================================
# Target rules for targets named sweep

# Build rule for target.
sweep: cmake_check_build_system
	$(MAKE) -f CMakeFiles/Makefile2 sweep
.PHONY : sweep

# fast build rule for target.
sweep/fast:
	$(MAKE) -f src/sweepc/CMakeFiles/sweep.dir/build.make src/sweepc/CMakeFiles/sweep.dir/build
.PHONY : sweep/fast

#=============================================================================
# Target rules for targets named cvodes

# Build rule for target.
cvodes: cmake_check_build_system
	$(MAKE) -f CMakeFiles/Makefile2 cvodes
.PHONY : cvodes

# fast build rule for target.
cvodes/fast:
	$(MAKE) -f src/odesolvers/cvodes/CMakeFiles/cvodes.dir/build.make src/odesolvers/cvodes/CMakeFiles/cvodes.dir/build
.PHONY : cvodes/fast

#=============================================================================
# Target rules for targets named cvode

# Build rule for target.
cvode: cmake_check_build_system
	$(MAKE) -f CMakeFiles/Makefile2 cvode
.PHONY : cvode

# fast build rule for target.
cvode/fast:
	$(MAKE) -f src/odesolvers/cvode/CMakeFiles/cvode.dir/build.make src/odesolvers/cvode/CMakeFiles/cvode.dir/build
.PHONY : cvode/fast

#=============================================================================
# Target rules for targets named radau

# Build rule for target.
radau: cmake_check_build_system
	$(MAKE) -f CMakeFiles/Makefile2 radau
.PHONY : radau

# fast build rule for target.
radau/fast:
	$(MAKE) -f src/odesolvers/radau/CMakeFiles/radau.dir/build.make src/odesolvers/radau/CMakeFiles/radau.dir/build
.PHONY : radau/fast

#=============================================================================
# Target rules for targets named camflow

# Build rule for target.
camflow: cmake_check_build_system
	$(MAKE) -f CMakeFiles/Makefile2 camflow
.PHONY : camflow

# fast build rule for target.
camflow/fast:
	$(MAKE) -f src/camflow/CMakeFiles/camflow.dir/build.make src/camflow/CMakeFiles/camflow.dir/build
.PHONY : camflow/fast

#=============================================================================
# Target rules for targets named mops

# Build rule for target.
mops: cmake_check_build_system
	$(MAKE) -f CMakeFiles/Makefile2 mops
.PHONY : mops

# fast build rule for target.
mops/fast:
	$(MAKE) -f src/mopsc/CMakeFiles/mops.dir/build.make src/mopsc/CMakeFiles/mops.dir/build
.PHONY : mops/fast

#=============================================================================
# Target rules for targets named brush

# Build rule for target.
brush: cmake_check_build_system
	$(MAKE) -f CMakeFiles/Makefile2 brush
.PHONY : brush

# fast build rule for target.
brush/fast:
	$(MAKE) -f src/brush/CMakeFiles/brush.dir/build.make src/brush/CMakeFiles/brush.dir/build
.PHONY : brush/fast

applications/solvers/brush/brush.o: applications/solvers/brush/brush.cpp.o
.PHONY : applications/solvers/brush/brush.o

# target to build an object file
applications/solvers/brush/brush.cpp.o:
	$(MAKE) -f CMakeFiles/brush-app.dir/build.make CMakeFiles/brush-app.dir/applications/solvers/brush/brush.cpp.o
.PHONY : applications/solvers/brush/brush.cpp.o

applications/solvers/brush/brush.i: applications/solvers/brush/brush.cpp.i
.PHONY : applications/solvers/brush/brush.i

# target to preprocess a source file
applications/solvers/brush/brush.cpp.i:
	$(MAKE) -f CMakeFiles/brush-app.dir/build.make CMakeFiles/brush-app.dir/applications/solvers/brush/brush.cpp.i
.PHONY : applications/solvers/brush/brush.cpp.i

applications/solvers/brush/brush.s: applications/solvers/brush/brush.cpp.s
.PHONY : applications/solvers/brush/brush.s

# target to generate assembly for a file
applications/solvers/brush/brush.cpp.s:
	$(MAKE) -f CMakeFiles/brush-app.dir/build.make CMakeFiles/brush-app.dir/applications/solvers/brush/brush.cpp.s
.PHONY : applications/solvers/brush/brush.cpp.s

applications/solvers/camflow/cam_kernel.o: applications/solvers/camflow/cam_kernel.cpp.o
.PHONY : applications/solvers/camflow/cam_kernel.o

# target to build an object file
applications/solvers/camflow/cam_kernel.cpp.o:
	$(MAKE) -f CMakeFiles/camflow-app.dir/build.make CMakeFiles/camflow-app.dir/applications/solvers/camflow/cam_kernel.cpp.o
.PHONY : applications/solvers/camflow/cam_kernel.cpp.o

applications/solvers/camflow/cam_kernel.i: applications/solvers/camflow/cam_kernel.cpp.i
.PHONY : applications/solvers/camflow/cam_kernel.i

# target to preprocess a source file
applications/solvers/camflow/cam_kernel.cpp.i:
	$(MAKE) -f CMakeFiles/camflow-app.dir/build.make CMakeFiles/camflow-app.dir/applications/solvers/camflow/cam_kernel.cpp.i
.PHONY : applications/solvers/camflow/cam_kernel.cpp.i

applications/solvers/camflow/cam_kernel.s: applications/solvers/camflow/cam_kernel.cpp.s
.PHONY : applications/solvers/camflow/cam_kernel.s

# target to generate assembly for a file
applications/solvers/camflow/cam_kernel.cpp.s:
	$(MAKE) -f CMakeFiles/camflow-app.dir/build.make CMakeFiles/camflow-app.dir/applications/solvers/camflow/cam_kernel.cpp.s
.PHONY : applications/solvers/camflow/cam_kernel.cpp.s

applications/solvers/mopsc/mops.o: applications/solvers/mopsc/mops.cpp.o
.PHONY : applications/solvers/mopsc/mops.o

# target to build an object file
applications/solvers/mopsc/mops.cpp.o:
	$(MAKE) -f CMakeFiles/mops-app.dir/build.make CMakeFiles/mops-app.dir/applications/solvers/mopsc/mops.cpp.o
.PHONY : applications/solvers/mopsc/mops.cpp.o

applications/solvers/mopsc/mops.i: applications/solvers/mopsc/mops.cpp.i
.PHONY : applications/solvers/mopsc/mops.i

# target to preprocess a source file
applications/solvers/mopsc/mops.cpp.i:
	$(MAKE) -f CMakeFiles/mops-app.dir/build.make CMakeFiles/mops-app.dir/applications/solvers/mopsc/mops.cpp.i
.PHONY : applications/solvers/mopsc/mops.cpp.i

applications/solvers/mopsc/mops.s: applications/solvers/mopsc/mops.cpp.s
.PHONY : applications/solvers/mopsc/mops.s

# target to generate assembly for a file
applications/solvers/mopsc/mops.cpp.s:
	$(MAKE) -f CMakeFiles/mops-app.dir/build.make CMakeFiles/mops-app.dir/applications/solvers/mopsc/mops.cpp.s
.PHONY : applications/solvers/mopsc/mops.cpp.s

applications/test-harnesses/brush/moonmd_test.o: applications/test-harnesses/brush/moonmd_test.cpp.o
.PHONY : applications/test-harnesses/brush/moonmd_test.o

# target to build an object file
applications/test-harnesses/brush/moonmd_test.cpp.o:
	$(MAKE) -f CMakeFiles/brush-moonmd-test.dir/build.make CMakeFiles/brush-moonmd-test.dir/applications/test-harnesses/brush/moonmd_test.cpp.o
.PHONY : applications/test-harnesses/brush/moonmd_test.cpp.o

applications/test-harnesses/brush/moonmd_test.i: applications/test-harnesses/brush/moonmd_test.cpp.i
.PHONY : applications/test-harnesses/brush/moonmd_test.i

# target to preprocess a source file
applications/test-harnesses/brush/moonmd_test.cpp.i:
	$(MAKE) -f CMakeFiles/brush-moonmd-test.dir/build.make CMakeFiles/brush-moonmd-test.dir/applications/test-harnesses/brush/moonmd_test.cpp.i
.PHONY : applications/test-harnesses/brush/moonmd_test.cpp.i

applications/test-harnesses/brush/moonmd_test.s: applications/test-harnesses/brush/moonmd_test.cpp.s
.PHONY : applications/test-harnesses/brush/moonmd_test.s

# target to generate assembly for a file
applications/test-harnesses/brush/moonmd_test.cpp.s:
	$(MAKE) -f CMakeFiles/brush-moonmd-test.dir/build.make CMakeFiles/brush-moonmd-test.dir/applications/test-harnesses/brush/moonmd_test.cpp.s
.PHONY : applications/test-harnesses/brush/moonmd_test.cpp.s

applications/test-harnesses/chemkinReader/chemkinReaderTest.o: applications/test-harnesses/chemkinReader/chemkinReaderTest.cpp.o
.PHONY : applications/test-harnesses/chemkinReader/chemkinReaderTest.o

# target to build an object file
applications/test-harnesses/chemkinReader/chemkinReaderTest.cpp.o:
	$(MAKE) -f CMakeFiles/chemkinReader-test.dir/build.make CMakeFiles/chemkinReader-test.dir/applications/test-harnesses/chemkinReader/chemkinReaderTest.cpp.o
.PHONY : applications/test-harnesses/chemkinReader/chemkinReaderTest.cpp.o

applications/test-harnesses/chemkinReader/chemkinReaderTest.i: applications/test-harnesses/chemkinReader/chemkinReaderTest.cpp.i
.PHONY : applications/test-harnesses/chemkinReader/chemkinReaderTest.i

# target to preprocess a source file
applications/test-harnesses/chemkinReader/chemkinReaderTest.cpp.i:
	$(MAKE) -f CMakeFiles/chemkinReader-test.dir/build.make CMakeFiles/chemkinReader-test.dir/applications/test-harnesses/chemkinReader/chemkinReaderTest.cpp.i
.PHONY : applications/test-harnesses/chemkinReader/chemkinReaderTest.cpp.i

applications/test-harnesses/chemkinReader/chemkinReaderTest.s: applications/test-harnesses/chemkinReader/chemkinReaderTest.cpp.s
.PHONY : applications/test-harnesses/chemkinReader/chemkinReaderTest.s

# target to generate assembly for a file
applications/test-harnesses/chemkinReader/chemkinReaderTest.cpp.s:
	$(MAKE) -f CMakeFiles/chemkinReader-test.dir/build.make CMakeFiles/chemkinReader-test.dir/applications/test-harnesses/chemkinReader/chemkinReaderTest.cpp.s
.PHONY : applications/test-harnesses/chemkinReader/chemkinReaderTest.cpp.s

applications/test-harnesses/sprogc/sprogc_test.o: applications/test-harnesses/sprogc/sprogc_test.cpp.o
.PHONY : applications/test-harnesses/sprogc/sprogc_test.o

# target to build an object file
applications/test-harnesses/sprogc/sprogc_test.cpp.o:
	$(MAKE) -f CMakeFiles/sprogc-test.dir/build.make CMakeFiles/sprogc-test.dir/applications/test-harnesses/sprogc/sprogc_test.cpp.o
.PHONY : applications/test-harnesses/sprogc/sprogc_test.cpp.o

applications/test-harnesses/sprogc/sprogc_test.i: applications/test-harnesses/sprogc/sprogc_test.cpp.i
.PHONY : applications/test-harnesses/sprogc/sprogc_test.i

# target to preprocess a source file
applications/test-harnesses/sprogc/sprogc_test.cpp.i:
	$(MAKE) -f CMakeFiles/sprogc-test.dir/build.make CMakeFiles/sprogc-test.dir/applications/test-harnesses/sprogc/sprogc_test.cpp.i
.PHONY : applications/test-harnesses/sprogc/sprogc_test.cpp.i

applications/test-harnesses/sprogc/sprogc_test.s: applications/test-harnesses/sprogc/sprogc_test.cpp.s
.PHONY : applications/test-harnesses/sprogc/sprogc_test.s

# target to generate assembly for a file
applications/test-harnesses/sprogc/sprogc_test.cpp.s:
	$(MAKE) -f CMakeFiles/sprogc-test.dir/build.make CMakeFiles/sprogc-test.dir/applications/test-harnesses/sprogc/sprogc_test.cpp.s
.PHONY : applications/test-harnesses/sprogc/sprogc_test.cpp.s

# Help Target
help:
	@echo "The following are some of the valid targets for this Makefile:"
	@echo "... all (the default if no target is provided)"
	@echo "... clean"
	@echo "... depend"
	@echo "... brush-app"
	@echo "... brush-moonmd-test"
	@echo "... camflow-app"
	@echo "... chemkinReader-test"
	@echo "... edit_cache"
	@echo "... mops-app"
	@echo "... rebuild_cache"
	@echo "... sprogc-test"
	@echo "... test"
	@echo "... chemkinReader"
	@echo "... comostrings"
	@echo "... camxml"
	@echo "... geometry"
	@echo "... sprog"
	@echo "... sweep"
	@echo "... cvodes"
	@echo "... cvode"
	@echo "... radau"
	@echo "... camflow"
	@echo "... mops"
	@echo "... brush"
	@echo "... applications/solvers/brush/brush.o"
	@echo "... applications/solvers/brush/brush.i"
	@echo "... applications/solvers/brush/brush.s"
	@echo "... applications/solvers/camflow/cam_kernel.o"
	@echo "... applications/solvers/camflow/cam_kernel.i"
	@echo "... applications/solvers/camflow/cam_kernel.s"
	@echo "... applications/solvers/mopsc/mops.o"
	@echo "... applications/solvers/mopsc/mops.i"
	@echo "... applications/solvers/mopsc/mops.s"
	@echo "... applications/test-harnesses/brush/moonmd_test.o"
	@echo "... applications/test-harnesses/brush/moonmd_test.i"
	@echo "... applications/test-harnesses/brush/moonmd_test.s"
	@echo "... applications/test-harnesses/chemkinReader/chemkinReaderTest.o"
	@echo "... applications/test-harnesses/chemkinReader/chemkinReaderTest.i"
	@echo "... applications/test-harnesses/chemkinReader/chemkinReaderTest.s"
	@echo "... applications/test-harnesses/sprogc/sprogc_test.o"
	@echo "... applications/test-harnesses/sprogc/sprogc_test.i"
	@echo "... applications/test-harnesses/sprogc/sprogc_test.s"
.PHONY : help



#=============================================================================
# Special targets to cleanup operation of make.

# Special rule to run CMake to check the build system integrity.
# No rule that depends on this can have commands that come from listfiles
# because they might be regenerated.
cmake_check_build_system:
	$(CMAKE_COMMAND) -H$(CMAKE_SOURCE_DIR) -B$(CMAKE_BINARY_DIR) --check-build-system CMakeFiles/Makefile.cmake 0
.PHONY : cmake_check_build_system

