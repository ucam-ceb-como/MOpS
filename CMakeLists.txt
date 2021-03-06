# The name of this project is "MOPSSUITE". CMakeLists files in this project can
# refer to the root source directory of the project as ${MOPSSUITE_SOURCE_DIR} and
# to the root binary directory of the project as ${MOPSSUITE_BINARY_DIR}.
cmake_minimum_required (VERSION 2.6)
project (MOPSSUITE)

IF(WIN32)
  SET(Boost_ADDITIONAL_VERSIONS "1.47" "1.47.0" "1.49" "1.49.0" "1.51" "1.51.0" "1.57.0")
  SET(BOOST_ROOT "C:\Users\gl413\Downloads\boost_1_57_0")
  SET(Boost_USE_STATIC_LIBS TRUE)
  SET(MSVC_BUILD_DIR "/\${CTEST_CONFIGURATION_TYPE}")
  SET(SHELL_INTERPRETER "bash")
  SET(CMAKE_MODULE_PATH ${PROJECT_BINARY_DIR}/../cmake/modules)
ELSE(WIN32)
  set(CMAKE_MODULE_PATH ${PROJECT_SOURCE_DIR}/cmake/modules)
  SET(BOOST_ROOT $ENV{HOME}/boost)
  SET(Boost_USE_MULTITHREADED OFF)
ENDIF(WIN32)
find_package(Boost 1.47.0 REQUIRED COMPONENTS regex serialization filesystem system program_options)
IF(Boost_FOUND)
    MESSAGE("Boost includes located in: ${Boost_INCLUDE_DIRS}")
    MESSAGE("Boost libraries located in: ${Boost_LIBRARIES}")
ENDIF(Boost_FOUND)

SET(LIBRARY_OUTPUT_PATH ${PROJECT_BINARY_DIR}/lib)
SET(EXECUTABLE_OUTPUT_PATH ${PROJECT_BINARY_DIR}/bin)

# find and setup openbabel
find_package(PkgConfig REQUIRED)
find_package(OpenBabel2 REQUIRED COMPONENTS forcefields)
include_directories(${OPENBABEL2_INCLUDE_DIR})
IF(OPENBABEL2_FOUND)
    MESSAGE("Openbabel includes located in: ${OPENBABEL2_INCLUDE_DIR}")
    MESSAGE("Openbabel libraries located in: ${OPENBABEL2_LIBRARIES}")
ENDIF(OPENBABEL2_FOUND)

enable_testing()

# Turn off line breaks in error messages so that IDEs can parse the errors properly
IF(CMAKE_COMPILER_IS_GNUCC)
  SET(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -fmessage-length=0")  
ENDIF(CMAKE_COMPILER_IS_GNUCC)
IF(CMAKE_COMPILER_IS_GNUCXX)
  SET(CMAKE_CXX_FLAGS_DEBUG "${CMAKE_CXX_FLAGS_DEBUG} -O0 -Wall -Wextra -Weffc++ -Wno-format-errors")

  # Add -pg to the compile and link flags to enable profiling
  SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++11" )
  SET(CMAKE_CXX_FLAGS_DEBUG "${CMAKE_CXX_FLAGS_DEBUG} -O0 -Wall")
  SET(CMAKE_LINK_FLAGS_DEBUG "${CMAKE_LINK_FLAGS_DEBUG}")
  SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -fmessage-length=0")
ENDIF(CMAKE_COMPILER_IS_GNUCXX)

# Add support for OpenMP, if it is available.  Comment out this block if your installation
# of OpenMP gives you problems.
INCLUDE(FindOpenMP)
IF(OPENMP_FOUND)
  SET(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${OpenMP_C_FLAGS}")
  SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")
ENDIF(OPENMP_FOUND)

######### The camflow application ##########################
add_executable(camflow-app ${MOPSSUITE_SOURCE_DIR}/applications/solvers/camflow/cam_kernel.cpp)
target_link_libraries(camflow-app camflow ${Boost_LIBRARIES})

add_executable(camflowGridGenerator ${MOPSSUITE_SOURCE_DIR}/applications/utilities/camflowGridGenerator/camflowGridGenerator.cpp)
target_link_libraries(camflowGridGenerator ${Boost_LIBRARIES})

# Regression tests for CamFlow
add_test(camflow.hydrogenBatchReactor ${SHELL_INTERPRETER} ${MOPSSUITE_SOURCE_DIR}/test/camflow/batchReactorRegress.sh ${EXECUTABLE_OUTPUT_PATH}${MSVC_BUILD_DIR}/camflow-app${EXE_SUFFIX} ${MOPSSUITE_SOURCE_DIR}/test/camflow/hydrogenBatchReactor)
add_test(camflow.hydrogenFlamelet ${SHELL_INTERPRETER} ${MOPSSUITE_SOURCE_DIR}/test/camflow/hydrogenFlamelet.sh ${EXECUTABLE_OUTPUT_PATH}${MSVC_BUILD_DIR}/camflow-app${EXE_SUFFIX} ${MOPSSUITE_SOURCE_DIR}/test/camflow/hydrogenFlamelet)


######### The Mops application ##########################
add_executable(mops-app ${MOPSSUITE_SOURCE_DIR}/applications/solvers/mopsc/mops.cpp)
target_link_libraries(mops-app mops ${Boost_LIBRARIES})

# Regression tests for Mops
add_test(mops.regress1 ${SHELL_INTERPRETER} ${MOPSSUITE_SOURCE_DIR}/test/mopsc/regress1.sh ${EXECUTABLE_OUTPUT_PATH}${MSVC_BUILD_DIR}/mops-app${EXE_SUFFIX} ${MOPSSUITE_SOURCE_DIR}/test/mopsc)

add_test(mops.regress2 ${SHELL_INTERPRETER} ${MOPSSUITE_SOURCE_DIR}/test/mopsc/regress2.sh ${EXECUTABLE_OUTPUT_PATH}${MSVC_BUILD_DIR}/mops-app${EXE_SUFFIX} ${MOPSSUITE_SOURCE_DIR}/test/mopsc)

add_test(mops.regress3 ${SHELL_INTERPRETER} ${MOPSSUITE_SOURCE_DIR}/test/mopsc/regress3.sh ${EXECUTABLE_OUTPUT_PATH}${MSVC_BUILD_DIR}/mops-app${EXE_SUFFIX} ${MOPSSUITE_SOURCE_DIR}/test/mopsc)

add_test(mops.regress4 ${SHELL_INTERPRETER} ${MOPSSUITE_SOURCE_DIR}/test/mopsc/regress4.sh ${EXECUTABLE_OUTPUT_PATH}${MSVC_BUILD_DIR}/mops-app${EXE_SUFFIX} ${MOPSSUITE_SOURCE_DIR}/test/mopsc)

add_test(mops.regress5 ${SHELL_INTERPRETER} ${MOPSSUITE_SOURCE_DIR}/test/mopsc/regress5.sh ${EXECUTABLE_OUTPUT_PATH}${MSVC_BUILD_DIR}/mops-app${EXE_SUFFIX} ${MOPSSUITE_SOURCE_DIR}/test/mopsc)

add_test(mops.regress6 ${SHELL_INTERPRETER} ${MOPSSUITE_SOURCE_DIR}/test/mopsc/regress6.sh ${EXECUTABLE_OUTPUT_PATH}${MSVC_BUILD_DIR}/mops-app${EXE_SUFFIX} ${MOPSSUITE_SOURCE_DIR}/test/mopsc)

add_test(mops.regress7 ${SHELL_INTERPRETER} ${MOPSSUITE_SOURCE_DIR}/test/mopsc/regress7.sh ${EXECUTABLE_OUTPUT_PATH}${MSVC_BUILD_DIR}/mops-app${EXE_SUFFIX} ${MOPSSUITE_SOURCE_DIR}/test/mopsc)

add_test(mops.regress8 ${SHELL_INTERPRETER} ${MOPSSUITE_SOURCE_DIR}/test/mopsc/regress8.sh ${EXECUTABLE_OUTPUT_PATH}${MSVC_BUILD_DIR}/mops-app${EXE_SUFFIX} ${MOPSSUITE_SOURCE_DIR}/test/mopsc/regress8)

add_test(mops.psr1 ${SHELL_INTERPRETER} ${MOPSSUITE_SOURCE_DIR}/test/mopsc/psr1.sh ${EXECUTABLE_OUTPUT_PATH}${MSVC_BUILD_DIR}/mops-app${EXE_SUFFIX} ${MOPSSUITE_SOURCE_DIR}/test/mopsc)

add_test(mops.psr2 ${SHELL_INTERPRETER} ${MOPSSUITE_SOURCE_DIR}/test/mopsc/psr2.sh ${EXECUTABLE_OUTPUT_PATH}${MSVC_BUILD_DIR}/mops-app${EXE_SUFFIX} ${MOPSSUITE_SOURCE_DIR}/test/mopsc/psr2)

add_test(mops.sinter1 ${SHELL_INTERPRETER} ${MOPSSUITE_SOURCE_DIR}/test/mopsc/sinter1.sh ${EXECUTABLE_OUTPUT_PATH}${MSVC_BUILD_DIR}/mops-app${EXE_SUFFIX} ${MOPSSUITE_SOURCE_DIR}/test/mopsc)

add_test(mops.pahtest1 ${SHELL_INTERPRETER} ${MOPSSUITE_SOURCE_DIR}/test/mopsc/pahtest1.sh ${EXECUTABLE_OUTPUT_PATH}${MSVC_BUILD_DIR}/mops-app${EXE_SUFFIX} ${MOPSSUITE_SOURCE_DIR}/test/mopsc)

add_test(mops.pahtest2 ${SHELL_INTERPRETER} ${MOPSSUITE_SOURCE_DIR}/test/mopsc/pahtest2.sh ${EXECUTABLE_OUTPUT_PATH}${MSVC_BUILD_DIR}/mops-app${EXE_SUFFIX} ${MOPSSUITE_SOURCE_DIR}/test/mopsc)

add_test(mops.pahtest3 ${SHELL_INTERPRETER} ${MOPSSUITE_SOURCE_DIR}/test/mopsc/pahtest3.sh ${EXECUTABLE_OUTPUT_PATH}${MSVC_BUILD_DIR}/mops-app${EXE_SUFFIX} ${MOPSSUITE_SOURCE_DIR}/test/mopsc)

add_test(mops.silica1 ${SHELL_INTERPRETER} ${MOPSSUITE_SOURCE_DIR}/test/mopsc/silica1.sh ${EXECUTABLE_OUTPUT_PATH}${MSVC_BUILD_DIR}/mops-app${EXE_SUFFIX} ${MOPSSUITE_SOURCE_DIR}/test/mopsc/silica1)

add_test(mops.silica2 ${SHELL_INTERPRETER} ${MOPSSUITE_SOURCE_DIR}/test/mopsc/silica2.sh ${EXECUTABLE_OUTPUT_PATH}${MSVC_BUILD_DIR}/mops-app${EXE_SUFFIX} ${MOPSSUITE_SOURCE_DIR}/test/mopsc/silica2)

add_test(mops.bintree1 ${SHELL_INTERPRETER} ${MOPSSUITE_SOURCE_DIR}/test/mopsc/bintree1.sh ${EXECUTABLE_OUTPUT_PATH}${MSVC_BUILD_DIR}/mops-app${EXE_SUFFIX} ${MOPSSUITE_SOURCE_DIR}/test/mopsc)

add_test(mops.bintree2a ${SHELL_INTERPRETER} ${MOPSSUITE_SOURCE_DIR}/test/mopsc/bintree2a.sh ${EXECUTABLE_OUTPUT_PATH}${MSVC_BUILD_DIR}/mops-app${EXE_SUFFIX} ${MOPSSUITE_SOURCE_DIR}/test/mopsc/bintree2)

add_test(mops.bintree2b ${SHELL_INTERPRETER} ${MOPSSUITE_SOURCE_DIR}/test/mopsc/bintree2b.sh ${EXECUTABLE_OUTPUT_PATH}${MSVC_BUILD_DIR}/mops-app${EXE_SUFFIX} ${MOPSSUITE_SOURCE_DIR}/test/mopsc/bintree2)

add_test(mops.bintree2c ${SHELL_INTERPRETER} ${MOPSSUITE_SOURCE_DIR}/test/mopsc/bintree2c.sh ${EXECUTABLE_OUTPUT_PATH}${MSVC_BUILD_DIR}/mops-app${EXE_SUFFIX} ${MOPSSUITE_SOURCE_DIR}/test/mopsc/bintree2)

add_test(mops.imaging1 ${SHELL_INTERPRETER} ${MOPSSUITE_SOURCE_DIR}/test/mopsc/imaging1.sh ${EXECUTABLE_OUTPUT_PATH}${MSVC_BUILD_DIR}/mops-app${EXE_SUFFIX} ${MOPSSUITE_SOURCE_DIR}/test/mopsc/imaging1)

add_test(mops.titania1 ${SHELL_INTERPRETER} ${MOPSSUITE_SOURCE_DIR}/test/mopsc/titania1.sh ${EXECUTABLE_OUTPUT_PATH}${MSVC_BUILD_DIR}/mops-app${EXE_SUFFIX} ${MOPSSUITE_SOURCE_DIR}/test/mopsc/titania1)

add_test(mops.titania3 ${SHELL_INTERPRETER} ${MOPSSUITE_SOURCE_DIR}/test/mopsc/titania3.sh ${EXECUTABLE_OUTPUT_PATH}${MSVC_BUILD_DIR}/mops-app${EXE_SUFFIX} ${MOPSSUITE_SOURCE_DIR}/test/mopsc/titania3)

add_test(mops.stagnation1 ${SHELL_INTERPRETER} ${MOPSSUITE_SOURCE_DIR}/test/mopsc/stagnation1.sh ${EXECUTABLE_OUTPUT_PATH}${MSVC_BUILD_DIR}/mops-app${EXE_SUFFIX} ${MOPSSUITE_SOURCE_DIR}/test/mopsc/stagnation1)

add_test(mops.stagnation2 ${SHELL_INTERPRETER} ${MOPSSUITE_SOURCE_DIR}/test/mopsc/stagnation2.sh ${EXECUTABLE_OUTPUT_PATH}${MSVC_BUILD_DIR}/mops-app${EXE_SUFFIX} ${MOPSSUITE_SOURCE_DIR}/test/mopsc/stagnation2)

add_test(mops.titaniaphase1 ${SHELL_INTERPRETER} ${MOPSSUITE_SOURCE_DIR}/test/mopsc/titaniaphase1.sh ${EXECUTABLE_OUTPUT_PATH}${MSVC_BUILD_DIR}/mops-app${EXE_SUFFIX} ${MOPSSUITE_SOURCE_DIR}/test/mopsc/titaniaphase1)

add_test(mops.titaniaphase2 ${SHELL_INTERPRETER} ${MOPSSUITE_SOURCE_DIR}/test/mopsc/titaniaphase2.sh ${EXECUTABLE_OUTPUT_PATH}${MSVC_BUILD_DIR}/mops-app${EXE_SUFFIX} ${MOPSSUITE_SOURCE_DIR}/test/mopsc/titaniaphase2)

add_test(mops.titaniaphase3 ${SHELL_INTERPRETER} ${MOPSSUITE_SOURCE_DIR}/test/mopsc/titaniaphase3.sh ${EXECUTABLE_OUTPUT_PATH}${MSVC_BUILD_DIR}/mops-app${EXE_SUFFIX} ${MOPSSUITE_SOURCE_DIR}/test/mopsc/titaniaphase3)

add_test(mops.titaniahybrid1 ${SHELL_INTERPRETER} ${MOPSSUITE_SOURCE_DIR}/test/mopsc/titaniahybrid1.sh ${EXECUTABLE_OUTPUT_PATH}${MSVC_BUILD_DIR}/mops-app${EXE_SUFFIX} ${MOPSSUITE_SOURCE_DIR}/test/mopsc/titaniahybrid1)

# Networking interface tests
add_test(mops.network1a ${SHELL_INTERPRETER} ${MOPSSUITE_SOURCE_DIR}/test/mopsc/network1.sh ${EXECUTABLE_OUTPUT_PATH}${MSVC_BUILD_DIR}/mops-app${EXE_SUFFIX} ${MOPSSUITE_SOURCE_DIR}/test/mopsc/network1 "a")

add_test(mops.network1b ${SHELL_INTERPRETER} ${MOPSSUITE_SOURCE_DIR}/test/mopsc/network1.sh ${EXECUTABLE_OUTPUT_PATH}${MSVC_BUILD_DIR}/mops-app${EXE_SUFFIX} ${MOPSSUITE_SOURCE_DIR}/test/mopsc/network1 "b")

add_test(mops.network1c ${SHELL_INTERPRETER} ${MOPSSUITE_SOURCE_DIR}/test/mopsc/network1.sh ${EXECUTABLE_OUTPUT_PATH}${MSVC_BUILD_DIR}/mops-app${EXE_SUFFIX} ${MOPSSUITE_SOURCE_DIR}/test/mopsc/network1 "c")

add_test(mops.network1d ${SHELL_INTERPRETER} ${MOPSSUITE_SOURCE_DIR}/test/mopsc/network1.sh ${EXECUTABLE_OUTPUT_PATH}${MSVC_BUILD_DIR}/mops-app${EXE_SUFFIX} ${MOPSSUITE_SOURCE_DIR}/test/mopsc/network1 "d")

add_test(mops.network2 ${SHELL_INTERPRETER} ${MOPSSUITE_SOURCE_DIR}/test/mopsc/network2.sh ${EXECUTABLE_OUTPUT_PATH}${MSVC_BUILD_DIR}/mops-app${EXE_SUFFIX} ${MOPSSUITE_SOURCE_DIR}/test/mopsc/network2)

########## The PAH-KMC Application ######################
add_executable(PAHkmc-app ${MOPSSUITE_SOURCE_DIR}/applications/solvers/PAHkmc/kmc_model.cpp)
target_link_libraries(PAHkmc-app sweep ${Boost_LIBRARIES} ${OPENBABEL2_LIBRARIES})
#target_link_libraries(PAHkmc-app sweep /home/gl413/Codes/Openbabel/install/lib/libopenbabel.so)

########## The Brush application #######################
add_executable(brush-app ${MOPSSUITE_SOURCE_DIR}/applications/solvers/brush/brush.cpp)
target_link_libraries(brush-app brush ${Boost_LIBRARIES})

# Regression tests for Brush
add_test(brush.regress1 ${SHELL_INTERPRETER} ${MOPSSUITE_SOURCE_DIR}/test/brush/regress1.sh ${EXECUTABLE_OUTPUT_PATH}${MSVC_BUILD_DIR}/brush-app${EXE_SUFFIX} ${MOPSSUITE_SOURCE_DIR}/test/brush)

add_test(brush.regress3 ${SHELL_INTERPRETER} ${MOPSSUITE_SOURCE_DIR}/test/brush/regress3.sh ${EXECUTABLE_OUTPUT_PATH}${MSVC_BUILD_DIR}/brush-app${EXE_SUFFIX} ${MOPSSUITE_SOURCE_DIR}/test/brush)

add_test(brush.regress4 ${SHELL_INTERPRETER} ${MOPSSUITE_SOURCE_DIR}/test/brush/regress4.sh ${EXECUTABLE_OUTPUT_PATH}${MSVC_BUILD_DIR}/brush-app${EXE_SUFFIX} ${MOPSSUITE_SOURCE_DIR}/test/brush)

add_test(brush.regress5 ${SHELL_INTERPRETER} ${MOPSSUITE_SOURCE_DIR}/test/brush/regress5.sh ${EXECUTABLE_OUTPUT_PATH}${MSVC_BUILD_DIR}/brush-app${EXE_SUFFIX} ${MOPSSUITE_SOURCE_DIR}/test/brush)

add_test(brush.regress6 ${SHELL_INTERPRETER} ${MOPSSUITE_SOURCE_DIR}/test/brush/regress6.sh ${EXECUTABLE_OUTPUT_PATH}${MSVC_BUILD_DIR}/brush-app${EXE_SUFFIX} ${MOPSSUITE_SOURCE_DIR}/test/brush)
########## A small program to exercise the Brush interface for moonmd #######################
add_executable(brush-moonmd-test ${MOPSSUITE_SOURCE_DIR}/applications/test-harnesses/brush/moonmd_test.cpp)
target_link_libraries(brush-moonmd-test brush ${Boost_LIBRARIES})

# Regression tests for Brush moon-md interface
add_test(brush.moonmd1 ${SHELL_INTERPRETER} ${MOPSSUITE_SOURCE_DIR}/test/brush/moonmd1.sh ${EXECUTABLE_OUTPUT_PATH}${MSVC_BUILD_DIR}/brush-moonmd-test${EXE_SUFFIX} ${MOPSSUITE_SOURCE_DIR}/test/brush)


########## Test Program for sprogc ######################
add_executable(sprogc-test ${MOPSSUITE_SOURCE_DIR}/applications/test-harnesses/sprogc/sprogc_test.cpp)
target_link_libraries(sprogc-test sprog ${Boost_LIBRARIES})

add_test(sprogc.regress1 ${SHELL_INTERPRETER} ${MOPSSUITE_SOURCE_DIR}/test/sprogc/regress1.sh ${EXECUTABLE_OUTPUT_PATH}${MSVC_BUILD_DIR}/sprogc-test${EXE_SUFFIX} ${MOPSSUITE_SOURCE_DIR}/test/sprogc)

add_test(sprogc.regress2 ${SHELL_INTERPRETER} ${MOPSSUITE_SOURCE_DIR}/test/sprogc/regress2.sh ${EXECUTABLE_OUTPUT_PATH}${MSVC_BUILD_DIR}/mops-app${EXE_SUFFIX} ${MOPSSUITE_SOURCE_DIR}/test/sprogc/regress2)

########## Test Program for chemkinReader ######################
add_executable(chemkinReader-test ${MOPSSUITE_SOURCE_DIR}/applications/test-harnesses/chemkinReader/chemkinReaderTest.cpp)
target_link_libraries(chemkinReader-test chemkinReader ${Boost_LIBRARIES})

add_test(chemkinReader.regress1 ${SHELL_INTERPRETER} ${MOPSSUITE_SOURCE_DIR}/test/io/chemkinReaderRegress.sh ${EXECUTABLE_OUTPUT_PATH}${MSVC_BUILD_DIR}/chemkinReader-test${EXE_SUFFIX} ${MOPSSUITE_SOURCE_DIR}/test/io)

########## Test Program for linear interpolator ######################
add_executable(linInterp-test ${MOPSSUITE_SOURCE_DIR}/applications/test-harnesses/utils/test_linear_interpolator.cpp)

add_test(NAME utils.lininterp1 COMMAND linInterp-test)

########## Test program for Sweep::FixedChemistry ##############
add_executable(sweepFixedMix-test ${MOPSSUITE_SOURCE_DIR}/applications/test-harnesses/sweepc/fixedmix_test.cpp)
target_link_libraries(sweepFixedMix-test brush ${Boost_LIBRARIES})
add_test(NAME sweep.FixedMix COMMAND ${SHELL_INTERPRETER} ${MOPSSUITE_SOURCE_DIR}/test/sweepc/fixedmix.sh $<TARGET_FILE:sweepFixedMix-test> ${MOPSSUITE_SOURCE_DIR}/test/sweepc/fixedmix)

# Subsidiary libraries for the solvers
add_subdirectory(${MOPSSUITE_SOURCE_DIR}/src/io/chemkinReader)
add_subdirectory(${MOPSSUITE_SOURCE_DIR}/src/io/comostrings)
add_subdirectory(${MOPSSUITE_SOURCE_DIR}/src/io/camxml)
add_subdirectory(${MOPSSUITE_SOURCE_DIR}/src/geometry)
add_subdirectory(${MOPSSUITE_SOURCE_DIR}/src/sprogc)
add_subdirectory(${MOPSSUITE_SOURCE_DIR}/src/sweepc)
add_subdirectory(${MOPSSUITE_SOURCE_DIR}/src/odesolvers/cvodes)
add_subdirectory(${MOPSSUITE_SOURCE_DIR}/src/odesolvers/cvode)
add_subdirectory(${MOPSSUITE_SOURCE_DIR}/src/odesolvers/radau)
add_subdirectory(${MOPSSUITE_SOURCE_DIR}/src/camflow)
add_subdirectory(${MOPSSUITE_SOURCE_DIR}/src/mopsc)
add_subdirectory(${MOPSSUITE_SOURCE_DIR}/src/brush)

# All the header files we have written
include_directories(${Boost_INCLUDE_DIRS}
                    ${MOPSSUITE_SOURCE_DIR}/src/utils/include
                    ${MOPSSUITE_SOURCE_DIR}/src/io/chemkinReader/include
                    ${MOPSSUITE_SOURCE_DIR}/src/io/comostrings/include
                    ${MOPSSUITE_SOURCE_DIR}/src/io/camxml/include
                    ${MOPSSUITE_SOURCE_DIR}/src/geometry/include
                    ${MOPSSUITE_SOURCE_DIR}/src/sprogc/include
                    ${MOPSSUITE_SOURCE_DIR}/src/sweepc/include
                    ${MOPSSUITE_SOURCE_DIR}/src/odesolvers/cvodes/include
                    ${MOPSSUITE_SOURCE_DIR}/src/odesolvers/cvodes/source/cvodes
                    ${MOPSSUITE_SOURCE_DIR}/src/odesolvers/cvode/include
                    ${MOPSSUITE_SOURCE_DIR}/src/odesolvers/radau/include
                    ${MOPSSUITE_SOURCE_DIR}/src/camflow/include
                    ${MOPSSUITE_SOURCE_DIR}/src/mopsc/include
                    ${MOPSSUITE_SOURCE_DIR}/src/brush/include)
					#/home/gl413/Codes/Openbabel/install/include/openbabel-2.0)
install(TARGETS mops-app DESTINATION ${PROJECT_SOURCE_DIR}/bin/try)