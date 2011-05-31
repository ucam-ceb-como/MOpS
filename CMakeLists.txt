# The name of this project is "MOPSSUITE". CMakeLists files in this project can
# refer to the root source directory of the project as ${MOPSSUITE_SOURCE_DIR} and
# to the root binary directory of the project as ${MOPSSUITE_BINARY_DIR}.
cmake_minimum_required (VERSION 2.6)
project (MOPSSUITE)

enable_testing()

# Turn off line breaks in error messages so that IDEs can parse the errors properly
IF(CMAKE_COMPILER_IS_GNUCC)
  SET(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -fmessage-length=0")
ENDIF(CMAKE_COMPILER_IS_GNUCC)
IF(CMAKE_COMPILER_IS_GNUCXX)
  SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -fmessage-length=0")
ENDIF(CMAKE_COMPILER_IS_GNUCXX)


######### The Mops application ##########################
add_executable(mops-app ${MOPSSUITE_SOURCE_DIR}/applications/solvers/mopsc/mops.cpp)
target_link_libraries(mops-app mops)

# Regression tests for Mops
add_test(mops.regress1 ${MOPSSUITE_SOURCE_DIR}/test/mopsc/regress1.sh ${MOPSSUITE_BINARY_DIR}/mops-app ${MOPSSUITE_SOURCE_DIR}/test/mopsc)

add_test(mops.regress2 ${MOPSSUITE_SOURCE_DIR}/test/mopsc/regress2.sh ${MOPSSUITE_BINARY_DIR}/mops-app ${MOPSSUITE_SOURCE_DIR}/test/mopsc)

add_test(mops.regress3 ${MOPSSUITE_SOURCE_DIR}/test/mopsc/regress3.sh ${MOPSSUITE_BINARY_DIR}/mops-app ${MOPSSUITE_SOURCE_DIR}/test/mopsc)

add_test(mops.regress4 ${MOPSSUITE_SOURCE_DIR}/test/mopsc/regress4.sh ${MOPSSUITE_BINARY_DIR}/mops-app ${MOPSSUITE_SOURCE_DIR}/test/mopsc)

add_test(mops.CSTRtest ${MOPSSUITE_SOURCE_DIR}/test/mopsc/CSTRRegress.sh ${MOPSSUITE_BINARY_DIR}/mops-app ${MOPSSUITE_SOURCE_DIR}/test/mopsc)

add_test(mops.pahtest ${MOPSSUITE_SOURCE_DIR}/test/mopsc/pahtest.sh ${MOPSSUITE_BINARY_DIR}/mops-app ${MOPSSUITE_SOURCE_DIR}/test/mopsc)

########## The Brush application #######################
add_executable(brush-app ${MOPSSUITE_SOURCE_DIR}/applications/solvers/brush/brush.cpp)
target_link_libraries(brush-app brush)

# Regression tests for Brush
add_test(brush.regress1 ${MOPSSUITE_SOURCE_DIR}/test/brush/regress1.sh ${MOPSSUITE_BINARY_DIR}/brush-app ${MOPSSUITE_SOURCE_DIR}/test/brush)

add_test(brush.regress2 ${MOPSSUITE_SOURCE_DIR}/test/brush/regress2.sh ${MOPSSUITE_BINARY_DIR}/brush-app ${MOPSSUITE_SOURCE_DIR}/test/brush)

add_test(brush.regress3 ${MOPSSUITE_SOURCE_DIR}/test/brush/regress3.sh ${MOPSSUITE_BINARY_DIR}/brush-app ${MOPSSUITE_SOURCE_DIR}/test/brush)

add_test(brush.regress4 ${MOPSSUITE_SOURCE_DIR}/test/brush/regress4.sh ${MOPSSUITE_BINARY_DIR}/brush-app ${MOPSSUITE_SOURCE_DIR}/test/brush)

########## A small program to exercise the Brush interface for moonmd #######################
add_executable(brush-moonmd-test ${MOPSSUITE_SOURCE_DIR}/applications/test-harnesses/brush/moonmd_test.cpp)
target_link_libraries(brush-moonmd-test brush)

# Regression tests for Brush moon-md interface
add_test(NAME brush.moonmd1 COMMAND ${MOPSSUITE_SOURCE_DIR}/test/brush/moonmd1.sh ${MOPSSUITE_BINARY_DIR}/brush-moonmd-test ${MOPSSUITE_SOURCE_DIR}/test/brush)

# Subsidiary libraries for the solvers
add_subdirectory(${MOPSSUITE_SOURCE_DIR}/src/io/comostrings)
add_subdirectory(${MOPSSUITE_SOURCE_DIR}/src/io/camxml)
add_subdirectory(${MOPSSUITE_SOURCE_DIR}/src/geometry)
add_subdirectory(${MOPSSUITE_SOURCE_DIR}/src/sprogc)
add_subdirectory(${MOPSSUITE_SOURCE_DIR}/src/sweepc)
add_subdirectory(${MOPSSUITE_SOURCE_DIR}/src/odesolvers/cvodes)
add_subdirectory(${MOPSSUITE_SOURCE_DIR}/src/mopsc)
add_subdirectory(${MOPSSUITE_SOURCE_DIR}/src/brush)

# All the header files we have written
include_directories(${MOPSSUITE_SOURCE_DIR}/src/utils/include
                    ${MOPSSUITE_SOURCE_DIR}/src/io/comostrings/include
                    ${MOPSSUITE_SOURCE_DIR}/src/io/camxml/include
                    ${MOPSSUITE_SOURCE_DIR}/src/geometry/include
                    ${MOPSSUITE_SOURCE_DIR}/src/sprogc/include
                    ${MOPSSUITE_SOURCE_DIR}/src/sweepc/include
                    ${MOPSSUITE_SOURCE_DIR}/src/odesolvers/cvodes/include
                    ${MOPSSUITE_SOURCE_DIR}/src/odesolvers/cvodes/source/cvodes
                    ${MOPSSUITE_SOURCE_DIR}/src/mopsc/include
                    ${MOPSSUITE_SOURCE_DIR}/src/brush/include)
