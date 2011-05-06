# The name of this project is "MOPSSUITE". CMakeLists files in this project can
# refer to the root source directory of the project as ${MOPSSUITE_SOURCE_DIR} and
# to the root binary directory of the project as ${MOPSSUITE_BINARY_DIR}.
cmake_minimum_required (VERSION 2.6)
project (MOPSSUITE)

enable_testing()

# Mops and the libraries it requires, the order in which the link libraries is listed
# is sensitive: the most basic libraries most come last.
add_executable(mops-app ${MOPSSUITE_SOURCE_DIR}/applications/solvers/mopsc/mops.cpp)
target_link_libraries(mops-app mops)

add_test(mops.runs mops-app)

# Brush and the libraries it requires
add_executable(brush-app ${MOPSSUITE_SOURCE_DIR}/applications/solvers/brush/brush.cpp)
target_link_libraries(brush-app brush)

add_test(brush.runs brush-app -h)
add_test(brush.regress1 ${MOPSSUITE_SOURCE_DIR}/test/brush/regress1.sh ${MOPSSUITE_BINARY_DIR}/brush-app ${MOPSSUITE_SOURCE_DIR}/test/brush)
add_test(brush.regress2 ${MOPSSUITE_SOURCE_DIR}/test/brush/regress2.sh ${MOPSSUITE_BINARY_DIR}/brush-app ${MOPSSUITE_SOURCE_DIR}/test/brush)
add_test(brush.regress3 ${MOPSSUITE_SOURCE_DIR}/test/brush/regress3.sh ${MOPSSUITE_BINARY_DIR}/brush-app ${MOPSSUITE_SOURCE_DIR}/test/brush)
add_test(brush.regress4 ${MOPSSUITE_SOURCE_DIR}/test/brush/regress4.sh ${MOPSSUITE_BINARY_DIR}/brush-app ${MOPSSUITE_SOURCE_DIR}/test/brush)

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
