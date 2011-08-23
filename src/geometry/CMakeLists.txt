# Specification file for the CMake build system
#
# Robert I A Patterson
# 29 April 2011
cmake_minimum_required(VERSION 2.6)

project (GEOMETRY)

file(GLOB_RECURSE INCS "include/*.h")
source_group("Header Files" FILES ${INCS})

add_library(geometry source/geometry1d.cpp
                     source/local_geometry1d.cpp
				     ${INCS})

include_directories(include 
                    ../io/camxml/include
                    ../sprogc/include)

# Linking against geometry will require linking against
# the following libraries to resolve links from within
# geometry.
target_link_libraries(geometry sprog camxml)
