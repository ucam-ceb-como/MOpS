# Specification file for the CMake build system

#
# Robert I A Patterson
# 29 April 2011
cmake_minimum_required(VERSION 2.6)

project (RADAU)

file(GLOB_RECURSE INCS "include/*.h")
source_group("Header Files" FILES ${INCS})

add_library(radau source/decsol.cpp
                  source/IntegratorT.cpp
                  source/StiffIntegratorT.cpp
				  ${INCS})

include_directories(include)
