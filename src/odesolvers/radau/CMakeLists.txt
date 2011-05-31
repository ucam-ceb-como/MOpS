# Specification file for the CMake build system

#
# Robert I A Patterson
# 29 April 2011
cmake_minimum_required(VERSION 2.6)

project (RADAU)

add_library(radau source/decsol.cpp
                  source/IntegratorT.cpp
                  source/StiffIntegratorT.cpp)

include_directories(include)
