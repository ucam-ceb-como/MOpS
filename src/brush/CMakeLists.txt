# Specification file for the CMake build system
#
# Robert I A Patterson
# 29 April 2011
cmake_minimum_required(VERSION 2.6)

project (BRUSH)

add_library(brush source/pred_corr_solver.cpp
                  source/reactor1d.cpp
                  source/reset_chemistry.cpp
                  source/settings_io.cpp
                  source/simulator.cpp
                  source/moonmd_interface.cpp)

include_directories(include 
                    ${Boost_INCLUDE_DIRS}
                    ../utils/include
                    ../io/chemkinReader/include
                    ../io/comostrings/include
                    ../io/camxml/include
                    ../geometry/include
                    ../sprogc/include
                    ../sweepc/include
                    ../odesolvers/cvodes/include
                    ../odesolvers/cvodes/source/cvodes
                    ../mopsc/include)

# To link against brush the following additional libraries will be required
# to resolve references from within brush.
target_link_libraries(brush mops sweep sprog cvodes camxml comostrings chemkinReader)
