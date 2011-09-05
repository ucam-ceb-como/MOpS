# Specification file for the CMake build system
#
# Robert I A Patterson
# 29 April 2011
cmake_minimum_required(VERSION 2.6)

project (SPROG)

file(GLOB_RECURSE INCS "include/*.h")
source_group("Header Files" FILES ${INCS})

add_library(sprog source/gpc_el_comp.cpp
                  source/gpc_element.cpp
                  source/gpc_gasphase.cpp
                  source/gpc_idealgas.cpp
                  source/gpc_mech.cpp
                  source/gpc_mech_io.cpp
                  source/gpc_mixture.cpp
                  source/gpc_mixture_factory.cpp
                  source/gpc_reaction.cpp
                  source/gpc_reaction_set.cpp
                  source/gpc_species.cpp
                  source/gpc_stoich.cpp
                  source/gpc_thermo.cpp
                  source/gpc_transport_factory.cpp
                  source/gpc_element_lib.cpp
                  ${INCS})

include_directories(include 
                    ${Boost_INCLUDE_DIRS}
                    ../utils/include
                    ../io/chemkinReader/include
                    ../io/comostrings/include
                    ../io/camxml/include)

# Linking against sprog will also require the following libraries
# to resolve references from within sprog.
target_link_libraries(sprog camxml comostrings chemkinReader)
