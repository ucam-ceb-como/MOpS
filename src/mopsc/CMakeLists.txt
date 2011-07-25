# Specification file for the CMake build system
#
# Robert I A Patterson
# 29 April 2011
cmake_minimum_required(VERSION 2.6)

project (MOPS)

add_library(mops source/cvodes_utils.cpp
                 source/loi_reduction.cpp
                 source/mops_flow_stream.cpp
                 source/mops_flux_postprocessor.cpp
                 source/mops_gpc_sensitivity.cpp
                 source/mops_mechanism.cpp
                 source/mops_mixture.cpp
                 source/mops_ode_solver.cpp
                 source/mops_particle_solver.cpp
                 source/mops_predcor_solver.cpp
                 source/mops_psr.cpp
                 source/mops_reactor.cpp
                 source/mops_reactor_factory.cpp
                 source/mops_rhs_func.cpp
                 source/mops_settings_io.cpp
                 source/mops_simplesplit_solver.cpp
                 source/mops_simulator.cpp
                 source/mops_solver.cpp
                 source/mops_src_terms.cpp
                 source/mops_strang_solver.cpp
                 source/mops_timeinterval.cpp
                 source/swp_flamesolver.cpp
                 source/swp_gas_profile.cpp)

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
                    ../odesolvers/cvodes/source/cvodes)

target_link_libraries(mops sweep sprog cvodes camxml comostrings chemkinReader)