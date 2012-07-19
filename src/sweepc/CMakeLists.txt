# Specification file for the CMake build system
#
# Robert I A Patterson
# 29 April 2011
cmake_minimum_required(VERSION 2.6)

project (SWEEP)

file(GLOB_RECURSE INCS "include/*.h")
source_group("Header Files" FILES ${INCS})

add_library(sweep source/swp_actsites_reaction.cpp
                  source/swp_addcoag.cpp
                  source/swp_sprog_idealgas_wrapper.cpp
                  source/swp_bintree_primary.cpp
                  source/swp_bintree_stats.cpp
                  source/swp_birth_process.cpp
                  source/swp_cell.cpp
                  source/swp_coagulation.cpp
                  source/swp_component.cpp
                  source/swp_condensation.cpp
                  source/swp_constant_inception.cpp
                  source/swp_constcoag.cpp
                  source/swp_death_process.cpp
                  source/swp_dimer_inception.cpp
                  source/swp_ensemble.cpp
                  source/swp_ensemble_stats.cpp
                  source/swp_gas_profile.cpp
                  source/swp_imgnode.cpp
                  source/swp_inception.cpp
                  source/swp_mechanism.cpp
                  source/swp_mech_parser.cpp
                  source/swp_model_factory.cpp
                  source/swp_PAH.cpp
                  source/swp_pah_inception.cpp
                  source/swp_PAH_primary.cpp
                  source/swp_PAH_stats.cpp
                  source/swp_particle.cpp
                  source/swp_particle_image.cpp
                  source/swp_particle_model.cpp
                  source/swp_particle_process.cpp
                  source/swp_particle_stats.cpp
                  source/swp_primary.cpp
                  source/swp_process.cpp
                  source/swp_process_factory.cpp
                  source/swp_silica_interparticle.cpp
                  source/swp_silica_primary.cpp
                  source/swp_silica_stats.cpp
                  source/swp_silicon_inception.cpp
                  source/swp_sintering_model.cpp
                  source/swp_solver.cpp
                  source/swp_subparticle.cpp
                  source/swp_surface_reaction.cpp
                  source/swp_surfvol_primary.cpp
                  source/swp_surfvolhydrogen_primary.cpp
                  source/swp_surfvol_stats.cpp
                  source/swp_surfvolhydrogen_stats.cpp
                  source/swp_tracker.cpp
                  source/swp_transcoag.cpp
                  source/swp_tree_add_weighted_cache.cpp
                  source/swp_tree_cache.cpp
                  source/swp_tree_transcoag_weighted_cache.cpp
                  source/swp_tree_weighted_cache.cpp
                  source/swp_weighted_addcoag.cpp
                  source/swp_weighted_constcoag.cpp
                  source/swp_weighted_transcoag.cpp
                  source/swp_kmc_gaspoint.cpp
                  source/swp_kmc_jump_process.cpp
                  source/swp_kmc_mech.cpp
                  source/swp_kmc_pah_process.cpp
                  source/swp_kmc_pah_structure.cpp
                  source/swp_kmc_reaction.cpp
                  source/swp_kmc_simulator.cpp
                  source/swp_kmc_structure_comp.cpp
                  ${INCS})

include_directories(include 
                    ${Boost_INCLUDE_DIRS} 
                    ../utils/include
                    ../io/chemkinReader/include
                    ../io/comostrings/include
                    ../io/camxml/include
                    ../geometry/include
                    ../sprogc/include)

# To link againast sweep the following libraries will also
# be required to resolve references from within sweep.
target_link_libraries(sweep sprog geometry camxml comostrings chemkinReader)
