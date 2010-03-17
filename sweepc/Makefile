# Compiler options.
export debug=0
export profile=0

# C++ compiler.
cpp=g++
#cpp=mpicxx -DUSE_MPI -DNONSPHERE
ifeq ($(debug),0)
  ifeq ($(profile),0)
    flag=-c -O2 -Iinclude -I../comostrings/include -I../camxml/include \
       -I../sprogc/include -I../geometry/include -I../utils/include -Wall
  else
    flag=-c -O2 -Iinclude -I../comostrings/include -I../camxml/include \
       -I../sprogc/include -I../geometry/include -I../utils/include -Wall -pg
  endif
else
  flag=-c -O0 -Iinclude -I../comostrings/include -I../camxml/include \
       -I../sprogc/include -I../geometry/include -I../utils/include -Wall \
       -g -pg -D_GLIBCXX_DEBUG
endif

# Source directories.

export path=.
src=$(path)/source
inc=$(path)/include
lib=$(path)/lib

# Additional libraries.

ifeq ($(debug),0)
  comostrings=../comostrings/lib/comostrings.a
  camxml=../camxml/lib/camxml.a
  sprog=../sprogc/lib/sprog.a
  geom=../geom/lib/geom.a
else
  comostrings=../comostrings/lib/comostrings_d.a
  camxml=../camxml/lib/camxml_d.a
  sprog=../sprogc/lib/sprog_d.a
  geom=../geom/lib/geom_d.a
endif

# Target objects.

obj=swp_cell.o swp_component.o swp_ensemble.o swp_ensemble_stats.o \
    swp_mech_parser.o swp_mechanism.o swp_subparticle.o swp_particle.o \
    swp_particle_cache.o swp_particle_stats.o swp_solver.o swp_tracker.o \
    swp_treenode.o swp_abf_model.o swp_actsites_model.o swp_aggmodel_cache.o \
    swp_submodel.o swp_submodel_cache.o swp_model_factory.o \
    swp_surfvol_cache.o swp_surfvol_primary.o swp_primary.o \
    swp_pripart_cache.o swp_pripart_primary.o swp_actsites_reaction.o \
    swp_coagulation.o swp_transcoag.o swp_addcoag.o swp_condensation.o swp_inception.o \
    swp_particle_process.o swp_process.o swp_process_factory.o \
    swp_surface_reaction.o swp_particle_model.o swp_imgnode.o \
    swp_particle_image.o swp_ensemble_image.o swp_ensemble_imgnode.o swp_surfvol_stats.o swp_pripart_stats.o \
    swp_sintering_model.o swp_arssc_model.o swp_arssc_inception.o \
    swp_death_process.o swp_birth_process.o swp_arssc_cache.o \
    swp_arssc_reaction.o swp_arssc_condensation.o swp_arssc_process.o \
    swp_PAH_primary.o swp_PAH_cache.o swp_PAH_stats.o swp_transport_process.o \
    swp_diffusion_process.o swp_advection_process.o \
    mt19937.o swp_dimer_inception.o swp_pah_inception.o \
    swp_molecule_evolution.o

# TARGET all: Compile to a static library (*.a).

all: $(obj)
	mkdir -p $(lib)
ifeq ($(debug),0)
	ar rc $(lib)/sweep.a $(obj)
else
	ar rc $(lib)/sweep_d.a $(obj)
endif

# Source code.

# Particle definition.

swp_component.o: $(src)/swp_component.cpp
	$(cpp) $(flag) $(src)/swp_component.cpp

swp_tracker.o: $(src)/swp_tracker.cpp
	$(cpp) $(flag) $(src)/swp_tracker.cpp

swp_particle_cache.o: $(src)/swp_particle_cache.cpp $(inc)/swp_particle_cache.h
	$(cpp) $(flag) $(src)/swp_particle_cache.cpp

swp_subparticle.o: $(src)/swp_subparticle.cpp $(inc)/swp_subparticle.h
	$(cpp) $(flag) $(src)/swp_subparticle.cpp

swp_particle.o: $(src)/swp_particle.cpp
	$(cpp) $(flag) $(src)/swp_particle.cpp

# Particle system & ensemble.

swp_cell.o: $(src)/swp_cell.cpp $(inc)/swp_mechanism.h
	$(cpp) $(flag) $(src)/swp_cell.cpp

swp_treenode.o: $(src)/swp_treenode.cpp $(inc)/swp_treenode.h
	$(cpp) $(flag) $(src)/swp_treenode.cpp

swp_ensemble.o: $(src)/swp_ensemble.cpp $(inc)/swp_ensemble.h $(inc)/swp_ensemble_inl.h ../utils/$(inc)/binary_tree.hpp
	$(cpp) $(flag) $(src)/swp_ensemble.cpp

# Statistics generation.

swp_particle_stats.o: $(src)/swp_particle_stats.cpp  $(inc)/swp_particle.h
	$(cpp) $(flag) $(src)/swp_particle_stats.cpp
	
swp_ensemble_stats.o: $(src)/swp_ensemble_stats.cpp  $(inc)/swp_model_factory.h \
                       $(inc)/swp_model_factory.h
	$(cpp) $(flag) $(src)/swp_ensemble_stats.cpp

swp_surfvol_stats.o: $(src)/swp_surfvol_stats.cpp
	$(cpp) $(flag) $(src)/swp_surfvol_stats.cpp

swp_pripart_stats.o: $(src)/swp_pripart_stats.cpp
	$(cpp) $(flag) $(src)/swp_pripart_stats.cpp

# Particle imaging.

swp_imgnode.o: $(src)/swp_imgnode.cpp
	$(cpp) $(flag) $(src)/swp_imgnode.cpp

swp_particle_image.o: $(src)/swp_particle_image.cpp
	$(cpp) $(flag) $(src)/swp_particle_image.cpp

swp_ensemble_imgnode.o: $(src)/swp_ensemble_imgnode.cpp
	$(cpp) $(flag) $(src)/swp_ensemble_imgnode.cpp

swp_ensemble_image.o: $(src)/swp_ensemble_image.cpp
	$(cpp) $(flag) $(src)/swp_ensemble_image.cpp

	
# Mechanism.

swp_particle_model.o: $(src)/swp_particle_model.cpp
	$(cpp) $(flag) $(src)/swp_particle_model.cpp

swp_mech_parser.o: $(src)/swp_mech_parser.cpp  $(inc)/swp_component.h
	$(cpp) $(flag) $(src)/swp_mech_parser.cpp

swp_mechanism.o: $(src)/swp_mechanism.cpp  $(inc)/swp_model_factory.h \
                  $(inc)/swp_model_factory.h  $(inc)/swp_model_factory.h
	$(cpp) $(flag) $(src)/swp_mechanism.cpp

# Active site models.

swp_abf_model.o: $(src)/swp_abf_model.cpp  $(inc)/swp_component.h
	$(cpp) $(flag) $(src)/swp_abf_model.cpp

swp_actsites_model.o: $(src)/swp_actsites_model.cpp
	$(cpp) $(flag) $(src)/swp_actsites_model.cpp

# Aggregation models.

swp_primary.o: $(src)/swp_primary.cpp
	$(cpp) $(flag) $(src)/swp_primary.cpp

swp_surfvol_primary.o: $(src)/swp_surfvol_primary.cpp
	$(cpp) $(flag) $(src)/swp_surfvol_primary.cpp

swp_pripart_primary.o: $(src)/swp_pripart_primary.cpp
	$(cpp) $(flag) $(src)/swp_pripart_primary.cpp

swp_aggmodel_cache.o: $(src)/swp_aggmodel_cache.cpp
	$(cpp) $(flag) $(src)/swp_aggmodel_cache.cpp

swp_surfvol_cache.o: $(src)/swp_surfvol_cache.cpp
	$(cpp) $(flag) $(src)/swp_surfvol_cache.cpp

swp_pripart_cache.o: $(src)/swp_pripart_cache.cpp
	$(cpp) $(flag) $(src)/swp_pripart_cache.cpp

swp_PAH_cache.o: $(src)/swp_PAH_cache.cpp
	$(cpp) $(flag) $(src)/swp_PAH_cache.cpp

swp_PAH_primary.o: $(src)/swp_PAH_primary.cpp
	$(cpp) $(flag) $(src)/swp_PAH_primary.cpp 

swp_PAH_stats.o: $(src)/swp_PAH_stats.cpp
	$(cpp) $(flag) $(src)/swp_PAH_stats.cpp

swp_molecule_evolution.o: $(src)/swp_molecule_evolution.cpp $(inc)/swp_molecule_evolution.h
	$(cpp) $(flag) $(src)/swp_molecule_evolution.cpp

# Sintering model.

swp_sintering_model.o: $(src)/swp_sintering_model.cpp
	$(cpp) $(flag) $(src)/swp_sintering_model.cpp

# Particle sub-models.

swp_submodel.o: $(src)/swp_submodel.cpp
	$(cpp) $(flag) $(src)/swp_submodel.cpp

swp_submodel_cache.o: $(src)/swp_submodel_cache.cpp
	$(cpp) $(flag) $(src)/swp_submodel_cache.cpp

swp_model_factory.o: $(src)/swp_model_factory.cpp
	$(cpp) $(flag) $(src)/swp_model_factory.cpp

swp_arssc_model.o: $(src)/swp_arssc_model.cpp
	$(cpp) $(flag) $(src)/swp_arssc_model.cpp

swp_arssc_cache.o: $(src)/swp_arssc_cache.cpp
	$(cpp) $(flag) $(src)/swp_arssc_cache.cpp

# Processes.

swp_actsites_reaction.o: $(src)/swp_actsites_reaction.cpp
	$(cpp) $(flag) $(src)/swp_actsites_reaction.cpp

swp_coagulation.o: $(src)/swp_coagulation.cpp $(inc)/swp_coagulation.h $(inc)/swp_mechanism.h
	$(cpp) $(flag) $(src)/swp_coagulation.cpp

swp_transcoag.o: $(src)/swp_transcoag.cpp $(inc)/swp_coagulation.h $(inc)/swp_transcoag.h \
                 $(inc)/swp_cell.h $(inc)/swp_mechanism.h $(inc)/swp_params.h
	$(cpp) $(flag) $(src)/swp_transcoag.cpp

swp_addcoag.o: $(src)/swp_addcoag.cpp $(inc)/swp_coagulation.h $(inc)/swp_addcoag.h \
               $(inc)/swp_cell.h $(inc)/swp_mechanism.h $(inc)/swp_params.h
	$(cpp) $(flag) $(src)/swp_addcoag.cpp

swp_condensation.o: $(src)/swp_condensation.cpp
	$(cpp) $(flag) $(src)/swp_condensation.cpp

swp_inception.o: $(src)/swp_inception.cpp
	$(cpp) $(flag) $(src)/swp_inception.cpp

swp_dimer_inception.o: $(src)/swp_dimer_inception.cpp
	$(cpp) $(flag) $(src)/swp_dimer_inception.cpp

swp_pah_inception.o: $(src)/swp_pah_inception.cpp
	$(cpp) $(flag) $(src)/swp_pah_inception.cpp

swp_particle_process.o: $(src)/swp_particle_process.cpp
	$(cpp) $(flag) $(src)/swp_particle_process.cpp

swp_process.o: $(src)/swp_process.cpp
	$(cpp) $(flag) $(src)/swp_process.cpp

swp_process_factory.o: $(src)/swp_process_factory.cpp $(inc)/swp_process_factory.h \
                       $(inc)/swp_transcoag.h $(inc)/swp_addcoag.h
	$(cpp) $(flag) $(src)/swp_process_factory.cpp

swp_surface_reaction.o: $(src)/swp_surface_reaction.cpp
	$(cpp) $(flag) $(src)/swp_surface_reaction.cpp

swp_arssc_inception.o: $(src)/swp_arssc_inception.cpp
	$(cpp) $(flag) $(src)/swp_arssc_inception.cpp

swp_arssc_reaction.o: $(src)/swp_arssc_reaction.cpp
	$(cpp) $(flag) $(src)/swp_arssc_reaction.cpp

swp_arssc_condensation.o: $(src)/swp_arssc_condensation.cpp
	$(cpp) $(flag) $(src)/swp_arssc_condensation.cpp

swp_arssc_process.o: $(src)/swp_arssc_process.cpp
	$(cpp) $(flag) $(src)/swp_arssc_process.cpp

swp_birth_process.o: $(src)/swp_birth_process.cpp
	$(cpp) $(flag) $(src)/swp_birth_process.cpp

swp_death_process.o: $(src)/swp_death_process.cpp
	$(cpp) $(flag) $(src)/swp_death_process.cpp

swp_transport_process.o: $(src)/swp_transport_process.cpp $(inc)/swp_transport_process.h
	$(cpp) $(flag) $(src)/swp_transport_process.cpp

swp_diffusion_process.o: $(src)/swp_diffusion_process.cpp $(inc)/swp_diffusion_process.h $(inc)/swp_transport_process.h
	$(cpp) $(flag) $(src)/swp_diffusion_process.cpp

swp_advection_process.o: $(src)/swp_advection_process.cpp $(inc)/swp_advection_process.h $(inc)/swp_transport_process.h
	$(cpp) $(flag) $(src)/swp_advection_process.cpp

# Solver.

swp_solver.o: $(src)/swp_solver.cpp
	$(cpp) $(flag) $(src)/swp_solver.cpp

# Random number generation

mt19937.o: $(src)/mt19937.cpp $(inc)/mt19937.h
	$(cpp) $(flag) $(src)/mt19937.cpp
	

# TARGET clean: leave just source code.

clean:
	rm -rf $(path)/*.o \#* *~
