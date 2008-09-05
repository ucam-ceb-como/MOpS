# Compiler options.
export debug=0
export profile=0

# C++ compiler.
cpp=g++
ifeq ($(debug),0)
  ifeq ($(profile),0)
    flag=-c -O2 -Iinclude -I../comostrings/include -I../camxml/include \
       -I../sprog/include -Wall
  else
    flag=-c -O2 -Iinclude -I../comostrings/include -I../camxml/include \
       -I../sprog/include -Wall -pg
  endif
else
  flag=-c -O0 -Iinclude -I../comostrings/include -I../camxml/include \
       -I../sprog/include -Wall -g -pg
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
  sprog=../sprog/lib/sprog.a
else
  comostrings=../comostrings/lib/comostrings_d.a
  camxml=../camxml/lib/camxml_d.a
  sprog=../sprog/lib/sprog_d.a
endif

# Target objects.

obj=swp_cell.o swp_component.o swp_ensemble.o swp_ensemble_stats.o \
    swp_mech_parser.o swp_mechanism.o swp_subparticle.o swp_particle.o \
    swp_particle_cache.o swp_particle_stats.o swp_solver.o swp_tracker.o \
    swp_treenode.o swp_abf_model.o swp_actsites_model.o swp_aggmodel_cache.o \
    swp_submodel.o swp_submodel_cache.o swp_model_factory.o \
    swp_surfvol_cache.o swp_surfvol_primary.o swp_primary.o \
    swp_pripart_cache.o swp_pripart_primary.o swp_actsites_reaction.o \
    swp_coagulation.o swp_condensation.o swp_inception.o \
    swp_particle_process.o swp_process.o swp_process_factory.o \
    swp_surface_reaction.o swp_particle_model.o swp_imgnode.o \
    swp_particle_image.o swp_surfvol_stats.o swp_pripart_stats.o \
    swp_sintering_model.o swp_arssc_model.o swp_arssc_inception.o \
    swp_death_process.o swp_birth_process.o swp_arssc_cache.o \
    swp_arssc_reaction.o swp_arssc_condensation.o swp_arssc_process.o

# TARGET all: Compile to a static library (*.a).

all: $(obj)
	mkdir -p $(lib)
ifeq ($(debug),0)
	ar rc $(lib)/sweep.a $(obj) $(comostrings) $(camxml) $(sprog)
else
	ar rc $(lib)/sweep_d.a $(obj) $(comostrings) $(camxml) $(sprog)
endif

# Source code.

# Particle definition.

swp_component.o: $(src)/swp_component.cpp
	$(cpp) $(flag) $(src)/swp_component.cpp

swp_tracker.o: $(src)/swp_tracker.cpp
	$(cpp) $(flag) $(src)/swp_tracker.cpp

swp_particle_cache.o: $(src)/swp_particle_cache.cpp
	$(cpp) $(flag) $(src)/swp_particle_cache.cpp

swp_subparticle.o: $(src)/swp_subparticle.cpp
	$(cpp) $(flag) $(src)/swp_subparticle.cpp

swp_particle.o: $(src)/swp_particle.cpp
	$(cpp) $(flag) $(src)/swp_particle.cpp

# Particle system & ensemble.

swp_cell.o: $(src)/swp_cell.cpp swp_mechanism.o
	$(cpp) $(flag) $(src)/swp_cell.cpp

swp_treenode.o: $(src)/swp_treenode.cpp
	$(cpp) $(flag) $(src)/swp_treenode.cpp

swp_ensemble.o: $(src)/swp_ensemble.cpp
	$(cpp) $(flag) $(src)/swp_ensemble.cpp

# Statistics generation.

swp_particle_stats.o: $(src)/swp_particle_stats.cpp swp_particle.o
	$(cpp) $(flag) $(src)/swp_particle_stats.cpp
	
swp_ensemble_stats.o: $(src)/swp_ensemble_stats.cpp swp_model_factory.o \
                      swp_mechanism.o
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
	
# Mechanism.

swp_particle_model.o: $(src)/swp_particle_model.cpp
	$(cpp) $(flag) $(src)/swp_particle_model.cpp

swp_mech_parser.o: $(src)/swp_mech_parser.cpp swp_component.o swp_tracker.o \
                   swp_inception.o swp_surface_reaction.o \
                   swp_actsites_reaction.o swp_condensation.o \
                   swp_abf_model.o $(camxml) $(comostrings)
	$(cpp) $(flag) $(src)/swp_mech_parser.cpp

swp_mechanism.o: $(src)/swp_mechanism.cpp swp_model_factory.o \
                 swp_process_factory.o swp_abf_model.o
	$(cpp) $(flag) $(src)/swp_mechanism.cpp

# Active site models.

swp_abf_model.o: $(src)/swp_abf_model.cpp swp_component.o
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

swp_coagulation.o: $(src)/swp_coagulation.cpp
	$(cpp) $(flag) $(src)/swp_coagulation.cpp

swp_condensation.o: $(src)/swp_condensation.cpp
	$(cpp) $(flag) $(src)/swp_condensation.cpp

swp_inception.o: $(src)/swp_inception.cpp
	$(cpp) $(flag) $(src)/swp_inception.cpp

swp_particle_process.o: $(src)/swp_particle_process.cpp
	$(cpp) $(flag) $(src)/swp_particle_process.cpp

swp_process.o: $(src)/swp_process.cpp
	$(cpp) $(flag) $(src)/swp_process.cpp

swp_process_factory.o: $(src)/swp_process_factory.cpp
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

# Solver.

swp_solver.o: $(src)/swp_solver.cpp
	$(cpp) $(flag) $(src)/swp_solver.cpp

	

# TARGET clean: leave just source code.

clean:
	rm -rf $(path)/*.o \#* *~
