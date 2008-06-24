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
       -I../sprog/include -Wall -ggdb -pg
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

obj=swp_cell.o swp_component.o swp_ensemble.o swp_ensemblestats.o \
    swp_mech_parser.o swp_mechanism.o swp_particle.o swp_particledata.o \
    swp_particlestats.o swp_solver.o swp_tracker.o swp_treenode.o \
    swp_abfmodel.o swp_activesites_model.o swp_coagmodel.o swp_coagmodeldata.o \
    swp_model.o swp_modeldata.o swp_modelfactory.o swp_pointcontactdata.o \
    swp_pointcontactmodel.o swp_primary.o swp_pripartdata.o \
    swp_pripartmodel.o swp_activesites_reaction.o swp_coagulation.o \
    swp_condensation.o swp_inception.o swp_particleprocess.o swp_process.o \
    swp_processfactory.o swp_surfacereaction.o

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

swp_particledata.o: $(src)/swp_particledata.cpp swp_model.o swp_component.o
	$(cpp) $(flag) $(src)/swp_particledata.cpp

swp_particle.o: $(src)/swp_particle.cpp swp_particledata.o swp_model.o
	$(cpp) $(flag) $(src)/swp_particle.cpp

# Particle system & ensemble.

swp_cell.o: $(src)/swp_cell.cpp swp_mechanism.o
	$(cpp) $(flag) $(src)/swp_cell.cpp

swp_treenode.o: $(src)/swp_treenode.cpp
	$(cpp) $(flag) $(src)/swp_treenode.cpp

swp_ensemble.o: $(src)/swp_ensemble.cpp
	$(cpp) $(flag) $(src)/swp_ensemble.cpp

# Statistics generation.

swp_particlestats.o: $(src)/swp_particlestats.cpp swp_particle.o
	$(cpp) $(flag) $(src)/swp_particlestats.cpp
	
swp_ensemblestats.o: $(src)/swp_ensemblestats.cpp swp_modelfactory.o \
                     swp_mechanism.o
	$(cpp) $(flag) $(src)/swp_ensemblestats.cpp

# Mechanism.

swp_mech_parser.o: $(src)/swp_mech_parser.cpp swp_component.o swp_tracker.o \
                   swp_inception.o swp_surfacereaction.o \
                   swp_activesites_reaction.o swp_condensation.o \
                   swp_abfmodel.o $(camxml) $(comostrings)
	$(cpp) $(flag) $(src)/swp_mech_parser.cpp

swp_mechanism.o: $(src)/swp_mechanism.cpp swp_modelfactory.o \
                 swp_processfactory.o swp_abfmodel.o
	$(cpp) $(flag) $(src)/swp_mechanism.cpp

# Particle models.

swp_abfmodel.o: $(src)/swp_abfmodel.cpp swp_component.o
	$(cpp) $(flag) $(src)/swp_abfmodel.cpp

swp_activesites_model.o: $(src)/swp_activesites_model.cpp
	$(cpp) $(flag) $(src)/swp_activesites_model.cpp

swp_coagmodel.o: $(src)/swp_coagmodel.cpp swp_particle.o
	$(cpp) $(flag) $(src)/swp_coagmodel.cpp

swp_coagmodeldata.o: $(src)/swp_coagmodeldata.cpp swp_coagmodel.o
	$(cpp) $(flag) $(src)/swp_coagmodeldata.cpp

swp_model.o: $(src)/swp_model.cpp
	$(cpp) $(flag) $(src)/swp_model.cpp

swp_modeldata.o: $(src)/swp_modeldata.cpp
	$(cpp) $(flag) $(src)/swp_modeldata.cpp

swp_modelfactory.o: $(src)/swp_modelfactory.cpp swp_coagmodeldata.o \
                    swp_pointcontactdata.o swp_pripartdata.o swp_coagmodel.o \
                    swp_particlestats.o swp_pointcontactmodel.o \
                    swp_pripartmodel.o swp_abfmodel.o
	$(cpp) $(flag) $(src)/swp_modelfactory.cpp

swp_pointcontactdata.o: $(src)/swp_pointcontactdata.cpp
	$(cpp) $(flag) $(src)/swp_pointcontactdata.cpp

swp_pointcontactmodel.o: $(src)/swp_pointcontactmodel.cpp
	$(cpp) $(flag) $(src)/swp_pointcontactmodel.cpp

swp_primary.o: $(src)/swp_primary.cpp
	$(cpp) $(flag) $(src)/swp_primary.cpp

swp_pripartdata.o: $(src)/swp_pripartdata.cpp
	$(cpp) $(flag) $(src)/swp_pripartdata.cpp

swp_pripartmodel.o: $(src)/swp_pripartmodel.cpp
	$(cpp) $(flag) $(src)/swp_pripartmodel.cpp

# Processes.

swp_activesites_reaction.o: $(src)/swp_activesites_reaction.cpp
	$(cpp) $(flag) $(src)/swp_activesites_reaction.cpp

swp_coagulation.o: $(src)/swp_coagulation.cpp
	$(cpp) $(flag) $(src)/swp_coagulation.cpp

swp_condensation.o: $(src)/swp_condensation.cpp
	$(cpp) $(flag) $(src)/swp_condensation.cpp

swp_inception.o: $(src)/swp_inception.cpp
	$(cpp) $(flag) $(src)/swp_inception.cpp

swp_particleprocess.o: $(src)/swp_particleprocess.cpp
	$(cpp) $(flag) $(src)/swp_particleprocess.cpp

swp_process.o: $(src)/swp_process.cpp
	$(cpp) $(flag) $(src)/swp_process.cpp

swp_processfactory.o: $(src)/swp_processfactory.cpp
	$(cpp) $(flag) $(src)/swp_processfactory.cpp

swp_surfacereaction.o: $(src)/swp_surfacereaction.cpp
	$(cpp) $(flag) $(src)/swp_surfacereaction.cpp

# Solver.

swp_solver.o: $(src)/swp_solver.cpp
	$(cpp) $(flag) $(src)/swp_solver.cpp

	

# TARGET clean: leave just source code.

clean:
	rm -rf $(path)/*.o \#* *~
