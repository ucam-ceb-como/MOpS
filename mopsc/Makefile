# Compiler options.
export debug=0
export profile=0

# C++ compiler.
cpp=g++
#cpp=mpicxx -DUSE_MPI
ifeq ($(debug),0)
  ifeq ($(profile),0)
    flag=-c -O2 -Iinclude -I../comostrings/include -I../camxml/include \
         -I../sprogc/include -I../sweepc/include -I../cvodes/include \
         -I../cvodes/source/cvodes -I../geometry/include -Wall 
    flink= 
  else
    flag=-c -O2 -Iinclude -I../comostrings/include -I../camxml/include \
         -I../sprogc/include -I../sweepc/include -I../cvodes/include \
         -I../cvodes/source/cvodes -I../geometry/include -Wall -pg
    flink=-pg
  endif
else
  flag=-c -O0 -Iinclude -I../comostrings/include -I../camxml/include \
       -I../sprogc/include -I../sweepc/include -I../cvodes/include \
       -I../cvodes/source/cvodes -I../geometry/include -Wall -ggdb -pg
  flink=-pg 
endif

# Source directories.

src=source
inc=include
lib=lib

# Additional libraries.

ifeq ($(debug),0)
  comostrings=../comostrings/lib/comostrings.a
  camxml=../camxml/lib/camxml.a
  sprog=../sprogc/lib/sprog.a
  sweep=../sweepc/lib/sweep.a
  cvodes=../cvodes/lib/cvodes.a
  geom=../geometry/lib/geometry.a
else
  comostrings=../comostrings/lib/comostrings_d.a
  camxml=../camxml/lib/camxml_d.a
  sprog=../sprogc/lib/sprog_d.a
  sweep=../sweepc/lib/sweep_d.a
  cvodes=../cvodes/lib/cvodes_d.a
  geom=../geometry/lib/geometry_d.a
endif
linklibs=$(comostrings) $(camxml) $(sprog) $(sweep) $(cvodes) $(geom)

# Target objects.

obj=mops_mechanism.o mops_mixture.o \
    mops_particle_solver.o mops_predcor_solver.o mops_psr.o mops_reactor.o \
    mops_reactor_factory.o  mops_ode_solver.o mops_settings_io.o mops_solver.o \
    mops_src_terms.o mops_strang_solver.o mops_timeinterval.o \
    mops_flow_stream.o mops_simulator.o swp_flamesolver.o swp_gas_profile.o mops_gpc_sensitivity.o \
    mops_flux_postprocessor.o mops_rhs_func.o mops_simplesplit_solver.o  swp_PAHsolver.o

# TARGET all: Compile to a static library (*.a).

all: comostringslib camxmllib sproglib sweeplib cvodeslib geomlib $(obj) mops.o
	mkdir -p bin
ifeq ($(debug),0)
  ifeq ($(profile),0)
		$(cpp) $(flink) $(obj) mops.o $(linklibs) -o bin/mops.x -static
  else
		$(cpp) $(flink) $(obj) mops.o $(linklibs) -o bin/mops_p.x -static
  endif 
else
	$(cpp) $(flink) $(obj) mops.o $(linklibs) -o bin/mops_d.x -static
endif

32bit: comostringslib camxmllib sproglib sweeplib cvodeslib $(obj)
	mkdir -p bin
ifeq ($(debug),0)
  ifeq ($(profile),0)
		$(cpp) $(flink) $(obj) mops.o $(linklibs) -o bin/mops32.x -static
  else
		$(cpp) $(flink) $(obj) mops.o $(linklibs) -o bin/mops32_p.x -static
  endif
else
		$(cpp) $(flink) $(obj) mops.o $(linklibs) -o bin/mops32_d.x -static
endif


# Required libraries.
comostringslib:
	$(MAKE) -f ../comostrings/Makefile -C ../comostrings debug=$(debug) \
  profile=$(profile)
camxmllib:
	$(MAKE) -f ../camxml/Makefile -C ../camxml debug=$(debug) profile=$(profile)
sproglib:
	$(MAKE) -f ../sprogc/Makefile -C ../sprogc debug=$(debug) profile=$(profile)
sweeplib:
	$(MAKE) -f ../sweepc/Makefile -C ../sweepc debug=$(debug) profile=$(profile)
cvodeslib:
	$(MAKE) -f ../cvodes/Makefile -C ../cvodes debug=$(debug) profile=$(profile)
geomlib:
	$(MAKE) -f ../geometry/Makefile -C ../geometry debug=$(debug) profile=$(profile)

# Mops library
mopslib: $(obj) comostringslib camxmllib sproglib sweeplib cvodeslib
	mkdir -p lib
ifeq ($(debug),0)
	ar -rcs $(lib)/mops.a $(obj)
else
	ar -rcs $(lib)/mops_d.a $(obj)
endif

# Driver program.

mops.o: $(src)/mops.cpp
	$(cpp) $(flag) $(src)/mops.cpp

mops_settings_io.o: $(src)/mops_settings_io.cpp
	$(cpp) $(flag) $(src)/mops_settings_io.cpp

# Combined mechanisms.

mops_mechanism.o: $(src)/mops_mechanism.cpp
	$(cpp) $(flag) $(src)/mops_mechanism.cpp

# Systems (mixtures and reactors).

mops_mixture.o: $(src)/mops_mixture.cpp
	$(cpp) $(flag) $(src)/mops_mixture.cpp

mops_flow_stream.o: $(src)/mops_flow_stream.cpp
	$(cpp) $(flag) $(src)/mops_flow_stream.cpp

mops_reactor.o: $(src)/mops_reactor.cpp
	$(cpp) $(flag) $(src)/mops_reactor.cpp

mops_psr.o: $(src)/mops_psr.cpp
	$(cpp) $(flag) $(src)/mops_psr.cpp

mops_reactor_factory.o: $(src)/mops_reactor_factory.cpp
	$(cpp) $(flag) $(src)/mops_reactor_factory.cpp

# Solvers.

mops_timeinterval.o: $(src)/mops_timeinterval.cpp
	$(cpp) $(flag) $(src)/mops_timeinterval.cpp

mops_ode_solver.o: $(src)/mops_ode_solver.cpp ../cvodes/nvector_serial.o
	$(cpp) $(flag) $(src)/mops_ode_solver.cpp

mops_solver.o: $(src)/mops_solver.cpp
	$(cpp) $(flag) $(src)/mops_solver.cpp
	
mops_particle_solver.o: $(src)/mops_particle_solver.cpp
	$(cpp) $(flag) $(src)/mops_particle_solver.cpp

mops_strang_solver.o: $(src)/mops_strang_solver.cpp
	$(cpp) $(flag) $(src)/mops_strang_solver.cpp

swp_flamesolver.o: $(src)/swp_flamesolver.cpp
	$(cpp) $(flag) $(src)/swp_flamesolver.cpp

swp_PAHsolver.o: $(src)/swp_PAHsolver.cpp
	$(cpp) $(flag) $(src)/swp_PAHsolver.cpp

mops_predcor_solver.o: $(src)/mops_predcor_solver.cpp
	$(cpp) $(flag) $(src)/mops_predcor_solver.cpp

mops_simulator.o: $(src)/mops_simulator.cpp
	$(cpp) $(flag) $(src)/mops_simulator.cpp
	
# sensistivity and flux analysis	

mops_flux_postprocessor.o: $(src)/mops_flux_postprocessor.cpp
	$(cpp) $(flag) $(src)/mops_flux_postprocessor.cpp

mops_gpc_sensitivity.o: $(src)/mops_gpc_sensitivity.cpp
	$(cpp) $(flag) $(src)/mops_gpc_sensitivity.cpp
  

mops_rhs_func.o: $(src)/mops_rhs_func.cpp
	$(cpp) $(flag) $(src)/mops_rhs_func.cpp

mops_simplesplit_solver.o: $(src)/mops_simplesplit_solver.cpp
	$(cpp) $(flag) $(src)/mops_simplesplit_solver.cpp



# Auxilliary classes.

mops_src_terms.o: $(src)/mops_src_terms.cpp
	$(cpp) $(flag) $(src)/mops_src_terms.cpp

swp_gas_profile.o: $(src)/swp_gas_profile.cpp
	$(cpp) $(flag) $(src)/swp_gas_profile.cpp


# TARGET clean: leave just source code.

cleanlibs:
	$(MAKE) clean -f ../comostrings/Makefile -C ../comostrings
	$(MAKE) clean -f ../camxml/Makefile -C ../camxml
	$(MAKE) clean -f ../sprogc/Makefile -C ../sprogc
	$(MAKE) clean -f ../sweepc/Makefile -C ../sweepc
	$(MAKE) clean -f ../cvodes/Makefile -C ../cvodes
	$(MAKE) clean -f ../geometry/Makefile -C ../geometry
	rm -rf *.o \#* *~

cleanparticles:
	$(MAKE) clean -f ../sweepc/Makefile -C ../sweepc
	rm -rf *.o \#* *~

clean:
	rm -rf *.o \#* *~ $(lib)/*.a

