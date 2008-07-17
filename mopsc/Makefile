# Compiler options.
export debug=0
export profile=0

# C++ compiler.
cpp=g++
ifeq ($(debug),0)
  ifeq ($(profile),0)
    flag=-c -O2 -Iinclude -I../comostrings/include -I../camxml/include \
         -I../sprog/include -I../sweep/include -I../cvode/include \
         -I../cvode/source -Wall
    flink= 
  else
    flag=-c -O2 -Iinclude -I../comostrings/include -I../camxml/include \
         -I../sprog/include -I../sweep/include -I../cvode/include \
         -I../cvode/source -Wall -pg
    flink=-pg
  endif
else
  flag=-c -O0 -Iinclude -I../comostrings/include -I../camxml/include \
       -I../sprog/include -I../sweep/include -I../cvode/include \
       -I../cvode/source -Wall -ggdb -pg
  flink=-pg 
endif

# Source directories.

src=source
inc=include

# Additional libraries.

ifeq ($(debug),0)
  comostrings=../comostrings/lib/comostrings.a
  camxml=../camxml/lib/camxml.a
  sprog=../sprog/lib/sprog.a
  sweep=../sweep/lib/sweep.a
  cvode=../cvode/lib/cvode.a
else
  comostrings=../comostrings/lib/comostrings_d.a
  camxml=../camxml/lib/camxml_d.a
  sprog=../sprog/lib/sprog_d.a
  sweep=../sweep/lib/sweep_d.a
  cvode=../cvode/lib/cvode_d.a
endif
linklibs=$(comostrings) $(camxml) $(sprog) $(sweep) $(cvode)

# Target objects.

obj=mops.o mops_mechanism.o mops_mixture.o \
    mops_particle_solver.o mops_predcor_solver.o mops_psr.o mops_reactor.o \
    mops_reactor_factory.o  mops_ode_solver.o mops_settings_io.o mops_solver.o \
    mops_src_terms.o mops_strang_solver.o mops_timeinterval.o \
    mops_flow_stream.o mops_simulator.o swp_flamesolver.o swp_gas_profile.o

# TARGET all: Compile to a static library (*.a).

all: comostringslib camxmllib sproglib sweeplib cvodelib $(obj)
	mkdir -p bin
ifeq ($(debug),0)
  ifeq ($(profile),0)
		$(cpp) $(flink) $(obj) $(linklibs) -o bin/mops.x
  else
		$(cpp) $(flink) $(obj) $(linklibs) -o bin/mops_p.x
  endif 
else
	$(cpp) $(flink) $(obj) $(linklibs) -o bin/mops_d.x
endif

# Required libraries.
comostringslib:
	$(MAKE) -f ../comostrings/Makefile -C ../comostrings debug=$(debug) \
  profile=$(profile)
camxmllib:
	$(MAKE) -f ../camxml/Makefile -C ../camxml debug=$(debug) profile=$(profile)
sproglib:
	$(MAKE) -f ../sprog/Makefile -C ../sprog debug=$(debug) profile=$(profile)
sweeplib:
	$(MAKE) -f ../sweep/Makefile -C ../sweep debug=$(debug) profile=$(profile)
cvodelib:
	$(MAKE) -f ../cvode/Makefile -C ../cvode debug=$(debug) profile=$(profile)

# Source code.

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

mops_ode_solver.o: $(src)/mops_ode_solver.cpp ../cvode/nvector_serial.o
	$(cpp) $(flag) $(src)/mops_ode_solver.cpp

mops_solver.o: $(src)/mops_solver.cpp
	$(cpp) $(flag) $(src)/mops_solver.cpp
	
mops_particle_solver.o: $(src)/mops_particle_solver.cpp
	$(cpp) $(flag) $(src)/mops_particle_solver.cpp

mops_strang_solver.o: $(src)/mops_strang_solver.cpp
	$(cpp) $(flag) $(src)/mops_strang_solver.cpp

swp_flamesolver.o: $(src)/swp_flamesolver.cpp
	$(cpp) $(flag) $(src)/swp_flamesolver.cpp

mops_predcor_solver.o: $(src)/mops_predcor_solver.cpp
	$(cpp) $(flag) $(src)/mops_predcor_solver.cpp

mops_simulator.o: $(src)/mops_simulator.cpp
	$(cpp) $(flag) $(src)/mops_simulator.cpp

# Auxilliary classes.

mops_src_terms.o: $(src)/mops_src_terms.cpp
	$(cpp) $(flag) $(src)/mops_src_terms.cpp

swp_gas_profile.o: $(src)/swp_gas_profile.cpp
	$(cpp) $(flag) $(src)/swp_gas_profile.cpp


# TARGET clean: leave just source code.

clean:
	$(MAKE) clean -f ../comostrings/Makefile -C ../comostrings
	$(MAKE) clean -f ../camxml/Makefile -C ../camxml
	$(MAKE) clean -f ../sprog/Makefile -C ../sprog
	$(MAKE) clean -f ../sweep/Makefile -C ../sweep
	$(MAKE) clean -f ../cvode/Makefile -C ../cvode
	rm -rf *.o \#* *~
