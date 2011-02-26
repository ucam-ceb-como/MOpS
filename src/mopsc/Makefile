# newMakefile
# Author: Robert I A Patterson <Robert.Patterson@wias-berlin.de>
# Copyright (C) 2011 Robert I A Patterson.

# Makefile for building the mops library as part of larger population
# balance solver projects 

# Licence:
#    This file is part of "mops".
#
#    mops is free software; you can redistribute it and/or
#    modify it under the terms of the GNU General Public License
#    as published by the Free Software Foundation; either version 2
#    of the License, or (at your option) any later version.
#
#    This program is distributed in the hope that it will be useful,
#    but WITHOUT ANY WARRANTY; without even the implied warranty of
#    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#    GNU General Public License for more details.
#
#    You should have received a copy of the GNU General Public License
#    along with this program; if not, write to the Free Software
#    Foundation, Inc., 59 Temple Place - Suite 330, Boston, MA  02111-1307, USA.
#
#  Contact:
#    Prof Markus Kraft
#    Dept of Chemical Engineering
#    University of Cambridge
#    New Museums Site
#    Pembroke Street
#    Cambridge
#    CB2 3RA
#    UK
#
#    Email:       mk306@cam.ac.uk
#    Website:     http://como.cheng.cam.ac.uk

# Compiler options.
export debug=1
export profile=0

export path=.
src=$(path)/source
inc=$(path)/include
librelease=../../lib/release
libdebug=../../lib/debug
libprofile=../../lib/profile

# C++ compiler.
cpp=g++
# Includes.
includes=-Iinclude -I../io/comostrings/include -I../utils/include -I../io/camxml/include \
         -I../sprogc/include -I../geometry/include -I../utils/include -I../sweepc/include \
         -I../odesolvers/cvodes/include -I../odesolvers/cvodes/source/cvodes

ifeq ($(debug),0)
  ifeq ($(profile),0)
    flag=-c -O2 $(includes) -Wall
    lib=$(librelease)
  else
    flag=-c -O2 $(includes)  -Wall -pg
    lib=$(libprofile)
  endif
else
  flag=-c -O0 $(includes) -Wall -ggdb -pg -D_GLIBCXX_DEBUG
  lib=$(libdebug)
endif


# Target objects.
obj=mops_mechanism.o mops_mixture.o \
    mops_particle_solver.o mops_predcor_solver.o mops_psr.o mops_reactor.o \
    mops_reactor_factory.o  mops_ode_solver.o mops_settings_io.o mops_solver.o \
    mops_src_terms.o mops_strang_solver.o mops_timeinterval.o \
    mops_flow_stream.o mops_simulator.o swp_flamesolver.o swp_gas_profile.o mops_gpc_sensitivity.o \
    mops_flux_postprocessor.o mops_rhs_func.o mops_simplesplit_solver.o \
    cvodes_utils.o loi_reduction.o

# TARGET all: Compile to a static library (*.a).
all: $(obj)
	mkdir -p $(lib)
	ar rc $(lib)/libmops.a $(obj)

# Source code.
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

mops_ode_solver.o: $(src)/mops_ode_solver.cpp
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
	
# sensistivity and flux analysis	

mops_flux_postprocessor.o: $(src)/mops_flux_postprocessor.cpp
	$(cpp) $(flag) $(src)/mops_flux_postprocessor.cpp

mops_gpc_sensitivity.o: $(src)/mops_gpc_sensitivity.cpp
	$(cpp) $(flag) $(src)/mops_gpc_sensitivity.cpp
  

mops_rhs_func.o: $(src)/mops_rhs_func.cpp
	$(cpp) $(flag) $(src)/mops_rhs_func.cpp

mops_simplesplit_solver.o: $(src)/mops_simplesplit_solver.cpp
	$(cpp) $(flag) $(src)/mops_simplesplit_solver.cpp

loi_reduction.o: $(src)/loi_reduction.cpp
	$(cpp) $(flag) $(src)/loi_reduction.cpp

# Auxilliary classes.

mops_src_terms.o: $(src)/mops_src_terms.cpp
	$(cpp) $(flag) $(src)/mops_src_terms.cpp

swp_gas_profile.o: $(src)/swp_gas_profile.cpp
	$(cpp) $(flag) $(src)/swp_gas_profile.cpp

cvodes_utils.o: $(src)/cvodes_utils.cpp
	$(cpp) $(flag) $(src)/cvodes_utils.cpp


# TARGET clean: leave just source code.

clean:
	rm -rf $(path)/*.o \#* *~

cleanlibs:
	rm -rf $(path)/*.o \#* *~ $(librelease)/libmops.a $(libdebug)/libmops.a $(libprofile)/libmops.a
