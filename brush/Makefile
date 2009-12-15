# Compiler options.
export debug=1
export profile=0

# C++ compiler.
cpp=g++
ifeq ($(debug),0)
  ifeq ($(profile),0)
    flag=-c -O2 -Iinclude -I../comostrings/include -I../camxml/include \
       -I../geometry/include -I../utils/include -I../sprogc/include \
       -I../sweepc/include -I../mopsc/include -Wall
    flink=
  else
    flag=-c -O2 -Iinclude -I../comostrings/include -I../camxml/include \
       -I../geometry/include -I../utils/include -I../sprogc/include \
       -I../sweepc/include -I../mopsc/include -Wall -pg
    flink=-pg
  endif
else
  flag=-c -O0 -Iinclude -I../comostrings/include -I../camxml/include \
       -I../geometry/include -I../utils/include -I../sprogc/include \
       -I../sweepc/include -I../mopsc/include -Wall -ggdb
  flink=
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
  sweep=../sweepc/lib/sweep.a
  mops=../mopsc/lib/mops.a
  cvodes=../cvodes/lib/cvodes.a
  geom=../geometry/lib/geometry.a
  brushlib=./lib/brush.a
else
  comostrings=../comostrings/lib/comostrings_d.a
  camxml=../camxml/lib/camxml_d.a
  sprog=../sprogc/lib/sprog_d.a
  sweep=../sweepc/lib/sweep_d.a
  mops=../mopsc/lib/mops_d.a
  cvodes=../cvodes/lib/cvodes_d.a
  geom=../geometry/lib/geometry_d.a
  brushlib=./lib/brush_d.a
endif
linklibs=$(mops) $(sweep) $(geom) $(sprog) $(camxml) $(comostrings) $(cvodes)
# Target objects.

obj=reactor1d.o pred_corr_solver.o reset_chemistry.o simulator.o settings_io.o

# Default target - executable
all: $(brushlib) brush.o
	mkdir -p bin
ifeq ($(debug),0)
	$(cpp) $(flink) brush.o $(brushlib) $(linklibs) -o bin/brush.x
else
	$(cpp) $(flink) brush.o $(brushlib) $(linklibs) -o bin/brush_d.x
endif

# Bring everything except the main function into a library
$(brushlib): $(obj)
	mkdir -p lib
ifeq ($(debug),0)
	ar -rc lib/brush.a $(obj)
else
	ar -rc lib/brush_d.a $(obj)
endif

# Build the objects that contain the code for the real work
brush.o: $(src)/brush.cpp
	$(cpp) $(flag) $(src)/brush.cpp

reactor1d.o: $(src)/reactor1d.cpp $(inc)/reactor1d.h $(inc)/brush_params.h
	$(cpp) $(flag) $(src)/reactor1d.cpp

pred_corr_solver.o: $(src)/pred_corr_solver.cpp $(inc)/pred_corr_solver.h $(inc)/brush_params.h
	$(cpp) $(flag) $(src)/pred_corr_solver.cpp

reset_chemistry.o: $(src)/reset_chemistry.cpp $(inc)/reset_chemistry.h $(inc)/reactor1d.h $(inc)/brush_params.h
	$(cpp) $(flag) $(src)/reset_chemistry.cpp

simulator.o: $(src)/simulator.cpp $(inc)/simulator.h $(inc)/brush_params.h
	$(cpp) $(flag) $(src)/simulator.cpp

settings_io.o: $(src)/settings_io.cpp $(inc)/settings_io.h $(inc)/brush_params.h
	$(cpp) $(flag) $(src)/settings_io.cpp

# leave just source code.
clean:
	rm -rf $(path)/*.o \#* *~ $(lib)/*.a $(path)/bin/*.x

sweepclean:
	$(MAKE) -f ../sweepc/Makefile -C ../sweepc clean
	$(MAKE) -f ../sweepc/Makefile -C ../sweepc debug=$(debug) profile=$(profile)

sweep:
	$(MAKE) -f ../sweepc/Makefile -C ../sweepc debug=$(debug) profile=$(profile)

mopsclean:
	$(MAKE) -f ../mopsc/Makefile -C ../mopsc clean
	$(MAKE) -f ../mopsc/Makefile -C ../mopsc debug=$(debug) profile=$(profile) mopslib

mopslib:
	$(MAKE) -f ../mopsc/Makefile -C ../mopsc debug=$(debug) profile=$(profile) mopslib

geom:
	$(MAKE) -f ../geometry/Makefile -C ../geometry debug=$(debug) profile=$(profile)

geomclean:
	$(MAKE) -f ../geometry/Makefile -C ../geometry clean
	$(MAKE) -f ../geometry/Makefile -C ../geometry debug=$(debug) profile=$(profile)

