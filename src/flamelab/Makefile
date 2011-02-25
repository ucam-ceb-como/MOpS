# Compiler options.
export debug=0
export profile=0

# C++ compiler.
cpp=g++
ifeq ($(debug),0)
  ifeq ($(profile),0)
    flag=-c -O2 -Iinclude -I../comostrings/include -I../camxml/include \
         -I../sprog/include -I../cvode/include \
         -I../cvode/source -Wall
    flink= 
  else
    flag=-c -O2 -Iinclude -I../comostrings/include -I../camxml/include \
         -I../sprog/include -I../cvode/include \
         -I../cvode/source -Wall -pg
    flink=-pg
  endif
else
  flag=-c -O2  -Iinclude -I../comostrings/include -I../camxml/include \
       -I../sprog/include -I../cvode/include \
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
  cvode=../cvode/lib/cvode.a
else
  comostrings=../comostrings/lib/comostrings_d.a
  camxml=../camxml/lib/camxml_d.a
  sprog=../sprog/lib/sprog_d.a
  cvode=../cvode/lib/cvode_d.a
endif
linklibs=$(comostrings) $(camxml) $(sprog) $(cvode)

# Target objects.

obj= flame.o  fl_solver_control.o \
	  fl_single_cell.o fl_reactor_controls.o fl_reactor.o \
	  fl_process_conditions.o fl_premix.o fl_io.o\
	  fl_initial.o fl_geometry.o fl_counter_diffusion.o\
	  fl_cell_interface.o 

# TARGET all: Compile to a static library (*.a).

all: comostringslib camxmllib sproglib cvodelib $(obj)
	mkdir -p bin
ifeq ($(debug),0)
  ifeq ($(profile),0)
		$(cpp) $(flink) $(obj) $(linklibs) -o bin/flame.x
  else
		$(cpp) $(flink) $(obj) $(linklibs) -o bin/flame_p.x
  endif 
else
	$(cpp) $(flink) $(obj) $(linklibs) -o bin/flame_d.x
endif

# Required libraries.
comostringslib:
	$(MAKE) -f ../comostrings/Makefile -C ../comostrings debug=$(debug) \
  profile=$(profile)
camxmllib:
	$(MAKE) -f ../camxml/Makefile -C ../camxml debug=$(debug) profile=$(profile)
sproglib:
	$(MAKE) -f ../sprog/Makefile -C ../sprog debug=$(debug) profile=$(profile)
cvodelib:
	$(MAKE) -f ../cvode/Makefile -C ../cvode debug=$(debug) profile=$(profile)

# Source code.

# Driver program.

flame.o: $(src)/flame.cpp
	$(cpp) $(flag) $(src)/flame.cpp

# Combined mechanisms.

fl_solver_control.o: $(src)/fl_solver_control.cpp
	$(cpp) $(flag) $(src)/fl_solver_control.cpp

# Systems (mixtures and reactors).

fl_single_cell.o: $(src)/fl_single_cell.cpp
	$(cpp) $(flag) $(src)/fl_single_cell.cpp

fl_reactor_controls.o: $(src)/fl_reactor_controls.cpp
	$(cpp) $(flag) $(src)/fl_reactor_controls.cpp

fl_reactor.o: $(src)/fl_reactor.cpp
	$(cpp) $(flag) $(src)/fl_reactor.cpp

fl_process_conditions.o: $(src)/fl_process_conditions.cpp
	$(cpp) $(flag) $(src)/fl_process_conditions.cpp

fl_premix.o: $(src)/fl_premix.cpp
	$(cpp) $(flag) $(src)/fl_premix.cpp

# Solvers.

fl_io.o: $(src)/fl_io.cpp
	$(cpp) $(flag) $(src)/fl_io.cpp

fl_initial.o: $(src)/fl_initial.cpp
	$(cpp) $(flag) $(src)/fl_initial.cpp

fl_geometry.o: $(src)/fl_geometry.cpp
	$(cpp) $(flag) $(src)/fl_geometry.cpp
	
fl_counter_diffusion.o: $(src)/fl_counter_diffusion.cpp
	$(cpp) $(flag) $(src)/fl_counter_diffusion.cpp

fl_cell_interface.o: $(src)/fl_cell_interface.cpp
	$(cpp) $(flag) $(src)/fl_cell_interface.cpp


# TARGET clean: leave just source code.

clean:
	$(MAKE) clean -f ../comostrings/Makefile -C ../comostrings
	$(MAKE) clean -f ../camxml/Makefile -C ../camxml
	$(MAKE) clean -f ../sprog/Makefile -C ../sprog
	$(MAKE) clean -f ../cvode/Makefile -C ../cvode
	rm -rf *.o \#* *~
