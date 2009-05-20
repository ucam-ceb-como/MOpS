# Compiler options.
export debug=0
export profile=0

# C++ compiler.
cpp=g++
ifeq ($(debug),0)
  ifeq ($(profile),0)
    flag=-c -O2 -Iinclude -I../comostrings/include -I../camxml/include \
         -I../sprog/include -I../cvode/include \
			-I ../ida/include -I ../kinsol/include -I../radau/include \
         -I../cvode/source -Wall
    flink= 
  else
    flag=-c -O2 -Iinclude -I../comostrings/include -I../camxml/include \
         -I../sprog/include -I../cvode/include \
			-I ../ida/include -I ../kinsol/include -I../radau/include \
         -I../cvode/source -Wall -pg
    flink=-pg
  endif
else
  flag=-c -O2  -Iinclude -I../comostrings/include -I../camxml/include \
       -I../sprog/include -I../cvode/include \
		 -I ../ida/include -I ../kinsol/include -I../radau/include \
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
  ida=../ida/lib/ida.a
  kinsol=../kinsol/lib/kinsol.a
  radau=../radau/lib/radau.a
else
  comostrings=../comostrings/lib/comostrings_d.a
  camxml=../camxml/lib/camxml_d.a
  sprog=../sprog/lib/sprog_d.a
  cvode=../cvode/lib/cvode_d.a
  ida=../ida/lib/ida_d.a
  kinsol=../kinsol/lib/kinsol_d.a
  radau=../radau/lib/radau_d.a
endif
linklibs=$(comostrings) $(camxml) $(sprog) $(cvode) $(ida) $(kinsol) $(radau)

# Target objects.

obj=   cam_kernel.o cam_read.o cam_setup.o cam_geometry.o \
		 cam_control.o cam_conveter.o cam_math.o cam_conc.o \
		 cam_admin.o cam_boundary.o cam_configuration.o cam_profile.o\
		 cam_reporter.o cam_models.o cam_plug.o cam_premix.o cam_residual.o \
		 stagflow.o \
		 cvode_wrapper.o ida_wrapper.o kinsol_wrapper.o radau_wrapper.o


# TARGET all: Compile to a static library (*.a).

all: comostringslib camxmllib sproglib cvodelib idalib kinlib radaulib $(obj)
	mkdir -p bin
ifeq ($(debug),0)
  ifeq ($(profile),0)
		$(cpp) $(flink) $(obj) $(linklibs) -o bin/camflow
  else
		$(cpp) $(flink) $(obj) $(linklibs) -o bin/camflow_p
  endif 
else
	$(cpp) $(flink) $(obj) $(linklibs) -o bin/camflow_d
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
idalib:
	$(MAKE) -f ../ida/Makefile -C ../ida debug=$(debug) profile=$(profile)
kinlib:
	$(MAKE) -f ../kinsol/Makefile -C ../kinsol debug=$(debug) profile=$(profile)
radaulib:
	$(MAKE) -f ../radau/Makefile -C ../radau debug=$(debug) profile=$(profile)

# Source code.

# Driver program.

cam_kernel.o: $(src)/cam_kernel.cpp
	$(cpp) $(flag) $(src)/cam_kernel.cpp
cam_read.o: $(src)/cam_read.cpp
	$(cpp) $(flag) $(src)/cam_read.cpp
cam_setup.o: $(src)/cam_setup.cpp
	$(cpp) $(flag) $(src)/cam_setup.cpp
cam_geometry.o: $(src)/cam_geometry.cpp
	$(cpp) $(flag) $(src)/cam_geometry.cpp
cam_control.o: $(src)/cam_control.cpp
	$(cpp) $(flag) $(src)/cam_control.cpp
cam_conveter.o: $(src)/cam_conveter.cpp
	$(cpp) $(flag) $(src)/cam_conveter.cpp
cam_math.o: $(src)/cam_math.cpp
	$(cpp) $(flag) $(src)/cam_math.cpp
cam_residual.o: $(src)/cam_residual.cpp
	$(cpp) $(flag) $(src)/cam_residual.cpp
cam_conc.o: $(src)/cam_conc.cpp
	$(cpp) $(flag) $(src)/cam_conc.cpp
cam_admin.o: $(src)/cam_admin.cpp
	$(cpp) $(flag) $(src)/cam_admin.cpp
cam_boundary.o: $(src)/cam_boundary.cpp
	$(cpp) $(flag) $(src)/cam_boundary.cpp
cam_configuration.o: $(src)/cam_configuration.cpp
	$(cpp) $(flag) $(src)/cam_configuration.cpp
cam_profile.o: $(src)/cam_profile.cpp
	$(cpp) $(flag) $(src)/cam_profile.cpp
cam_reporter.o: $(src)/cam_reporter.cpp
	$(cpp) $(flag) $(src)/cam_reporter.cpp
cam_models.o: $(src)/cam_models.cpp
	$(cpp) $(flag) $(src)/cam_models.cpp
cam_plug.o: $(src)/cam_plug.cpp
	$(cpp) $(flag) $(src)/cam_plug.cpp
cam_premix.o: $(src)/cam_premix.cpp
	$(cpp) $(flag) $(src)/cam_premix.cpp
stagflow.o: $(src)/stagflow.cpp
	$(cpp) $(flag) $(src)/stagflow.cpp
cvode_wrapper.o: $(src)/cvode_wrapper.cpp
	$(cpp) $(flag) $(src)/cvode_wrapper.cpp
ida_wrapper.o: $(src)/ida_wrapper.cpp
	$(cpp) $(flag) $(src)/ida_wrapper.cpp
kinsol_wrapper.o: $(src)/kinsol_wrapper.cpp
	$(cpp) $(flag) $(src)/kinsol_wrapper.cpp
radau_wrapper.o: $(src)/radau_wrapper.cpp
	$(cpp) $(flag) $(src)/radau_wrapper.cpp


# Combined mechanisms.


# TARGET clean: leave just source code.

clean:
	$(MAKE) clean -f ../comostrings/Makefile -C ../comostrings
	$(MAKE) clean -f ../camxml/Makefile -C ../camxml
	$(MAKE) clean -f ../sprog/Makefile -C ../sprog
	$(MAKE) clean -f ../cvode/Makefile -C ../cvode
	$(MAKE) clean -f ../ida/Makefile -C ../ida
	$(MAKE) clean -f ../kinsol/Makefile -C ../kinsol
	$(MAKE) clean -f ../radau/Makefile -C ../radau
	rm -rf *.o \#* *~
