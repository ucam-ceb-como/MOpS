# Compiler options.
export debug=1
export profile=0

# C++ compiler.
cpp=g++
ifeq ($(debug),0)
  ifeq ($(profile),0)
    flag=-c -O2 -Iinclude -Iinclude -I../comostrings/include -I../camxml/include \
       -I../sprogc/include -I../cvode/include -I../radau/include -I../brush/include \
       -I../geometry/include  -I../mopsc/include \
       -I../sweepc/include -I../utils/include \
       -I../cvode/source -Wall
    flink= 
  else
    flag=-c -O2 -Iinclude -I../comostrings/include -I../camxml/include \
       -I../sprogc/include -I../cvode/include -I../radau/include -I../brush/include \
       -I../geometry/include  -I../mopsc/include \
       -I../sweepc/include -I../utils/include \
       -I../cvode/source -Wall -pg
    flink=-pg
  endif
else
  flag=-c -O0 -Iinclude -I../comostrings/include -I../camxml/include \
       -I../sprogc/include -I../cvode/include -I../radau/include -I../brush/include \
       -I../geometry/include  -I../mopsc/include \
       -I../sweepc/include -I../utils/include \
       -I../cvode/source -I../radau/source -Wall -ggdb -D_GLIBCXX_DEBUG
  flink=-pg 
endif

# Source directories.

src=source
inc=include

export path=.
lib=$(path)/lib

# Additional libraries.

ifeq ($(debug),0)
  comostrings=../comostrings/lib/comostrings.a
  camxml=../camxml/lib/camxml.a
  sprog=../sprogc/lib/sprog.a
  cvode=../cvode/lib/cvode.a
  radau=../radau/lib/radau.a
else
  comostrings=../comostrings/lib/comostrings_d.a
  camxml=../camxml/lib/camxml_d.a
  sprog=../sprogc/lib/sprog_d.a
  cvode=../cvode/lib/cvode_d.a
  radau=../radau/lib/radau_d.a
endif
linklibs=$(comostrings) $(camxml) $(sprog) $(cvode) $(radau) $(sweep) $(mops) $(brush)

# Target objects.

obj=   cam_read.o cam_setup.o cam_geometry.o \
		 cam_control.o cam_conveter.o cam_math.o cam_conc.o \
		 cam_admin.o cam_boundary.o cam_configuration.o cam_profile.o\
		 cam_reporter.o cam_models.o cam_plug.o cam_premix.o cam_residual.o \
		 stagflow.o cam_soot.o batch.o flamelet.o interface.o\
		 cvode_wrapper.o radau_wrapper.o

objlocal=	cam_kernel.o

# TARGET binary: Compile to a  binary.

binary: comostringslib camxmllib sproglib cvodelib radaulib $(obj) $(objlocal)
	mkdir -p bin
ifeq ($(debug),0)
  ifeq ($(profile),0)
		$(cpp) $(flink) $(obj) $(objlocal) $(linklibs) -o bin/camflow
  else
		$(cpp) $(flink) $(obj) $(objlocal) $(linklibs) -o bin/camflow_p
  endif 
else
	$(cpp) $(flink) $(obj) $(objlocal) $(linklibs) -o bin/camflow_d
endif

library: $(obj) soot_flamelet.o
	mkdir -p $(lib)
ifeq ($(debug),0)
	ar rc $(lib)/camflow.a $(obj) soot_flamelet.o $(linklibs)
else
	ar rc $(lib)/camflow_d.a $(obj) soot_flamelet.o $(linklibs) 
endif


# Required libraries.
comostringslib:
	$(MAKE) -f ../comostrings/Makefile -C ../comostrings debug=$(debug) \
  profile=$(profile)
camxmllib:
	$(MAKE) -f ../camxml/Makefile -C ../camxml debug=$(debug) profile=$(profile)
sproglib:
	$(MAKE) -f ../sprogc/Makefile -C ../sprogc debug=$(debug) profile=$(profile)
cvodelib:
	$(MAKE) -f ../cvode/Makefile -C ../cvode debug=$(debug) profile=$(profile)
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
batch.o: $(src)/batch.cpp
	$(cpp) $(flag) $(src)/batch.cpp
cam_soot.o: $(src)/cam_soot.cpp
	$(cpp) $(flag) $(src)/cam_soot.cpp
interface.o: $(src)/interface.cpp
	$(cpp) $(flag) $(src)/interface.cpp
cvode_wrapper.o: $(src)/cvode_wrapper.cpp
	$(cpp) $(flag) $(src)/cvode_wrapper.cpp
radau_wrapper.o: $(src)/radau_wrapper.cpp
	$(cpp) $(flag) $(src)/radau_wrapper.cpp
flamelet.o: $(src)/flamelet.cpp
	$(cpp) $(flag) $(src)/flamelet.cpp
soot_flamelet.o: $(src)/soot_flamelet.cpp
	$(cpp) $(flag) $(src)/soot_flamelet.cpp


# Combined mechanisms.


# TARGET clean: leave just source code.
clean:
	rm -rf *.o \#* *~

cleanlibs:
	$(MAKE) clean -f ../comostrings/Makefile -C ../comostrings
	$(MAKE) clean -f ../camxml/Makefile -C ../camxml
	$(MAKE) clean -f ../sprogc/Makefile -C ../sprogc
	$(MAKE) clean -f ../cvode/Makefile -C ../cvode
	rm -rf *.o \#* *~

