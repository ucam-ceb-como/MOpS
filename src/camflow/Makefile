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
includes=-Iinclude -I../io/comostrings/include -I../io/camxml/include \
       -I../sprogc/include -I../odesolvers/cvode/include -I../odesolvers/radau/include \
       -I../geometry/include -I../utils/include

ifeq ($(debug),0)
  ifeq ($(profile),0)
    flag=-c -O2 $(includes) -Wall
    flink= 
    lib=$(librelease)
  else
    flag=-c -O2 $(includes) -Wall -pg
    flink=-pg
    lib=$(libprofile)
  endif
else
  flag=-c -O0 $(includes) -Wall -ggdb -pg -D_GLIBCXX_DEBUG
  flink=-pg
  lib=$(libdebug)
endif

# Target objects.
obj=    cam_read.o cam_setup.o cam_geometry.o \
		cam_control.o cam_conveter.o cam_math.o cam_conc.o cam_radiation.o\
		cam_admin.o cam_boundary.o cam_configuration.o cam_profile.o\
		cam_reporter.o cam_models.o cam_residual.o \
		cam_soot.o cam_sdr.o flamelet.o batch.o \
		cvode_wrapper.o radau_wrapper.o

all: $(obj) 
	mkdir -p $(lib)
	ar rc $(lib)/libcamflow.a $(obj)

# Source code.

# Driver program.
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
cam_radiation.o: $(src)/cam_radiation.cpp
	$(cpp) $(flag) $(src)/cam_radiation.cpp
cam_models.o: $(src)/cam_models.cpp
	$(cpp) $(flag) $(src)/cam_models.cpp
cam_soot.o: $(src)/cam_soot.cpp
	$(cpp) $(flag) $(src)/cam_soot.cpp
cvode_wrapper.o: $(src)/cvode_wrapper.cpp
	$(cpp) $(flag) $(src)/cvode_wrapper.cpp
radau_wrapper.o: $(src)/radau_wrapper.cpp
	$(cpp) $(flag) $(src)/radau_wrapper.cpp
cam_sdr.o: $(src)/cam_sdr.cpp
	$(cpp) $(flag) $(src)/cam_sdr.cpp
flamelet.o: $(src)/flamelet.cpp
	$(cpp) $(flag) $(src)/flamelet.cpp
batch.o: $(src)/batch.cpp
	$(cpp) $(flag) $(src)/batch.cpp
	
# Combined mechanisms.

# TARGET clean: leave just source code.
clean:
	rm -rf *.o \#* *~

cleanlibs:
	rm -rf *.o \#* *~ $(librelease)/libcamflow.a $(libdebug)/libcamflow.a $(libprofile)/libcamflow.a
