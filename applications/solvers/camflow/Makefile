# Compiler options.
export debug=1
export profile=0

export path=.
librelease=../../../lib/release
libdebug=../../../lib/debug
binrelease=../../../bin/release
bindebug=../../../bin/debug
binprofile=../../../bin/profile
source=../../../src

# C++ compiler.
cpp=g++
# Includes.
includes=-I$(source)/io/comostrings/include -I$(source)/io/camxml/include \
       -I$(source)/sprogc/include -I$(source)/odesolvers/cvode/include -I$(source)/odesolvers/radau/include \
       -I$(source)/geometry/include -I$(source)/utils/include -I$(source)/camflow/include 

ifeq ($(debug),0)
  ifeq ($(profile),0)
    flag=-c -O2 $(includes) -Wall
    flink= 
    lib=$(librelease)
    bin=$(binrelease)
  else
    flag=-c -O2 $(includes) -Wall -pg
    flink=-pg
    lib=$(librelease)
    bin=$(binprofile)
  endif
else
  flag=-c -O0 $(includes) -Wall -ggdb -pg -D_GLIBCXX_DEBUG
  flink=-pg
  lib=$(libdebug)
  bin=$(bindebug)
endif

# Additional libraries.

camflow=$(lib)/libcamflow.a
comostrings=$(lib)/libcomostrings.a
camxml=$(lib)/libcamxml.a
sprog=$(lib)/libsprog.a
cvode=$(lib)/libcvode.a
radau=$(lib)/libradau.a

linklibs=$(camflow) $(comostrings) $(camxml) $(sprog) $(cvode) $(radau)

# TARGET binary: Compile to a  binary.
objlocal=cam_kernel.o

all: $(objlocal)
	mkdir -p $(bin)
	$(cpp) $(flink) $(objlocal) $(linklibs) -o $(bin)/camflow


# Driver program.

cam_kernel.o: cam_kernel.cpp
	$(cpp) $(flag) cam_kernel.cpp

# TARGET clean: leave just source code.
clean:
	rm -rf *.o \#* *~
