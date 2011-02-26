# Compiler options.
export debug=1
export profile=0

export path=.
librelease=../../../lib/release
libdebug=../../../lib/debug
libprofile=../../../lib/profile
binrelease=../../../bin/release
bindebug=../../../bin/debug
binprofile=../../../bin/profile
source=../../../src

# C++ compiler.
cpp=g++
#cpp=mpicxx -DUSE_MPI -DNONSPHERE

# Includes.
includes=-I$(source)/mopsc/include -I$(source)/io/comostrings/include -I$(source)/io/camxml/include \
         -I$(source)/sprogc/include -I$(source)/sweepc/include -I$(source)/odesolvers/cvodes/include \
         -I$(source)/odesolvers/cvodes/source/cvodes -I$(source)/geometry/include -I$(source)/utils/include

ifeq ($(debug),0)
  ifeq ($(profile),0)
    flag=-c -O2 $(includes) -Wall 
    flink= 
    lib=$(librelease)
    bin=$(binrelease)
  else
    flag=-c -O2 $(includes) -Wall -pg
    flink=-pg
    lib=$(libprofile)
    bin=$(binprofile)
  endif
else
  flag=-c -O0 $(includes) -Wall -ggdb -pg -D_GLIBCXX_DEBUG
  flink=-pg 
  lib=$(libdebug)
  bin=$(bindebug)
endif

# Additional libraries.
mops=$(lib)/libmops.a
comostrings=$(lib)/libcomostrings.a
camxml=$(lib)/libcamxml.a
sprog=$(lib)/libsprog.a
sweep=$(lib)/libsweep.a
cvodes=$(lib)/libcvodes.a
geom=$(lib)/libgeometry.a

linklibs= $(mops) $(comostrings) $(camxml) $(sprog) $(sweep) $(cvodes) $(geom)

# TARGET binary: Compile to a  binary.
objlocal=mops.o

all: $(objlocal)
	mkdir -p $(bin)
	$(cpp) $(flink) $(objlocal) $(linklibs) -o $(bin)/mops

32bit: $(objlocal)
	mkdir -p $(bin)
	$(cpp) $(flink) $(objlocal) $(linklibs) -o $(bin)/mops -static

# Driver program.
mops.o: mops.cpp
	$(cpp) $(flag) mops.cpp

# TARGET clean: leave just source code.
clean:
	rm -rf *.o \#* *~
