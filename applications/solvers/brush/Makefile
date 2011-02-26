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
includes=-I$(source)/brush/include -I$(source)/io/comostrings/include -I$(source)/io/camxml/include \
         -I$(source)/geometry/include -I$(source)/utils/include -I$(source)/sprogc/include \
         -I$(source)/sweepc/include -I$(source)/mopsc/include 

ifeq ($(debug),0)
  ifeq ($(profile),0)
    flag=-c -O2 $(includes) -Wall
    flink=
    lib=$(librelease)
    bin=$(binrelease)
  else
    flag=-c -O2 (includes) -Wall -pg
    flink=-pg
    lib=$(librelease)
    bin=$(binprofile)
  endif
else
  flag=-c -O0 (includes) -Wall -ggdb -D_GLIBCXX_DEBUG
  flink=
  lib=$(libdebug)
  bin=$(bindebug)
endif

# Additional libraries.
comostrings=$(lib)/libcomostrings.a
camxml=$(lib)/libcamxml.a
sprog=$(lib)/libsprog.a
sweep=$(lib)/libsweep.a
mops=$(lib)/libmops.a
cvodes=$(lib)/libcvodes.a
geom=$(lib)/libgeometry.a
brushlib=$(lib)/libbrush.a

linklibs=$(mops) $(sweep) $(geom) $(sprog) $(camxml) $(comostrings) $(cvodes)

# Default target - executable
all: brush.o
	mkdir -p $(bin)
	$(cpp) $(flink) brush.o $(brushlib) $(linklibs) -o $(bin)/brush

# Build the objects that contain the code for the real work
brush.o: brush.cpp
	$(cpp) $(flag) brush.cpp

# leave just source code.
clean:
	rm -rf $(path)/*.o \#* *~
