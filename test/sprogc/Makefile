# Compiler options.
export debug=1
export profile=0

export path=.
librelease=../../lib/release
libdebug=../../lib/debug
libprofile=../../lib/profile
binrelease=../../bin/release
bindebug=../../bin/debug
binprofile=../../bin/profile
source=../../src

# C++ compiler.
cpp=g++
# Includes.
includes=-I$(source)/io/comostrings/include \
         -I$(source)/utils/include \
         -I$(source)/sprogc/include

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

comostrings=$(lib)/libcomostrings.a
sprog=$(lib)/libsprog.a

linklibs=$(sprog) $(comostrings)

# TARGET binary: Compile to a binary.
objlocal=sprogc_test_main.o

all: $(objlocal)
	$(cpp) $(flink) $(objlocal) $(linklibs) -o sprogc_test_main


# Driver program.
sprogc_test_main.o: sprogc_test_main.cpp
	$(cpp) $(flag) sprogc_test_main.cpp

# TARGET clean: leave just source code.
clean:
	rm -rf *.o \#* *~
