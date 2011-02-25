# Compiler options.
export debug=1
export profile=0

export path=.
src=$(path)/source
inc=$(path)/include
librelease=../../lib/release
libdebug=../../lib/debug

# C++ compiler.
cpp=g++
# Includes.
includes=-Iinclude -I../io/camxml/include -I../sprogc/include 

ifeq ($(debug),0)
  ifeq ($(profile),0)
    flag=-c -O2 $(includes) -Wall
    flink=
    lib=$(librelease)
  else
    flag=-c -O2 $(includes) -Wall -pg
    flink=-pg
    lib=$(librelease)
  endif
else
  flag=-c -O0 $(includes) -Wall -ggdb -D_GLIBCXX_DEBUG
  flink=
  lib=$(libdebug)
endif

# Target objects.
obj=geometry1d.o local_geometry1d.o

geomlib=$(lib)/libgeometry.a

# Bring everything  into a library
all: $(obj)
	mkdir -p $(lib)
	ar -rc $(lib)/libgeometry.a $(obj)

# Build the objects that contain the code for the real work

geometry1d.o: $(src)/geometry1d.cpp $(inc)/geometry1d.h $(inc)/types.h
	$(cpp) $(flag) $(src)/geometry1d.cpp

local_geometry1d.o: $(src)/local_geometry1d.cpp $(inc)/local_geometry1d.h $(inc)/geometry1d.h $(inc)/types.h
	$(cpp) $(flag) $(src)/local_geometry1d.cpp


# leave just source code.
clean:
	rm -rf $(path)/*.o \#* *~

cleanlibs:
	rm -rf $(path)/*.o \#* *~ $(librelease)/libgeometry.a $(libdebug)/libgeometry.a

rebuild:
	rm -rf $(path)/*.o \#* *~ $(librelease)/libgeometry.a $(libdebug)/libgeometry.a
	make geomlib

