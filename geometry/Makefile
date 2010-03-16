# Compiler options.
export debug=1
export profile=0

# C++ compiler.
cpp=g++
ifeq ($(debug),0)
  ifeq ($(profile),0)
    flag=-c -O2 -Iinclude -I../camxml/include \
       -I../sprogc/include -Wall
    flink=
  else
    flag=-c -O2 -Iinclude -I../camxml/include \
       -I../sprogc/include -Wall -pg
    flink=-pg
  endif
else
  flag=-c -O0 -Iinclude -I../camxml/include \
       -I../sprogc/include -Wall -ggdb -D_GLIBCXX_DEBUG
  flink=
endif


ifeq ($(debug),0)
  geomlib=./lib/geometry.a
else
  geomlib=./lib/geometry_d.a
endif
# Source directories.

export path=.
src=$(path)/source
inc=$(path)/include
lib=$(path)/lib

obj=geometry1d.o local_geometry1d.o


# Bring everything  into a library
$(geomlib): $(obj)
	mkdir -p lib
ifeq ($(debug),0)
	ar -rc lib/geometry.a $(obj)
else
	ar -rc lib/geometry_d.a $(obj)
endif

# Build the objects that contain the code for the real work

geometry1d.o: $(src)/geometry1d.cpp $(inc)/geometry1d.h $(inc)/types.h
	$(cpp) $(flag) $(src)/geometry1d.cpp

local_geometry1d.o: $(src)/local_geometry1d.cpp $(inc)/local_geometry1d.h $(inc)/geometry1d.h $(inc)/types.h
	$(cpp) $(flag) $(src)/local_geometry1d.cpp


# leave just source code.
clean:
	rm -rf $(path)/*.o \#* *~ $(lib)/*.a

rebuild:
	rm -rf $(path)/*.o \#* *~ $(lib)/*.a
	make geomlib

