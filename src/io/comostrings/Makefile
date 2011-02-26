# Compiler options.
export debug=0
export profile=0

export path=.
src=$(path)/source
inc=$(path)/include
librelease=../../../lib/release
libdebug=../../../lib/debug
libprofile=../../../lib/profile

# C++ compiler.
cpp=g++
# Includes.
includes=-Iinclude

ifeq ($(debug),0)
  ifeq ($(profile),0)
    flag=-c -O2 $(includes) -Wall
    lib=$(librelease)
  else
    flag=-c -O2 $(includes) -Wall -pg
    lib=$(libprofile)
  endif
else
  flag=-c -O0 $(includes) -Wall -ggdb -pg -D_GLIBCXX_DEBUG
  lib=$(libdebug)
endif

# Target objects.
obj=string_functions.o csv_io.o console_io.o data_io.o

# TARGET all: Compile to a static library (*.a).
all: $(obj)
	mkdir -p $(lib)
	ar rc $(lib)/libcomostrings.a $(obj)

# Source code.

string_functions.o: $(src)/string_functions.cpp
	$(cpp) $(flag) $(src)/string_functions.cpp

csv_io.o: $(src)/csv_io.cpp string_functions.o
	$(cpp) $(flag) $(src)/csv_io.cpp

console_io.o: $(src)/console_io.cpp
	$(cpp) $(flag) $(src)/console_io.cpp

data_io.o: $(src)/data_io.cpp
	$(cpp) $(flag) $(src)/data_io.cpp
# TARGET clean: leave just source code.

clean:
	rm -rf $(path)/*.o \#* *~

cleanlibs:
	rm -rf $(path)/*.o \#* *~ $(librelease)/libcomostrings.a $(libdebug)/libcomostrings.a $(libprofile)/libcomostrings.a
