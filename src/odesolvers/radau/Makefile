# Compiler options.
export debug=1
export profile=0

# C++ compiler.
cpp=g++
ifeq ($(debug),0)
  ifeq ($(profile),0)
    flag=-c -O2 -Iinclude -Wall
  else
    flag=-c -O2 -Iinclude -Wall -pg
  endif
else
  flag=-c -O0 -Iinclude -Wall -ggdb -pg
endif

# Source directories.

export path=.
src=$(path)/source
inc=$(path)/include
lib=$(path)/lib
# Target objects.

obj=decsol.o StiffIntegratorT.o IntegratorT.o

# TARGET all: Compile to a static library (*.a).

all: $(obj)
	mkdir -p $(lib)
ifeq ($(debug),0)
	ar rc $(lib)/radau.a $(obj)
else
	ar rc $(lib)/radau_d.a $(obj)
endif

# Source code.

# RADAU

decsol.o: $(src)/decsol.cpp
	$(cpp) $(flag) $(src)/decsol.cpp

IntegratorT.o: $(src)/IntegratorT.cpp
	$(cpp) $(flag) $(src)/IntegratorT.cpp

StiffIntegratorT.o: $(src)/StiffIntegratorT.cpp
	$(cpp) $(flag) $(src)/StiffIntegratorT.cpp


# TARGET clean: leave just source code.

clean:
	rm -rf $(path)/*.o
