# Compiler options.
export debug=0
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
  flag=-c -O0 -Iinclude -Wall -ggdb -pg -D_GLIBCXX_DEBUG
endif

# Source directories.

export path=.
src=$(path)/source
inc=$(path)include
lib=$(path)/lib
# Target objects.

obj=string_functions.o csv_io.o console_io.o data_io.o

# TARGET all: Compile to a static library (*.a).

all: $(obj)
	mkdir -p $(lib)
ifeq ($(debug),0)
	ar rc $(lib)/comostrings.a $(obj)
else
	ar rc $(lib)/comostrings_d.a $(obj)
endif

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
	rm -rf $(path)/*.o  
