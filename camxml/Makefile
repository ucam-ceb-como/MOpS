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

obj=stringf.o unicodef.o xmlattribute.o xmlelement.o xmldocument.o

# TARGET all: Compile to a static library (*.a).

all: $(obj)
	mkdir -p $(lib)
ifeq ($(debug),0)
	ar rc $(lib)/camxml.a $(obj)
else
	ar rc $(lib)/camxml_d.a $(obj)
endif

# Source code.
stringf.o: $(src)/stringf.cpp
	$(cpp) $(flag) $(src)/stringf.cpp

unicodef.o: $(src)/unicodef.cpp
	$(cpp) $(flag) $(src)/unicodef.cpp

xmlattribute.o: $(src)/xmlattribute.cpp unicodef.o
	$(cpp) $(flag) $(src)/xmlattribute.cpp

xmlelement.o: $(src)/xmlelement.cpp unicodef.o stringf.o xmlattribute.o
	$(cpp) $(flag) $(src)/xmlelement.cpp

xmldocument.o: $(src)/xmldocument.cpp unicodef.o xmlattribute.o xmlelement.o
	$(cpp) $(flag) $(src)/xmldocument.cpp

# TARGET clean: leave just source code.

clean:
	rm -rf $(path)/*.o \#* *~
