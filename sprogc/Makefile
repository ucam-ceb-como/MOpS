# Compiler options.
export debug=0
export profile=0

# C++ compiler.
cpp=g++
ifeq ($(debug),0)
  ifeq ($(profile),0)
    flag=-c -O2 -Iinclude -I../comostrings/include -Wall
  else
    flag=-c -O2 -Iinclude -I../comostrings/include -Wall -pg
  endif
else
  flag=-c -O0 -Iinclude -I../comostrings/include -Wall -ggdb -pg
endif

# Source directories.

export path=.
src=$(path)/source
inc=$(path)/include
lib=$(path)/lib

# Additional libraries.

ifeq ($(debug),0)
  comostrings=../comostrings/lib/comostrings.a
else
  comostrings=../comostrings/lib/comostrings_d.a
endif

# Target objects.

obj=gpc_el_comp.o gpc_element.o gpc_element_lib.o gpc_gasphase.o \
    gpc_idealgas.o gpc_mech.o gpc_mech_io.o gpc_mixture.o \
    gpc_mixture_factory.o gpc_reaction.o gpc_reaction_set.o gpc_species.o \
    gpc_stoich.o gpc_thermo.o gpc_string.o gpc_transport_data.o gpc_transport_factory.o

# TARGET all: Compile to a static library (*.a).

all: $(obj)
	mkdir -p $(lib)
ifeq ($(debug),0)
	ar rc $(lib)/sprog.a $(obj) $(comostrings)
else
	ar rc $(lib)/sprog_d.a $(obj) $(comostrings)
endif

# Source code.

# Chemical elements.

gpc_el_comp.o: $(src)/gpc_el_comp.cpp gpc_element.o
	$(cpp) $(flag) $(src)/gpc_el_comp.cpp

gpc_element.o: $(src)/gpc_element.cpp
	$(cpp) $(flag) $(src)/gpc_element.cpp

gpc_element_lib.o: $(inc)/gpc_element_lib.cpp gpc_element.o
	$(cpp) $(flag) $(inc)/gpc_element_lib.cpp

# Chemical species.

gpc_species.o: $(src)/gpc_species.cpp gpc_el_comp.o
	$(cpp) $(flag) $(src)/gpc_species.cpp

gpc_stoich.o: $(src)/gpc_stoich.cpp
	$(cpp) $(flag) $(src)/gpc_stoich.cpp

gpc_thermo.o: $(src)/gpc_thermo.cpp
	$(cpp) $(flag) $(src)/gpc_thermo.cpp

# Mixtures (incl. gases).

gpc_mixture.o: $(src)/gpc_mixture.cpp
	$(cpp) $(flag) $(src)/gpc_mixture.cpp

gpc_gasphase.o: $(src)/gpc_gasphase.cpp 
	$(cpp) $(flag) $(src)/gpc_gasphase.cpp

gpc_idealgas.o: $(src)/gpc_idealgas.cpp
	$(cpp) $(flag) $(src)/gpc_idealgas.cpp
	
gpc_mixture_factory.o: $(src)/gpc_mixture_factory.cpp gpc_mixture.o \
                       gpc_gasphase.o gpc_idealgas.o
	$(cpp) $(flag) $(src)/gpc_mixture_factory.cpp

# Chemical reactions.

gpc_reaction.o: $(src)/gpc_reaction.cpp gpc_stoich.o
	$(cpp) $(flag) $(src)/gpc_reaction.cpp

gpc_reaction_set.o: $(src)/gpc_reaction_set.cpp gpc_reaction.o gpc_stoich.o
	$(cpp) $(flag) $(src)/gpc_reaction_set.cpp

# Mechanisms/Models.

gpc_mech.o: $(src)/gpc_mech.cpp
	$(cpp) $(flag) $(src)/gpc_mech.cpp

gpc_mech_io.o: $(src)/gpc_mech_io.cpp gpc_mech.o
	$(cpp) $(flag) $(src)/gpc_mech_io.cpp

gpc_string.o: $(src)/gpc_string.cpp
	$(cpp) $(flag) $(src)/gpc_string.cpp

# transport models
gpc_transport_data.o: $(src)/gpc_transport_data.cpp
	$(cpp) $(flag) $(src)/gpc_transport_data.cpp

gpc_transport_factory.o: $(src)/gpc_transport_factory.cpp
	$(cpp) $(flag) $(src)/gpc_transport_factory.cpp

# TARGET clean: leave just source code.

clean:
	rm -rf $(path)/*.o \#* *~
