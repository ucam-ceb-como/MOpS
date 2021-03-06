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
ifeq ($(debug),0)
  ifeq ($(profile),0)
    flag=-c -O2 -Iinclude -Wall
    lib=$(librelease)
  else
    flag=-c -O2 -Iinclude -Wall -pg
    lib=$(libprofile)
  endif
else
  flag=-c -O0 -Iinclude -Wall -ggdb -pg
  lib=$(libdebug)
endif

# Target objects.
obj=ida.o ida_band.o ida_bbdpre.o ida_dense.o ida_ic.o \
    ida_io.o ida_spbcgs.o ida_spgmr.o ida_spils.o \
    ida_sptfqmr.o fnvector_serial.o nvector_serial.o sundials_band.o \
    sundials_dense.o sundials_iterative.o sundials_math.o sundials_nvector.o \
    sundials_smalldense.o sundials_spbcgs.o sundials_spgmr.o sundials_sptfqmr.o

# TARGET all: Compile to a static library (*.a).
all: $(obj)
	mkdir -p $(lib)
	ar rc $(lib)/libida.a $(obj)

# Source code.

# CVODE.

ida.o: $(src)/ida.c
	$(cpp) $(flag) $(src)/ida.c

ida_band.o: $(src)/ida_band.c
	$(cpp) $(flag) $(src)/ida_band.c

ida_bbdpre.o: $(src)/ida_bbdpre.c
	$(cpp) $(flag) $(src)/ida_bbdpre.c

ida_dense.o: $(src)/ida_dense.c
	$(cpp) $(flag) $(src)/ida_dense.c

ida_ic.o: $(src)/ida_ic.c
	$(cpp) $(flag) $(src)/ida_ic.c

ida_io.o: $(src)/ida_io.c
	$(cpp) $(flag) $(src)/ida_io.c

ida_spbcgs.o: $(src)/ida_spbcgs.c
	$(cpp) $(flag) $(src)/ida_spbcgs.c

ida_spgmr.o: $(src)/ida_spgmr.c
	$(cpp) $(flag) $(src)/ida_spgmr.c

ida_spils.o: $(src)/ida_spils.c
	$(cpp) $(flag) $(src)/ida_spils.c

ida_sptfqmr.o: $(src)/ida_sptfqmr.c
	$(cpp) $(flag) $(src)/ida_sptfqmr.c

# NVECTOR.

fnvector_serial.o: $(src)/nvec_ser/fnvector_serial.c
	$(cpp) $(flag) $(src)/nvec_ser/fnvector_serial.c

nvector_serial.o: $(src)/nvec_ser/nvector_serial.c
	$(cpp) $(flag) $(src)/nvec_ser/nvector_serial.c

# SUNDIALS.

sundials_band.o: $(src)/sundials/sundials_band.c
	$(cpp) $(flag) $(src)/sundials/sundials_band.c

sundials_dense.o: $(src)/sundials/sundials_dense.c
	$(cpp) $(flag) $(src)/sundials/sundials_dense.c

sundials_iterative.o: $(src)/sundials/sundials_iterative.c
	$(cpp) $(flag) $(src)/sundials/sundials_iterative.c

sundials_math.o: $(src)/sundials/sundials_math.c
	$(cpp) $(flag) $(src)/sundials/sundials_math.c

sundials_nvector.o: $(src)/sundials/sundials_nvector.c
	$(cpp) $(flag) $(src)/sundials/sundials_nvector.c

sundials_smalldense.o: $(src)/sundials/sundials_smalldense.c
	$(cpp) $(flag) $(src)/sundials/sundials_smalldense.c

sundials_spbcgs.o: $(src)/sundials/sundials_spbcgs.c
	$(cpp) $(flag) $(src)/sundials/sundials_spbcgs.c

sundials_spgmr.o: $(src)/sundials/sundials_spgmr.c
	$(cpp) $(flag) $(src)/sundials/sundials_spgmr.c

sundials_sptfqmr.o: $(src)/sundials/sundials_sptfqmr.c
	$(cpp) $(flag) $(src)/sundials/sundials_sptfqmr.c



# TARGET clean: leave just source code.

clean:
	rm -rf $(path)/*.o \#* *~

cleanlibs:
	rm -rf $(path)/*.o \#* *~ $(librelease)/libida.a $(libdebug)/libida.a $(libprofile)/libida.a
