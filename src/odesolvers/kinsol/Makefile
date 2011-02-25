# Compiler options.
export debug=0
export profile=0

export path=.
src=$(path)/source
inc=$(path)/include
librelease=../../../lib/release
libdebug=../../../lib/debug

# C++ compiler.
cpp=g++
ifeq ($(debug),0)
  ifeq ($(profile),0)
    flag=-c -O2 -Iinclude -Wall
    lib=$(librelease)
  else
    flag=-c -O2 -Iinclude -Wall -pg
    lib=$(librelease)
  endif
else
  flag=-c -O0 -Iinclude -Wall -ggdb -pg
  lib=$(libdebug)
endif

# Target objects.
obj=kinsol.o kinsol_band.o kinsol_bbdpre.o kinsol_dense.o \
    kinsol_io.o kinsol_spbcgs.o kinsol_spgmr.o kinsol_spils.o \
    kinsol_sptfqmr.o fnvector_serial.o nvector_serial.o sundials_band.o \
    sundials_dense.o sundials_iterative.o sundials_math.o sundials_nvector.o \
    sundials_smalldense.o sundials_spbcgs.o sundials_spgmr.o sundials_sptfqmr.o

# TARGET all: Compile to a static library (*.a).
all: $(obj)
	mkdir -p $(lib)
	ar rc $(lib)/libkinsol.a $(obj)

# Source code.

# CVODE.

kinsol.o: $(src)/kinsol.c
	$(cpp) $(flag) $(src)/kinsol.c

kinsol_band.o: $(src)/kinsol_band.c
	$(cpp) $(flag) $(src)/kinsol_band.c

kinsol_bandpre.o: $(src)/kinsol_bandpre.c
	$(cpp) $(flag) $(src)/kinsol_bandpre.c

kinsol_bbdpre.o: $(src)/kinsol_bbdpre.c
	$(cpp) $(flag) $(src)/kinsol_bbdpre.c

kinsol_dense.o: $(src)/kinsol_dense.c
	$(cpp) $(flag) $(src)/kinsol_dense.c

kinsol_diag.o: $(src)/kinsol_diag.c
	$(cpp) $(flag) $(src)/kinsol_diag.c

kinsol_io.o: $(src)/kinsol_io.c
	$(cpp) $(flag) $(src)/kinsol_io.c

kinsol_spbcgs.o: $(src)/kinsol_spbcgs.c
	$(cpp) $(flag) $(src)/kinsol_spbcgs.c

kinsol_spgmr.o: $(src)/kinsol_spgmr.c
	$(cpp) $(flag) $(src)/kinsol_spgmr.c

kinsol_spils.o: $(src)/kinsol_spils.c
	$(cpp) $(flag) $(src)/kinsol_spils.c

kinsol_sptfqmr.o: $(src)/kinsol_sptfqmr.c
	$(cpp) $(flag) $(src)/kinsol_sptfqmr.c

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
	rm -rf $(path)/*.o \#* *~ $(librelease)/libkinsol.a $(libdebug)/libkinsol.a
