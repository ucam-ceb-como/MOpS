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
  flag=-c -O0 -Iinclude -Wall -ggdb -pg
endif

# Source directories.

export path=.
src=$(path)/source
inc=$(path)/include
lib=$(path)/lib

# Target objects.

obj=cvodea.o cvodea_io.o cvodes.o cvodes_band.o cvodes_bandpre.o \
    cvodes_bbdpre.o cvodes_dense.o cvodes_diag.o cvodes_direct.o cvodes_io.o cvodes_lapack.io \
    cvodes_spbcgs.o cvodes_spgmr.o cvodes_spils.o cvodes_sptfqmr.o \
    fnvector_serial.o nvector_serial.o sundials_band.o sundials_dense.o \
    sundials_iterative.o sundials_math.o sundials_nvector.o sundials_direct.o \
    sundials_spbcgs.o sundials_spgmr.o sundials_sptfqmr.o


# TARGET all: Compile to a static library (*.a).

all: $(obj)
	mkdir -p $(lib)
ifeq ($(debug),0)
	ar rc $(lib)/cvodes.a $(obj)
else
	ar rc $(lib)/cvodes_d.a $(obj)
endif

# Source code.

# CVODE.

cvodes.o: $(src)/cvodes/cvodes.c
	$(cpp) $(flag) $(src)/cvodes/cvodes.c

cvodea.o: $(src)/cvodes/cvodea.c
	$(cpp) $(flag) $(src)/cvodes/cvodea.c

cvodea_io.o: $(src)/cvodes/cvodea_io.c
	$(cpp) $(flag) $(src)/cvodes/cvodea_io.c

cvodes_band.o: $(src)/cvodes/cvodes_band.c
	$(cpp) $(flag) $(src)/cvodes/cvodes_band.c

cvodes_bandpre.o: $(src)/cvodes/cvodes_bandpre.c
	$(cpp) $(flag) $(src)/cvodes/cvodes_bandpre.c

cvodes_bbdpre.o: $(src)/cvodes/cvodes_bbdpre.c
	$(cpp) $(flag) $(src)/cvodes/cvodes_bbdpre.c

cvodes_dense.o: $(src)/cvodes/cvodes_dense.c
	$(cpp) $(flag) $(src)/cvodes/cvodes_dense.c

cvodes_diag.o: $(src)/cvodes/cvodes_diag.c
	$(cpp) $(flag) $(src)/cvodes/cvodes_diag.c

cvodes_diag.o: $(src)/cvodes/cvodes_direct.c
	$(cpp) $(flag) $(src)/cvodes/cvodes_direct.c

cvodes_io.o: $(src)/cvodes/cvodes_io.c
	$(cpp) $(flag) $(src)/cvodes/cvodes_io.c

cvodes_spbcgs.o: $(src)/cvodes/cvodes_lapack.c
	$(cpp) $(flag) $(src)/cvodes/cvodes_lapack.c

cvodes_spbcgs.o: $(src)/cvodes/cvodes_spbcgs.c
	$(cpp) $(flag) $(src)/cvodes/cvodes_spbcgs.c

cvodes_spgmr.o: $(src)/cvodes/cvodes_spgmr.c
	$(cpp) $(flag) $(src)/cvodes/cvodes_spgmr.c

cvodes_spils.o: $(src)/cvodes/cvodes_spils.c
	$(cpp) $(flag) $(src)/cvodes/cvodes_spils.c

cvodes_sptfqmr.o: $(src)/cvodes/cvodes_sptfqmr.c
	$(cpp) $(flag) $(src)/cvodes/cvodes_sptfqmr.c

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

sundials_smalldense.o: $(src)/sundials/sundials_direct.c
	$(cpp) $(flag) $(src)/sundials/sundials_direct.c

sundials_nvector.o: $(src)/sundials/sundials_nvector.c
	$(cpp) $(flag) $(src)/sundials/sundials_nvector.c

sundials_spbcgs.o: $(src)/sundials/sundials_spbcgs.c
	$(cpp) $(flag) $(src)/sundials/sundials_spbcgs.c

sundials_spgmr.o: $(src)/sundials/sundials_spgmr.c
	$(cpp) $(flag) $(src)/sundials/sundials_spgmr.c

sundials_sptfqmr.o: $(src)/sundials/sundials_sptfqmr.c
	$(cpp) $(flag) $(src)/sundials/sundials_sptfqmr.c


# TARGET clean: leave just source code.

clean:
	rm -rf $(path)/*.o \#* *~
