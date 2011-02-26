# Compiler options.
export debug=1
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
obj=cvode.o cvode_band.o cvode_bandpre.o cvode_bbdpre.o cvode_dense.o \
    cvode_diag.o cvode_io.o cvode_spbcgs.o cvode_spgmr.o cvode_spils.o \
    cvode_sptfqmr.o fnvector_serial.o nvector_serial.o sundials_band.o \
    sundials_dense.o sundials_iterative.o sundials_math.o sundials_nvector.o \
    sundials_smalldense.o sundials_spbcgs.o sundials_spgmr.o sundials_sptfqmr.o

# TARGET all: Compile to a static library (*.a).
all: $(obj)
	mkdir -p $(lib)
	ar rc $(lib)/libcvode.a $(obj)

# Source code.

# CVODE.

cvode.o: $(src)/cvode.c
	$(cpp) $(flag) $(src)/cvode.c

cvode_band.o: $(src)/cvode_band.c
	$(cpp) $(flag) $(src)/cvode_band.c

cvode_bandpre.o: $(src)/cvode_bandpre.c
	$(cpp) $(flag) $(src)/cvode_bandpre.c

cvode_bbdpre.o: $(src)/cvode_bbdpre.c
	$(cpp) $(flag) $(src)/cvode_bbdpre.c

cvode_dense.o: $(src)/cvode_dense.c
	$(cpp) $(flag) $(src)/cvode_dense.c

cvode_diag.o: $(src)/cvode_diag.c
	$(cpp) $(flag) $(src)/cvode_diag.c

cvode_io.o: $(src)/cvode_io.c
	$(cpp) $(flag) $(src)/cvode_io.c

cvode_spbcgs.o: $(src)/cvode_spbcgs.c
	$(cpp) $(flag) $(src)/cvode_spbcgs.c

cvode_spgmr.o: $(src)/cvode_spgmr.c
	$(cpp) $(flag) $(src)/cvode_spgmr.c

cvode_spils.o: $(src)/cvode_spils.c
	$(cpp) $(flag) $(src)/cvode_spils.c

cvode_sptfqmr.o: $(src)/cvode_sptfqmr.c
	$(cpp) $(flag) $(src)/cvode_sptfqmr.c

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
	rm -rf $(path)/*.o \#* *~ $(librelease)/libcvode.a $(libdebug)/libcvode.a $(libprofile)/libcvode.a
