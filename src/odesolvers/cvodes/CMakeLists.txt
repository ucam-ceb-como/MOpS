# Specification file for the CMake build system

#
# Robert I A Patterson
# 29 April 2011
cmake_minimum_required(VERSION 2.6)

project (CVODES)

add_library(cvodes source/cvodes/cvodea.c
                   source/cvodes/cvodea_io.c
                   source/cvodes/cvodes_band.c
                   source/cvodes/cvodes_bandpre.c
                   source/cvodes/cvodes_bbdpre.c
                   source/cvodes/cvodes.c
                   source/cvodes/cvodes_dense.c
                   source/cvodes/cvodes_diag.c
                   source/cvodes/cvodes_direct.c
                   source/cvodes/cvodes_io.c
                   source/cvodes/cvodes_lapack.c
                   source/cvodes/cvodes_spbcgs.c
                   source/cvodes/cvodes_spgmr.c
                   source/cvodes/cvodes_spils.c
                   source/cvodes/cvodes_sptfqmr.c
                   source/nvec_ser/fnvector_serial.c
                   source/nvec_ser/nvector_serial.c
                   source/sundials/sundials_band.c
                   source/sundials/sundials_dense.c
                   source/sundials/sundials_direct.c
                   source/sundials/sundials_iterative.c
                   source/sundials/sundials_math.c
                   source/sundials/sundials_nvector.c
                   source/sundials/sundials_spbcgs.c
                   source/sundials/sundials_spgmr.c
                   source/sundials/sundials_sptfqmr.c)

include_directories(include)
