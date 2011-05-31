# Specification file for the CMake build system

#
# Robert I A Patterson
# 29 April 2011
cmake_minimum_required(VERSION 2.6)

project (CVODE)

add_library(cvode source/cvode.c
                  source/cvode_band.c
                  source/cvode_bandpre.c
                  source/cvode_bbdpre.c
                  source/cvode_dense.c
                  source/cvode_diag.c
                  source/cvode_io.c
                  source/cvode_spbcgs.c
                  source/cvode_spgmr.c
                  source/cvode_spils.c
                  source/cvode_sptfqmr.c
                  source/nvec_ser/fnvector_serial.c
                  source/nvec_ser/nvector_serial.c
                  source/sundials/sundials_band.c
                  source/sundials/sundials_dense.c
                  source/sundials/sundials_iterative.c
                  source/sundials/sundials_math.c
                  source/sundials/sundials_nvector.c
                  source/sundials/sundials_smalldense.c
                  source/sundials/sundials_spbcgs.c
                  source/sundials/sundials_spgmr.c
                  source/sundials/sundials_sptfqmr.c)

include_directories(include)
