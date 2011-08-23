# Specification file for the CMake build system
#
# Robert I A Patterson
# 28 April 2011

project (COMOSTRINGS)

file(GLOB_RECURSE INCS "include/*.h")
source_group("Header Files" FILES ${INCS})

add_library(comostrings source/console_io.cpp
                        source/csv_io.cpp
                        source/data_io.cpp
                        source/string_functions.cpp
                        ${INCS})

include_directories (include)
