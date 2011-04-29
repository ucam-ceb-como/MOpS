# Specification file for the CMake cross platform build system

project (CAMXML)

add_library(camxml source/stringf.cpp
                   source/unicodef.cpp
                   source/xmlattribute.cpp
                   source/xmldocument.cpp
                   source/xmlelement.cpp)

include_directories (include ../comostrings/include)

# Anyone linking against this library will also need to link against the comostrings
# library, because it is references from camxml.
target_link_libraries(camxml comostrings)